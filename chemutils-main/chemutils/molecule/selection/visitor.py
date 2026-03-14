# type: ignore

import re
from collections.abc import Callable
from typing import Any

from openeye import oechem

from chemutils.molecule.graph import iterate_substructure_matches

from .generated.SelectionLanguageParser import SelectionLanguageParser
from .generated.SelectionLanguageVisitor import SelectionLanguageVisitor

# Maximum error in a float for equality
# In future, replace with something that respects significant figures
MAX_DIFFERENCE = 1e-5


class SelectionVisitor(SelectionLanguageVisitor):
    def __init__(self, mol: oechem.OEMolBase):
        self.mol = mol
        self.atoms = list(mol.GetAtoms())
        self.atom_idx_to_atom: dict[int, oechem.OEAtomBase] = {
            atom.GetIdx(): atom for atom in self.atoms
        }
        self.atom_idx_to_chain_idx = {}
        self.atom_idx_to_residue_idx = {}

        _, self.atom_idx_to_molecule_idx = oechem.OEDetermineComponents(mol)
        _, self.atom_idx_to_ring_system_idx = oechem.OEDetermineRingSystems(mol)

        view = oechem.OEHierView(self.mol)
        res_idx = 0
        for chain_idx, chain in enumerate(view.GetChains()):
            for frag in chain.GetFragments():
                for res in frag.GetResidues():
                    for atom in res.GetAtoms():
                        idx = atom.GetIdx()
                        self.atom_idx_to_chain_idx[idx] = chain_idx
                        self.atom_idx_to_residue_idx[idx] = res_idx
                    res_idx += 1

    def visitFull_expr(self, ctx: SelectionLanguageParser.Full_exprContext) -> oechem.OEAtomBondSet:
        indices = super().visit(ctx.expr())
        atombondset = oechem.OEAtomBondSet()
        for idx in indices:
            atombondset.AddAtom(self.atom_idx_to_atom[idx])
        return atombondset

    def visitParenthesized_expr(
        self, ctx: SelectionLanguageParser.Parenthesized_exprContext
    ) -> set[int]:
        return self.visit(ctx.expr())

    def visitNot(self, ctx: SelectionLanguageParser.NotContext) -> set[int]:
        """Handle the `not` query, selecting atoms that do not satisfy the subexpression."""
        return set(self.atom_idx_to_atom) - self.visit(ctx.expr())

    def visitSame(self, ctx: SelectionLanguageParser.SameContext) -> set[int]:
        """Handle the `same chain as` query, expanding the selection to every atom in the same chain."""
        atom_property = self.visit(ctx.atom_property())
        sel: set[int] = self.visit(ctx.expr())
        values = {prop for idx in sel if (prop := atom_property(idx)) is not None}
        return {idx for idx in self.atom_idx_to_chain_idx if atom_property(idx) in values}

    def visitAtom_property(
        self, ctx: SelectionLanguageParser.Atom_propertyContext
    ) -> Callable[[int], Any]:
        match ctx.getText():
            case "chain":
                return lambda idx: self.atom_idx_to_chain_idx[idx]
            case "residue":
                return lambda idx: self.atom_idx_to_residue_idx[idx]
            case "molecule":
                return lambda idx: self.atom_idx_to_molecule_idx[idx]
            case "ring_system":
                return (
                    lambda idx: ring_id
                    if (ring_id := self.atom_idx_to_ring_system_idx[idx]) > 0
                    else None
                )
            case "name":
                return lambda idx: self.atom_idx_to_atom[idx].GetName().strip()
            case "element":
                return lambda idx: self.atom_idx_to_atom[idx].GetAtomicNum()
            case "formal_charge":
                return lambda idx: self.atom_idx_to_atom[idx].GetFormalCharge()
            case "degree":
                return lambda idx: self.atom_idx_to_atom[idx].GetDegree()
            case "hvy_degree":
                return lambda idx: self.atom_idx_to_atom[idx].GetHvyDegree()
            case "valence":
                return lambda idx: self.atom_idx_to_atom[idx].GetValence()
            case "hcount":
                return lambda idx: self.atom_idx_to_atom[idx].GetTotalHCount()
            case "resname":
                return (
                    lambda idx: oechem.OEAtomGetResidue(self.atom_idx_to_atom[idx])
                    .GetName()
                    .strip()
                )
            case "resnum":
                return lambda idx: oechem.OEAtomGetResidue(
                    self.atom_idx_to_atom[idx]
                ).GetResidueNumber()
            case "altloc":
                return (
                    lambda idx: oechem.OEAtomGetResidue(self.atom_idx_to_atom[idx])
                    .GetAlternateLocation()
                    .strip()
                )
            case "b_factor":
                return lambda idx: oechem.OEAtomGetResidue(self.atom_idx_to_atom[idx]).GetBFactor()
            case "occ":
                return lambda idx: oechem.OEAtomGetResidue(
                    self.atom_idx_to_atom[idx]
                ).GetOccupancy()
            case _:
                raise ValueError(f"Failed to handle atom property {ctx.getText()}")

    def visitBonded_to(self, ctx):
        sel: set[int] = self.visit(ctx.expr())
        new_sel = set()
        for idx in sel:
            atom = self.atom_idx_to_atom[idx]
            for nbr in atom.GetAtoms():
                new_sel.add(nbr.GetIdx())
        return new_sel - sel

    def visitWithin(self, ctx: SelectionLanguageParser.WithinContext) -> set[int]:
        if self.mol.GetDimension() != 3:
            raise ValueError("Molecule must be 3D to use `within` query.")
        radius = self.visit(ctx.float_())
        target = self.visit(ctx.expr())
        new_sel = set()
        for target_idx in target:
            for nbr in oechem.OEGetNearestNbrs(self.mol, self.atom_idx_to_atom[target_idx], radius):
                new_sel.add(nbr.GetBgn().GetIdx())
                new_sel.add(nbr.GetEnd().GetIdx())
        return new_sel

    def visitAround(self, ctx: SelectionLanguageParser.AroundContext) -> set[int]:
        if self.mol.GetDimension() != 3:
            raise ValueError("Molecule must be 3D to use `within` query.")
        radius = self.visit(ctx.float_())
        target = self.visit(ctx.expr())
        new_sel = set()
        for target_idx in target:
            for nbr in oechem.OEGetNearestNbrs(self.mol, self.atom_idx_to_atom[target_idx], radius):
                new_sel.add(nbr.GetBgn().GetIdx())
                new_sel.add(nbr.GetEnd().GetIdx())
        return new_sel - target

    def visitBeyond(self, ctx: SelectionLanguageParser.BeyondContext) -> set[int]:
        """Handle the `beyond dist of expr2` query, selecting atoms in `expr1` that are beyond `dist`` of an atom of `expr2`."""
        if self.mol.GetDimension() != 3:
            raise ValueError("Molecule must be 3D to use `beyond` query.")
        radius = self.visit(ctx.float_())
        sel: set[int] = set(self.atom_idx_to_atom)
        target: set[int] = self.visit(ctx.expr())
        for target_idx in target:
            for nbr in oechem.OEGetNearestNbrs(self.mol, self.atom_idx_to_atom[target_idx], radius):
                if (bgn_idx := nbr.GetBgn().GetIdx()) in sel:
                    sel.remove(bgn_idx)
                if (end_idx := nbr.GetEnd().GetIdx()) in sel:
                    sel.remove(end_idx)
        return sel

    def visitAnd(self, ctx: SelectionLanguageParser.AndContext) -> set[int]:
        """Handle the `and` query, selecting atoms that satisfy all of the subexpressions."""
        constituent_selections = [self.visit(expr) for expr in ctx.expr()]
        return set.intersection(*constituent_selections)

    def visitXor(self, ctx: SelectionLanguageParser.XorContext) -> set[int]:
        """Handler the `xor` query, selecting atoms that satisfy exclusively one of the subexpressions."""
        constituent_selections = [self.visit(expr) for expr in ctx.expr()]
        return set.symmetric_difference(*constituent_selections)

    def visitOr(self, ctx: SelectionLanguageParser.OrContext) -> set[int]:
        """Handler the `or` query, selecting atoms that satisfy any of the subexpressions."""
        constituent_selections = [self.visit(expr) for expr in ctx.expr()]
        return set.union(*constituent_selections)

    def visitAll(self, ctx: SelectionLanguageParser.AllContext) -> set[int]:
        """Handle the `all` keyword, selecting all atoms."""
        return set(self.atom_idx_to_atom)

    def visitNone(self, ctx: SelectionLanguageParser.NoneContext) -> set[int]:
        """Handle the `none` keyword, selecting no atoms."""
        return set()

    def visitChain(self, ctx: SelectionLanguageParser.ChainContext) -> set[int]:
        return self.evaluate_string_atom_property(
            ctx, lambda atom: oechem.OEAtomGetResidue(atom).GetChainID()
        )

    def visitResidue_name(self, ctx: SelectionLanguageParser.Residue_nameContext) -> set[int]:
        return self.evaluate_string_atom_property(
            ctx, lambda atom: oechem.OEAtomGetResidue(atom).GetName()
        )

    def visitResidue_number(self, ctx: SelectionLanguageParser.Residue_nameContext) -> set[int]:
        return self.evaluate_numeric_atom_property(
            ctx, lambda atom: oechem.OEAtomGetResidue(atom).GetResidueNumber()
        )

    def visitAlternate_location(
        self, ctx: SelectionLanguageParser.Alternate_locationContext
    ) -> set[int]:
        return self.evaluate_string_atom_property(
            ctx, lambda atom: oechem.OEAtomGetResidue(atom).GetAlternateLocation().strip()
        )

    def visitAtom_serial(self, ctx: SelectionLanguageParser.Atom_serialContext) -> set[int]:
        return self.evaluate_numeric_atom_property(
            ctx, lambda atom: oechem.OEAtomGetResidue(atom).GetSerialNumber()
        )

    def visitB_factor(self, ctx: SelectionLanguageParser.B_factorContext) -> set[int]:
        return self.evaluate_numeric_atom_property(
            ctx, lambda atom: oechem.OEAtomGetResidue(atom).GetBFactor()
        )

    def visitOccupancy(self, ctx: SelectionLanguageParser.OccupancyContext) -> set[int]:
        return self.evaluate_numeric_atom_property(
            ctx, lambda atom: oechem.OEAtomGetResidue(atom).GetOccupancy()
        )

    def visitHetatm(self, ctx: SelectionLanguageParser.HetatmContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if oechem.OEIsHetAtom()(atom)}

    def visitIndex(self, ctx: SelectionLanguageParser.IndexContext) -> set[int]:
        return self.evaluate_numeric_atom_property(ctx, lambda atom: atom.GetIdx())

    def visitElement(self, ctx: SelectionLanguageParser.ElementContext) -> set[int]:
        elements = [identifier.getText() for identifier in ctx.IDENTIFIER()]
        atomic_nums = {oechem.OEGetAtomicNum(element) for element in elements}
        return {atom.GetIdx() for atom in self.atoms if atom.GetAtomicNum() in atomic_nums}

    def visitAtom_name(self, ctx: SelectionLanguageParser.Atom_nameContext) -> set[int]:
        return self.evaluate_string_atom_property(ctx, lambda atom: atom.GetName())

    def visitFormal_charge(self, ctx: SelectionLanguageParser.Formal_chargeContext) -> set[int]:
        return self.evaluate_numeric_atom_property(ctx, lambda atom: atom.GetFormalCharge())

    def visitDegree(self, ctx: SelectionLanguageParser.DegreeContext) -> set[int]:
        return self.evaluate_numeric_atom_property(ctx, lambda atom: atom.GetDegree())

    def visitHeavy_degree(self, ctx: SelectionLanguageParser.Heavy_degreeContext) -> set[int]:
        return self.evaluate_numeric_atom_property(ctx, lambda atom: atom.GetHvyDegree())

    def visitValence(self, ctx: SelectionLanguageParser.ValenceContext) -> set[int]:
        return self.evaluate_numeric_atom_property(ctx, lambda atom: atom.GetValence())

    def visitHydrogen_count(self, ctx: SelectionLanguageParser.Hydrogen_countContext) -> set[int]:
        return self.evaluate_numeric_atom_property(ctx, lambda atom: atom.GetTotalHCount())

    def visitRing_size(self, ctx: SelectionLanguageParser.Ring_sizeContext) -> set[int]:
        integer_slices = self.visit(ctx.integer_slices())
        ring_sizes = set()
        for slice in integer_slices:
            if isinstance(slice, int):
                ring_sizes.add(slice)
            else:
                start, end = slice
                ring_sizes.update(range(start, end + 1))
        return {
            idx
            for idx, atom in self.atom_idx_to_atom.items()
            if any(oechem.OEAtomIsInRingSize(atom, n) for n in ring_sizes)
        }

    def visitAromatic(self, ctx: SelectionLanguageParser.AromaticContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if atom.IsAromatic()}

    def visitAliphatic(self, ctx: SelectionLanguageParser.AliphaticContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if not atom.IsAromatic()}

    def visitHalogen(self, ctx: SelectionLanguageParser.HalogenContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if atom.IsHalogen()}

    def visitMetal(self, ctx: SelectionLanguageParser.MetalContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if atom.IsMetal()}

    def visitIn_ring(self, ctx: SelectionLanguageParser.In_ringContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if atom.IsInRing()}

    def visitWater(self, ctx: SelectionLanguageParser.In_ringContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if oechem.OEIsWater()(atom)}

    def visitBackbone(self, ctx: SelectionLanguageParser.BackboneContext) -> set[int]:
        return {
            idx for idx, atom in self.atom_idx_to_atom.items() if oechem.OEIsBackboneAtom()(atom)
        }

    def visitChiral(self, ctx: SelectionLanguageParser.ChiralContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if atom.IsChiral()}

    def visitAlpha_carbon(self, ctx: SelectionLanguageParser.Alpha_carbonContext) -> set[int]:
        return {idx for idx, atom in self.atom_idx_to_atom.items() if oechem.OEIsCAlpha()(atom)}

    def visitSmarts_expr(self, ctx):
        smarts = ctx.getText()
        sel = set()
        for match in iterate_substructure_matches(target=self.mol, pattern=smarts):
            for atom in match.GetTargetAtoms():
                sel.add(atom.GetIdx())
        return sel

    def visitMolecule(self, ctx):
        # `molecule {smarts}` selects whole molecules that match a given SMARTS.
        smarts = ctx.smarts().getText()
        sel = set()
        for match in iterate_substructure_matches(target=self.mol, pattern=smarts):
            atom_idxs = {atom.GetIdx() for atom in match.GetTargetAtoms()}
            molecule_idxs = {self.atom_idx_to_molecule_idx[idx] for idx in atom_idxs}
            if len(molecule_idxs) > 1:
                continue
            (molecule_idx,) = molecule_idxs
            if sum(
                molecule_idx == molecule_idx2 for molecule_idx2 in self.atom_idx_to_molecule_idx
            ) != len(atom_idxs):
                continue
            sel.update(atom_idxs)
        return sel

    def evaluate_string_atom_property(self, ctx: Any, property: Callable[[oechem.OEAtomBase], str]):
        string_selector = ctx.string_selector()
        if string_selector:
            predicate = self.visit(string_selector)
            # Check that the property is in the set of identifiers
            return {idx for idx, atom in self.atom_idx_to_atom.items() if predicate(property(atom))}
        else:
            # Some queries allow no identifiers, which means 'does this identifier exist'
            # For example, `altloc`, to mean 'does this atom have any altloc'
            return {idx for idx, atom in self.atom_idx_to_atom.items() if property(atom).strip()}

    def evaluate_numeric_atom_property(
        self, ctx: Any, property: Callable[[oechem.OEAtomBase], float]
    ):
        if hasattr(ctx, "integer_selector"):
            number_selector = self.visit(ctx.integer_selector())
        elif hasattr(ctx, "float_selector"):
            number_selector = self.visit(ctx.float_selector())
        else:
            raise ValueError("Property does not support integer_selector or float_selector!")

        def evaluate(atom, number_selector=number_selector):
            atom_property = property(atom)
            return number_selector(atom_property)

        return {atom.GetIdx() for atom in self.mol.GetAtoms() if evaluate(atom)}

    def visitString_selector(self, ctx):
        identifiers = {identifier.getText().strip() for identifier in ctx.IDENTIFIER()} | {
            identifier.getText().strip() for identifier in ctx.POSITIVE_INTEGER()
        }
        strings = {string.getText()[1:-1] for string in ctx.STRING()}

        def predicate(t, identifiers=identifiers, strings=strings):
            return any(t.strip() == identifier for identifier in identifiers) or any(
                t == string for string in strings
            )

        return predicate

    def visitInteger_selector(
        self, ctx: SelectionLanguageParser.Integer_selectorContext
    ) -> Callable[[int], bool]:
        if ctx.integer_slices():
            slices = self.visit(ctx.integer_slices())

            def predicate(t, slices=slices):
                for slice in slices:
                    if isinstance(slice, int):
                        if t == slice:
                            return True
                    else:
                        start, end = slice
                        if t >= start and t <= end:
                            return True

            return predicate
        elif ctx.integer_comparison():
            return self.visit(ctx.integer_comparison())
        else:
            raise ValueError("Invalid integer selector")

    def visitFloat_selector(
        self, ctx: SelectionLanguageParser.Float_selectorContext
    ) -> Callable[[float], bool]:
        if ctx.float_slices():
            slices = self.visit(ctx.float_slices())

            def predicate(t, slices=slices):
                for slice in slices:
                    if isinstance(slice, float):
                        if abs(t - slice) < MAX_DIFFERENCE:
                            return True
                    else:
                        start, end = slice
                        if t >= start and t <= end:
                            return True

            return predicate
        elif ctx.float_comparison():
            return self.visit(ctx.float_comparison())
        elif ctx.integer_selector():
            return self.visit(ctx.integer_selector())
        else:
            raise ValueError("Invalid float selector")

    def visitInteger_slices(
        self, ctx: SelectionLanguageParser.Integer_slicesContext
    ) -> list[int | tuple[int, int]]:
        slices: list[int | tuple[int, int]] = [
            *(self.visit(integer) for integer in ctx.integer()),
            *(self.visit(integer_range) for integer_range in ctx.integer_range()),
        ]
        return slices

    def visitInteger(self, ctx: SelectionLanguageParser.IntegerContext) -> int:
        return int(ctx.getText())

    def visitInteger_range(
        self, ctx: SelectionLanguageParser.Integer_rangeContext
    ) -> tuple[int, int]:
        range = ctx.getText().split(":") if ":" in ctx.getText() else ctx.getText().split("-")
        return int(range[0]), int(range[1])

    def visitInteger_comparison(
        self, ctx: SelectionLanguageParser.Integer_comparisonContext
    ) -> Callable[[int], bool]:
        target: int = self.visit(ctx.integer())
        operator = ctx.OPERATOR().getText()
        match operator:
            case "<":
                return lambda arg: arg < target
            case "<=":
                return lambda arg: arg <= target
            case ">":
                return lambda arg: arg > target
            case ">=":
                return lambda arg: arg >= target
            case "=":
                return lambda arg: arg == target
            case "!=":
                return lambda arg: arg != target

    def visitFloat_slices(
        self, ctx: SelectionLanguageParser.Float_slicesContext
    ) -> list[float | tuple[float, float]]:
        slices: list[float | tuple[float, float]] = [
            *(self.visit(float_) for float_ in ctx.float_()),
            *(self.visit(float_range) for float_range in ctx.float_range()),
        ]
        return slices

    def visitFloat(self, ctx: SelectionLanguageParser.FloatContext) -> float:
        return float(ctx.getText())

    def visitFloat_range(
        self, ctx: SelectionLanguageParser.Float_rangeContext
    ) -> tuple[float, float]:
        range = ctx.getText().split(":") if ":" in ctx.getText() else ctx.getText().split("-")
        return float(range[0]), float(range[1])

    def visitFloat_comparison(
        self, ctx: SelectionLanguageParser.Float_comparisonContext
    ) -> Callable[[float], bool]:
        target: float = self.visit(ctx.float_())
        operator = ctx.OPERATOR().getText()
        match operator:
            case "<":
                return lambda arg: arg < target
            case "<=":
                return lambda arg: arg <= target
            case ">":
                return lambda arg: arg > target
            case ">=":
                return lambda arg: arg >= target
            case "=":
                return lambda arg: abs(arg - target) < MAX_DIFFERENCE
            case "!=":
                return lambda arg: abs(arg - target) > MAX_DIFFERENCE

    def visitOpeneye_atom_expr(self, ctx):
        value = re.search(r"\d+", ctx.getText())
        return {int(value.group(0))}

    def visitOpeneye_atombondset_expr(self, ctx):
        return set.union(*(self.visit(atom) for atom in ctx.openeye_atom_expr()))
