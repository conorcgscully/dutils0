from openeye import oechem

from chemutils.molecule.selection import subset_molecule
from chemutils.molecule.selection.subset import add_internal_bonds_to_selection

# SMARTS pattern for marking atoms
# Quotes taken from Ertl 2017
# `mark all heteroatoms in a molecule, including halogens`
HETEROATOM_SMARTS = "!#6"
# `[carbon] atoms connected by non-aromatic double or triple bond to any heteroatom`
# `[carbon] atoms in nonaromatic carbon-carbon double or triple bonds`
CARBON_MULTIBOND_SMARTS = "#6$(*=,#*)"
# `acetal carbons, i.e. sp3 carbons connected to two or more oxygens, nitrogens or sulfurs; these O, N or S atoms must have only single bonds`
CARBON_ACETAL_SMARTS = "#6$(*([OX2,NX3,SX2!$(*=*)])[OX2,NX3,SX2!$(*=*)])!$(*=*)"
# `all atoms in oxirane, aziridine and thiirane rings (such rings are traditionally considered to be functional groups due to their high reactivity).`
CARBON_HETEROTHREECYCLE_SMARTS = "C$(*1[O,N,S]C1)"
FG_SMARTS = f"[{HETEROATOM_SMARTS},{CARBON_MULTIBOND_SMARTS},{CARBON_ACETAL_SMARTS},{CARBON_HETEROTHREECYCLE_SMARTS}]"

SUBSEARCH = oechem.OESubSearch(FG_SMARTS)
assert SUBSEARCH.IsValid()


def get_ertl_functional_groups(oemol: oechem.OEGraphMol, /) -> oechem.OEAtomBondSetVector:
    """
    Get the subsets of a molecule that represent functional groups.

    This is based on the method of Ertl, which is described in the following paper:

    Ertl, P. An algorithm to identify functional groups in organic molecules. J Cheminform 9, 36 (2017)
    https://doi.org/10.1186/s13321-017-0225-z

    This groups connected heteroatoms and double-bonded atoms, assigning to each
    a pseudosmiles label that represents the functional group.

    This function returns a vector of these subsets as `OEAtomBondSet` objects. Each one
    has a role of the form `ertl:pseudosmiles:...`.

    Args:
        oemol: The molecule to analyze.

    Returns:
        Functional groups as an `oechem.OEAtomBondSetVector`.
    """
    subset = oechem.OEAtomBondSet()
    # Mark all atoms that match Ertl's criteria
    for match in SUBSEARCH.Match(oemol):
        for atom in match.GetTargetAtoms():
            subset.AddAtom(atom)
    # Add bonds between marked atoms, to group them into functional groups
    subset = add_internal_bonds_to_selection(subset)
    # Aromatic bonds are not included
    for bond in list(subset.GetBonds()):
        if bond.IsAromatic():
            subset.RemoveBond(bond)

    # Create a molecular subset with the marked atoms and bonds
    # This will give a fragmented molecule, where each submolecule is
    # an Ertl functional group
    atom_mapping: dict[oechem.OEAtomBase, oechem.OEAtomBase] = {}
    bond_mapping: dict[oechem.OEBondBase, oechem.OEBondBase] = {}
    submol = subset_molecule(
        oemol,
        selection=subset,
        include_internal_bonds=False,
        bond_handling="hydrogens",
        atom_mapping=atom_mapping,
        bond_mapping=bond_mapping,
    )
    atom_mapping_reversed = {v: k for k, v in atom_mapping.items()}
    bond_mapping_reversed = {v: k for k, v in bond_mapping.items()}

    # Resolve atomic environment (adding R-groups describing the environment)
    for src_atom, dest_atom in list(atom_mapping.items()):
        # Convert aromatic atoms to R-groups with names `[Car]`, `[Nar]`, etc.
        if dest_atom.IsAromatic():
            dest_atom.SetName(oechem.OEGetAtomicSymbol(dest_atom.GetAtomicNum()) + "ar")
            dest_atom.SetAtomicNum(0)
            continue
        # Non-carbonyl carbons don't have environment atoms
        # `environments on carbon atoms are deleted, the only exception are substituents on carbonyl that are retained (to distinguish between aldehydes and ketones)`
        if dest_atom.GetAtomicNum() == 6 and not oechem.OEHasDoubleBondO(dest_atom):
            continue
        # Don't need to convert hydrogens to R-groups if there's no hydrogens
        if dest_atom.GetImplicitHCount() == 0:
            continue
        # Determine how many hydrogens should *not* be replaced by R-groups
        maintained_hydrogens = 0
        # `hydrogens on the -OH groups`
        # `hydrogens on the simple amines and thiols (i.e. FGs with just single central N or S atom)`
        if (
            (dest_atom.GetAtomicNum() == 6 and oechem.OEHasDoubleBondO(dest_atom))
            or (dest_atom.GetAtomicNum() == 8)
            or (dest_atom.GetAtomicNum() in [7, 16] and dest_atom.GetHvyDegree() == 0)
        ):
            maintained_hydrogens = src_atom.GetImplicitHCount()

        # Generally, `R` is used for these R-group
        label = "R"
        # For amines and alcohols with a single carbon neighbour, replace the C with `[Car]` or `[Cal]`
        # `environments on single atomic N or O FGs with one carbon connected, where this carbon is retained also with its type (aliphatic or aromatic)`
        if (
            dest_atom.GetAtomicNum() in [7, 8]
            and dest_atom.GetHvyDegree() == 0
            and src_atom.GetHvyDegree() == 1
        ):
            label = "Car" if next(iter(src_atom.GetAtoms())).IsAromatic() else "Cal"

        # Replace all hydrogens (above `maintained_hydrogens`) with R-groups
        for _ in range(maintained_hydrogens, dest_atom.GetImplicitHCount()):
            rgroup = submol.NewAtom(0)
            rgroup.SetName(label)
            bond = submol.NewBond(dest_atom, rgroup, 1)
        dest_atom.SetImplicitHCount(maintained_hydrogens)

    # Get pseudosmiles for each functional group
    full_smiles = oechem.OECreateSmiString(
        submol, oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_SuperAtoms
    )
    smiles_list = full_smiles.split(".")

    _, atom_idx_to_component_idx = oechem.OEDetermineComponents(submol)

    # Divide up the molecule into oechem.OEAtomBondSet's for each functional group
    subsets_to_atombondset: dict[int, oechem.OEAtomBondSet] = {}
    for atom in oechem.OEGetSmiStringOrder(
        submol, oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_SuperAtoms
    ):
        component_idx = atom_idx_to_component_idx[atom.GetIdx()]
        if component_idx not in subsets_to_atombondset:
            subset_smiles = smiles_list[len(subsets_to_atombondset)]
            subsets_to_atombondset[component_idx] = oechem.OEAtomBondSet()
            subsets_to_atombondset[component_idx].AddRole(f"ertl:pseudosmiles:{subset_smiles}")
        if (orig_atom := atom_mapping_reversed.get(atom)) is not None:
            subsets_to_atombondset[component_idx].AddAtom(orig_atom)
    for bond in submol.GetBonds():
        atom = bond.GetBgn()
        component_idx = atom_idx_to_component_idx[atom.GetIdx()]
        if (orig_bond := bond_mapping_reversed.get(bond)) is not None:
            subsets_to_atombondset[component_idx].AddBond(orig_bond)

    return oechem.OEAtomBondSetVector(list(subsets_to_atombondset.values()))
