from __future__ import annotations

from collections import defaultdict
from collections.abc import MutableMapping
from typing import TYPE_CHECKING, Literal, overload

from openeye import oechem

from chemutils.molecule import oemol_from_smiles, smiles_from_oemol

from .jvm.chemaxon.marvin.calculations import TautomerizationPlugin, pKaPlugin
from .jvm.chemaxon.struc import BondType
from .molecule import ConvertibleToCXMol, as_cxmol, smiles_from_cxmol

if TYPE_CHECKING:
    from .jvm.chemaxon.struc import MolAtom, MolBond, Molecule


def get_chemaxon_bond_order(bond: MolBond, /) -> int:
    # Not a constant so that no java communication occurs
    CXBONDTYPE_TO_ORDER = {
        BondType.SINGLE: 1,
        BondType.DOUBLE: 2,
        BondType.TRIPLE: 3,
    }

    return CXBONDTYPE_TO_ORDER[bond.getBondType()]


@overload
def _fetch(
    *,
    cxmol: Molecule,
    input: str,
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] | None = ...,
    bond_mapping: dict[oechem.OEBondBase, MolBond] | None = ...,
) -> str: ...


@overload
def _fetch(
    *,
    cxmol: Molecule,
    input: Molecule,
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] | None = ...,
    bond_mapping: dict[oechem.OEBondBase, MolBond] | None = ...,
) -> Molecule: ...


@overload
def _fetch(
    *,
    cxmol: Molecule,
    input: oechem.OEMolBase,
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] | None = ...,
    bond_mapping: dict[oechem.OEBondBase, MolBond] | None = ...,
) -> oechem.OEMolBase: ...


def _fetch(
    *,
    cxmol: Molecule,
    input: str | oechem.OEMolBase | Molecule,
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] | None = None,
    bond_mapping: dict[oechem.OEBondBase, MolBond] | None = None,
) -> str | oechem.OEMolBase | Molecule:
    if isinstance(input, str):
        return canonical_smiles_from_cxmol(cxmol)
    elif isinstance(input, oechem.OEMolBase):
        if atom_mapping is None or bond_mapping is None:
            raise RuntimeError("Cannot convert from oemol to cxmol without mappings.")
        cxmol.dearomatize()
        if [atom.GetIdx() for atom in input.GetAtoms()] != list(range(input.NumAtoms())):
            raise RuntimeError("Cannot use OpenEye molecule with non-sequential atom indices.")
        returned: oechem.OEMolBase = input.CreateCopy()
        for oeatom, cxatom in atom_mapping.items():
            oeatom = returned.GetAtom(oechem.OEHasAtomIdx(oeatom.GetIdx()))
            cxatom = cxmol.getAtom(cxatom.getParent().indexOf(cxatom))
            oeatom.SetFormalCharge(cxatom.getCharge())
            oeatom.SetImplicitHCount(
                cxatom.getImplicitHcount() + cxatom.getExplicitHcount() - oeatom.GetExplicitHCount()
            )
        for oebond, cxbond in bond_mapping.items():
            oebond = returned.GetBond(oechem.OEHasBondIdx(oebond.GetIdx()))
            cxbond = cxmol.getBond(cxbond.getParent().indexOf(cxbond))
            oebond.SetOrder(get_chemaxon_bond_order(cxbond))
        oechem.OEAssignAromaticFlags(returned)
        return returned
    return cxmol


def canonical_smiles_from_cxmol(cxmol: Molecule) -> str:
    """Get a OpenEye standardised SMILES from a ChemAxon molecule."""
    cxmol.dearomatize()
    smiles = smiles_from_cxmol(cxmol)
    oemol = oemol_from_smiles(smiles)
    oechem.OEAssignAromaticFlags(oemol)
    return smiles_from_oemol(oemol)


@overload
def _get_tautomers(
    input: str,
    /,
    *,
    distribution: Literal[False],
    canonical: bool = False,
    major: bool = False,
    normal: bool = False,
    pH: float | None = None,
) -> list[str]: ...


@overload
def _get_tautomers(
    input: str,
    /,
    *,
    distribution: Literal[True],
    canonical: bool = False,
    major: bool = False,
    normal: bool = False,
    pH: float | None = None,
) -> dict[str, float]: ...


@overload
def _get_tautomers(
    input: Molecule,
    /,
    *,
    distribution: Literal[False],
    canonical: bool = False,
    major: bool = False,
    normal: bool = False,
    pH: float | None = None,
) -> list[Molecule]: ...


@overload
def _get_tautomers(
    input: Molecule,
    /,
    *,
    distribution: Literal[True],
    canonical: bool = False,
    major: bool = False,
    normal: bool = False,
    pH: float | None = None,
) -> dict[Molecule, float]: ...


@overload
def _get_tautomers(
    input: oechem.OEMolBase,
    /,
    *,
    distribution: Literal[False],
    canonical: bool = False,
    major: bool = False,
    normal: bool = False,
    pH: float | None = None,
) -> list[oechem.OEMolBase]: ...


@overload
def _get_tautomers(
    input: oechem.OEMolBase,
    /,
    *,
    distribution: Literal[True],
    canonical: bool = False,
    major: bool = False,
    normal: bool = False,
    pH: float | None = None,
) -> dict[oechem.OEMolBase, float]: ...


def _get_tautomers(
    input: ConvertibleToCXMol,
    /,
    *,
    distribution: bool,
    canonical: bool = False,
    major: bool = False,
    normal: bool = False,
    pH: float | None = None,
) -> (
    list[str]
    | list[oechem.OEMolBase]
    | list[Molecule]
    | dict[str, float]
    | dict[oechem.OEMolBase, float]
    | dict[Molecule, float]
):
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] = {}
    bond_mapping: dict[oechem.OEAtomBase, MolBond] = {}
    cxmol = as_cxmol(input, atom_mapping=atom_mapping, bond_mapping=bond_mapping)

    plugin = TautomerizationPlugin()

    plugin.setMolecule(cxmol)
    if pH is not None:
        plugin.setpH(float(pH))
    plugin.setTakeCanonicalForm(canonical)
    plugin.setTakeMajorTautomer(major)
    plugin.setNormalTautomerGenerationMode(normal)
    plugin.setDominantTautomerDistributionCalculation(distribution)

    plugin.run()

    if distribution:
        return {
            _fetch(
                cxmol=plugin.getStructure(i),
                input=input,
                atom_mapping=atom_mapping,
                bond_mapping=bond_mapping,
            ): distr
            for i in range(plugin.getStructureCount())
            if (distr := plugin.getDominantTautomerDistribution(i) / 100.0) > 0.00005
        }
    else:
        return [
            _fetch(
                cxmol=plugin.getStructure(i),
                input=input,
                atom_mapping=atom_mapping,
                bond_mapping=bond_mapping,
            )
            for i in range(plugin.getStructureCount())
        ]


@overload
def get_canonical_tautomer(cxmol: str, /, *, normal: bool = False) -> str: ...
@overload
def get_canonical_tautomer(cxmol: Molecule, /, *, normal: bool = False) -> Molecule: ...
@overload
def get_canonical_tautomer(
    cxmol: oechem.OEMolBase, /, *, normal: bool = False
) -> oechem.OEMolBase: ...


def get_canonical_tautomer(
    cxmol: ConvertibleToCXMol, /, *, normal: bool = False
) -> str | oechem.OEMolBase | Molecule:
    """
    Get the canonical tautomer for a given molecule.

    This corresponds to the `cxcalc canonicaltautomer ...` command.

    Args:
        cxmol: SMILES, OpenEye molecule or ChemAxon molecule.
        normal: Should normal tautomers be generated. This corresponds to the `--normal` argument.

    Returns:
        OpenEye-standardised SMILES of the canonical tautomer.
    """
    (tautomer,) = _get_tautomers(cxmol, distribution=False, canonical=True, normal=normal)
    return tautomer


@overload
def get_major_tautomer(cxmol: str, /) -> str: ...
@overload
def get_major_tautomer(cxmol: Molecule, /) -> Molecule: ...
@overload
def get_major_tautomer(cxmol: oechem.OEMolBase, /) -> oechem.OEMolBase: ...


def get_major_tautomer(cxmol: ConvertibleToCXMol, /) -> ConvertibleToCXMol:
    """
    Get the major (most prevalent) tautomer for a given molecule.

    This corresponds to the `cxcalc majortautomer ...` command.

    Args:
        cxmol: SMILES, OpenEye molecule or ChemAxon molecule.

    Returns:
        OpenEye-standardised SMILES of the major tautomer.
    """
    (tautomer,) = _get_tautomers(cxmol, distribution=False, major=True)
    return tautomer


@overload
def get_tautomer_distribution(cxmol: str, /) -> dict[str, float]: ...
@overload
def get_tautomer_distribution(cxmol: Molecule, /) -> dict[Molecule, float]: ...
@overload
def get_tautomer_distribution(cxmol: oechem.OEMolBase, /) -> dict[oechem.OEMolBase, float]: ...


def get_tautomer_distribution(
    cxmol: ConvertibleToCXMol, /
) -> dict[str, float] | dict[Molecule, float] | dict[oechem.OEMolBase, float]:
    """
    Get the distribution of tautomers for a given molecule.

    This corresponds to the `cxcalc dominanttautomerdistribution ...` command.

    Args:
        cxmol: SMILES, OpenEye molecule or ChemAxon molecule.

    Returns:
        Mapping of OpenEye-standardised SMILES to their fraction within the distribution.
    """
    return _get_tautomers(cxmol, distribution=True)


@overload
def get_major_tautomer_and_protomer(cxmol: str, /, *, pH: float = 7.4) -> str: ...
@overload
def get_major_tautomer_and_protomer(cxmol: Molecule, /, *, pH: float = 7.4) -> Molecule: ...
@overload
def get_major_tautomer_and_protomer(
    cxmol: oechem.OEMolBase, /, *, pH: float = 7.4
) -> oechem.OEMolBase: ...


def get_major_tautomer_and_protomer(
    cxmol: ConvertibleToCXMol, /, *, pH: float = 7.4
) -> str | oechem.OEMolBase | Molecule:
    """
    Get the major (most prevalent) tautomer for a given molecule, including protomers.

    This corresponds to the `cxcalc majortautomer --pH <pH> ...` command.

    Args:
        cxmol: SMILES, OpenEye molecule or ChemAxon molecule.
        pH: pH to calculate the major tautomer at.

    Returns:
        OpenEye-standardised SMILES of the major tautomer.
    """
    (tautomer,) = _get_tautomers(cxmol, distribution=False, major=True, pH=pH)
    return tautomer


def get_tautomer_and_protomer_distribution(
    cxmol: ConvertibleToCXMol, /, *, pH: float = 7.4
) -> dict[str, float]:
    """
    Get the distribution of tautomers for a given molecule, including protomers.

    This corresponds to the `cxcalc dominanttautomerdistribution ...` command.

    Args:
        cxmol: SMILES, OpenEye molecule or ChemAxon molecule.
        pH: pH to calculate the distribution at.

    Returns:
        Mapping of OpenEye-standardised SMILES to their fraction within the distribution.
    """
    return _get_tautomers(cxmol, distribution=True, pH=pH)


@overload
def get_protomer_distribution(
    input: str,
    /,
    *,
    pH: float = 7.4,
    temperature: float | None = None,
) -> dict[str, float]: ...


@overload
def get_protomer_distribution(
    input: Molecule,
    /,
    *,
    pH: float = 7.4,
    temperature: float | None = None,
) -> dict[Molecule, float]: ...


@overload
def get_protomer_distribution(
    input: oechem.OEMolBase,
    /,
    *,
    pH: float = 7.4,
    temperature: float | None = None,
) -> dict[oechem.OEMolBase, float]: ...


def get_protomer_distribution(
    input: ConvertibleToCXMol,
    /,
    *,
    pH: float = 7.4,
    temperature: float | None = None,
) -> dict[str, float] | dict[oechem.OEMolBase, float] | dict[Molecule, float]:
    """
    Get the distribution of protomers (charge states) at a given pH.

    This does not include tautomerization. To include this, use `get_tautomer_and_protomer_distribution`.

    Args:
        input: Either a SMILES or an OpenEye molecule.
        pH: pH to calculate microspecies distribution at.
        temperature: Temperature to calculate at, in Kelvin.

    Returns:
        Mapping of microspecies (either as SMILES or OpenEye molecules, depending on input)
        to the calculated fraction that this microspecies exists.
    """
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] = {}
    bond_mapping: dict[oechem.OEBondBase, MolBond] = {}
    cxmol = as_cxmol(input, atom_mapping=atom_mapping, bond_mapping=bond_mapping)

    plugin = pKaPlugin()
    if temperature:
        plugin.setTemperature(float(temperature))

    plugin.setpH(float(pH))
    plugin.setMolecule(cxmol)

    assert plugin.isMsCalc()
    plugin.run()

    result: MutableMapping[ConvertibleToCXMol, float] = defaultdict(lambda: 0)

    for i in range(plugin.getMsCount()):
        # Deduplicate by SMILES
        cxmol = plugin.getMsMolecule(i)
        key = _fetch(cxmol=cxmol, input=input, atom_mapping=atom_mapping, bond_mapping=bond_mapping)
        if (distr := plugin.getSingleMsDistribution(i) / 100.0) >= 0.00005:
            result[key] += distr

    return dict(sorted(result.items(), key=lambda item: item[1], reverse=True))


def get_major_protomer(cxmol: ConvertibleToCXMol, /, *, pH: float = 7.4) -> str:
    """
    Get the major (most prevalent) protomer for a given molecule.

    Args:
        cxmol: SMILES, OpenEye molecule or ChemAxon molecule.
        pH: pH to calculate the major tautomer at.

    Returns:
        OpenEye-standardised SMILES of the major protomer.
    """
    # Order by percentage, going from highest to lowest
    protomers = sorted(
        get_protomer_distribution(cxmol, pH=pH).items(), key=lambda kv: kv[1], reverse=True
    )
    return protomers[0][0]  # type: ignore
