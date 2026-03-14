from typing import Union

from openeye import oechem

from chemutils.molecule import read_molecule_str
from chemutils.molecule.smiles import smiles_from_oemol

from .jvm import isjavainstance
from .jvm.chemaxon.formats import MolExporter, MolImporter
from .jvm.chemaxon.struc import MolAtom, MolBond, Molecule

# Need to use string because Molecule is not a real type
ConvertibleToCXMol = Union[str, oechem.OEMolBase, "Molecule"]


def cxmol_from_oemol(
    oemol: oechem.OEMolBase,
    /,
    *,
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] | None = None,
    bond_mapping: dict[oechem.OEBondBase, MolBond] | None = None,
) -> Molecule:
    # Don't generate random 2D coordinates when doing SDF's, it confuses Chemaxon
    atom_order: list[oechem.OEAtomBase] = []
    smiles = smiles_from_oemol(oemol, atom_order=atom_order)
    cxmol = MolImporter.importMol(
        smiles,
        "SMILES",
    )

    # Get a key that can sort bonds by their occurance in a SMILES
    def bond_key(index1: int, index2: int) -> tuple[int, int]:
        return (max(index1, index2), max(index1, index2) - min(index1, index2))

    if atom_mapping is not None:
        if len(atom_mapping) != 0:
            raise ValueError("`atom_mapping` argument must be an empty dictionary.")
        for oeatom, cxatom in zip(atom_order, cxmol.getAtomArray(), strict=True):
            atom_mapping[oeatom] = cxatom
    if bond_mapping is not None:
        bonds = list(oemol.GetBonds())
        bonds = sorted(
            bonds,
            key=lambda bond: bond_key(
                atom_order.index(bond.GetBgn()), atom_order.index(bond.GetEnd())
            ),
        )
        if len(bond_mapping) != 0:
            raise ValueError("`bond_mapping` argument must be an empty dictionary.")
        for oebond, cxbond in zip(bonds, cxmol.getBondArray(), strict=True):
            bond_mapping[oebond] = cxbond
    return cxmol


def oemol_from_cxmol(cxmol: Molecule, /) -> oechem.OEMolBase:
    return read_molecule_str(MolExporter.exportToFormat(cxmol, "SMILES"), format="SMI")


def cxmol_from_smiles(smiles: str, /) -> Molecule:
    return MolImporter.importMol(smiles, "SMILES")


def smiles_from_cxmol(cxmol: Molecule, /) -> str:
    return MolExporter.exportToFormat(cxmol, "SMILES")


def as_cxmol(
    input: ConvertibleToCXMol,
    /,
    *,
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] | None = None,
    bond_mapping: dict[oechem.OEBondBase, MolBond] | None = None,
) -> Molecule:
    """
    Convert an object to a Chemaxon molecule if it can be.

    This can convert both OpenEye molecules and smiles trings to Chemaxon molecules.

    Args:
        input: Either a Chemaxon molecule or something that can be converted to one.
        atom_mapping: A dictionary to fill with atom mappings.
        bond_mapping: A dictionary to fill with bond mappings.

    Returns:
        Chemaxon molecule.

    Raises:
        ValueError: Cannot convert the given object to a Chemaxon molecule.
    """
    if isinstance(input, oechem.OEMolBase):
        return cxmol_from_oemol(input, atom_mapping=atom_mapping, bond_mapping=bond_mapping)
    elif isinstance(input, str):
        return cxmol_from_smiles(input)
    elif isjavainstance(input, Molecule):
        return input
    else:
        raise ValueError(f"Cannot convert {input} to a Chemaxon Molecule.")
