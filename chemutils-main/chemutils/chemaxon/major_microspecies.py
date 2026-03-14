from typing import overload

from openeye import oechem

from .jvm.chemaxon.marvin.calculations import MajorMicrospeciesPlugin
from .jvm.chemaxon.struc import MolAtom
from .molecule import cxmol_from_oemol, cxmol_from_smiles
from .tautomer import canonical_smiles_from_cxmol


@overload
def get_major_microspecies(
    mol: str,
    /,
    *,
    pH: float = 7.4,
) -> str: ...


@overload
def get_major_microspecies(
    mol: oechem.OEMolBase,
    /,
    *,
    pH: float = 7.4,
) -> oechem.OEMolBase: ...


def get_major_microspecies(
    mol: str | oechem.OEMolBase,
    /,
    *,
    pH: float = 7.4,
) -> str | oechem.OEMolBase:
    """
    Get the major microspecies (protomer) present at a given pH.

    If an OpenEye molecule is provided, the returned OpenEye molecule will be a copy with the hydrogens
    and charges adjusted.

    Args:
        mol: Either a SMILES string or an OpenEye molecule.
        pH: pH to calculate the major microspecies at.

    Returns:
        SMILES string or OpenEye molecule, corresponding to the type of the `mol` input.
    """
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] = {}
    if isinstance(mol, str):
        cxmol = cxmol_from_smiles(mol)
    else:
        returned: oechem.OEMolBase = mol.CreateCopy()
        cxmol = cxmol_from_oemol(returned, atom_mapping=atom_mapping)

    plugin = MajorMicrospeciesPlugin()
    plugin.setpH(float(pH))
    plugin.setMolecule(cxmol)

    plugin.run()

    result_cxmol = plugin.getMajorMicrospecies()

    if isinstance(mol, str):
        return canonical_smiles_from_cxmol(result_cxmol)
    else:
        for cxatom_index, oeatom in enumerate(atom_mapping.keys()):
            cxatom = result_cxmol.getAtom(cxatom_index)
            oeatom.SetFormalCharge(cxatom.getCharge())
            oeatom.SetImplicitHCount(
                cxatom.getImplicitHcount() + cxatom.getExplicitHcount() - oeatom.GetExplicitHCount()
            )
        return returned
