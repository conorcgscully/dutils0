from openeye import oechem, oemolprop
from rdkit import Chem

from chemutils.chemaxon.logp import get_logp
from chemutils.molecule.charge import remove_formal_charges
from chemutils.rdmol_from_oemol import as_rdmol

from .property import MolecularProperty


def get_oexlogp(mol: oechem.OEGraphMol) -> float | None:
    """
    Calculate logP using OpenEye's modified version of XLogP.

    logP measures the partition coefficient of a molecule in octanol and water. The XLogP method was
    originally descibred by Wang, Fu and Lai (1997), and consists of atomic contributions with some
    intramolecular corrections.

    The OpenEye implementation makes some modifications as descibed in
    https://docs.eyesopen.com/toolkits/cpp/molproptk/molprops.html#molprops-logp, hence the result
    of this function should be referred to as OEXLogP, and not XLogP.

    Args:
        mol: Molecule to compute OEXLogP for.

    Returns:
        OEXLogP value, or None if it could not be computed for this molecule.
    """
    if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        mol = mol.CreateCopy()
        remove_formal_charges(mol)
    atomXLogP = oechem.OEFloatArray(mol.GetMaxAtomIdx())
    result = oemolprop.OEGetXLogPResult(mol, atomXLogP)
    if result.IsValid():
        return result.GetValue()  # type: ignore
    return None


def get_wildman_crippen_logp(mol: Chem.Mol | oechem.OEMolBase) -> float:
    """
    Calculate logP using the approach of Wildman and Crippen (1999).

    This is a per-atom approach to calculating logP, with each atom assigned to a
    specific type which has a contribution to the final value.

    Hydrogens are ignored for this calculation.
    """
    mol = as_rdmol(mol)
    return Chem.Crippen.MolLogP(mol)  # type: ignore


LogP = MolecularProperty("logp", cxmol_func=get_logp)

OEXLogP = MolecularProperty(id="oexlogp", oemol_func=get_oexlogp)

WildmanCrippenLogP = MolecularProperty(
    id="wildmancrippen_logp", rdkit_func=get_wildman_crippen_logp
)
