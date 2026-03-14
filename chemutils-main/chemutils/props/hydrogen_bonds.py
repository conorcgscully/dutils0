from openeye import oechem, oemolprop

from chemutils.molecule.hydrogen_bonds import (
    get_num_acceptor_lone_pairs,
    get_num_donor_hydrogens,
    get_num_hydrogen_acceptors,
    get_num_hydrogen_donors,
)

from .property import MolecularProperty


def get_lipinski_hydrogen_donors(mol: oechem.OEGraphMol) -> int:
    """
    Get the number of hydrogen bond donors based on the definition of Lipinski (1997).

    Lipinski defined a hydrogen donor as any nitrogen or oxygen with an attached hydrogen.
    """
    return oemolprop.OEGetLipinskiDonorCount(mol)  # type: ignore


def get_lipinski_hydrogen_acceptors(mol: oechem.OEGraphMol) -> int:
    """
    Get the number of hydrogen bond donors based on the definition of Lipinski (1997).

    Lipinski defined a hydrogen acceptor as any nitrogen or oxygen.
    """
    return oemolprop.OEGetLipinskiAcceptorCount(mol)  # type: ignore


NumberLipinskiHydrogenDonors = MolecularProperty(
    id="num_lipinski_hydrogen_donors", oemol_func=get_lipinski_hydrogen_donors
)

NumberLipinskiHydrogenAcceptors = MolecularProperty(
    id="num_lipinski_hydrogen_acceptors", oemol_func=get_lipinski_hydrogen_acceptors
)

NumberHydrogenDonors = MolecularProperty(
    id="num_hydrogen_donors", oemol_func=get_num_hydrogen_donors
)

NumberDonorHydrogens = MolecularProperty(
    id="num_donor_hydrogens", oemol_func=get_num_donor_hydrogens
)

NumberHydrogenAcceptors = MolecularProperty(
    id="num_hydrogen_acceptors", oemol_func=get_num_hydrogen_acceptors
)

NumberAcceptorLonePairs = MolecularProperty(
    id="num_acceptor_lone_pairs", oemol_func=get_num_acceptor_lone_pairs
)
