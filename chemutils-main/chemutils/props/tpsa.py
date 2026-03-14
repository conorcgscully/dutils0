from openeye import oechem, oemolprop

from chemutils.props.property import MolecularProperty


def get_2d_tpsa(mol: oechem.OEGraphMol) -> float:
    """
    Get the two-dimensional topological polar surface area, in angstroms squared.

    This is based on the method of Ertl, Rohde & Selzer (2000), using per-fragment contributions.

    References:
        Ertl, Rohde & Selzer, J. Med. Chem. 2000, 43, 20, 3714-3717
    """
    return oemolprop.OEGet2dPSA(mol)  # type: ignore


TwoDimensionalTPSA = MolecularProperty(id="tpsa_2d", oemol_func=get_2d_tpsa)
