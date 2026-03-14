import numpy as np
import numpy.typing as npt
from openeye import oechem

# SMARTS describing symmetric groups not catched by openeye -> charge delocalization
SYM_SMARTS = (
    # Carboxylic acids(deprotonated or protonated)
    ["[CX3](=O)[OX1H0-,OX2H1]"],
    # Phosphoric acid(deprotonated or protonated, both representations, no substituents - pure P(=O)(OH)3, should basically never happen)
    [
        "P(=[OX1])([$([OX2H]),$([OX1-])])([$([OX2H]),$([OX1-])])([$([OX2H]),$([OX1-])])",
        "[P+]([OX1-])([$([OX2H]),$([OX1-])])([$([OX2H]),$([OX1-])])([$([OX2H]),$([OX1-])])",
    ],
    # Phosphoric acid (deprotonated or protonated, both representations, one substituent)
    [
        "[PX4](=[OX1])([$([OX2H]),$([OX1-])])([$([OX2H]),$([OX1-])])",
        "[P+X4]([OX1-])([$([OX2H]),$([OX1-])])([$([OX2H]),$([OX1-])])",
    ],
    # Phosphoric acid (deprotonated or protonated, both representations, two substituents)
    ["[PX4](=[OX1])([$([OX2H]),$([OX1-])])", "[P+X4]([OX1-])([$([OX2H]),$([OX1-])])"],
    # Boronic acid (deprotonated or protonated)
    ["[BX3]([$([OX2H]),$([OX1-])])([$([OX2H]),$([OX1-])])"],
    # Amidinium #TODO this does not handle the case where the amidinium group is symmetric due to identical nitrogen substituents. No idea how to do this with SMARTS...
    ["[NH2][CX3]=[NH2+]"],
    # Guanidinium #TODO this does not handle the case of a pure Guanidinium ion, neglecting to re-symmetrize the last nitrogen. However, that is a single-molecule case.
    ["[CX3](=[NH2+])([NH2])([NX3])"],
    # nitro (also covering nitrates)
    ["[NX3](=[OX1])(=[OX1])", "[NX3+]([OX1-])(=[OX1])"],
    # TODO consider doing enumeration with OEGetReasonableTautomers
)

SYM_IDXS = ((1, 2), (1, 2, 3, 4), (1, 2, 3), (1, 2), (1, 2), (0, 2), (1, 2), (1, 2))


def order_symmetry_classes(classes: npt.NDArray[np.int64], /) -> npt.NDArray[np.int64]:
    """
    Reindex symmetry classes such that they are indexed in order.

    For example, the symmetry classes [3, 1, 5, 2, 3] will be converted to [1, 2, 3, 4, 1].

    Args:
        classes: A NumPy array of shape (N, ) containing a set of integers 0..M, where M <= N.

    Returns:
        Array of equal size with symmetry classes reordered.
    """
    _, unique_index, unique_inverse = np.unique(classes, return_index=True, return_inverse=True)
    return np.argsort(np.argsort(unique_index))[unique_inverse] + 1


def get_symmetry_classes(mol: oechem.OEGraphMol, /) -> npt.NDArray[np.int64]:
    """
    Compute the symmetry class of each atom in a molecule.

    Two atoms share a symmetry class if they can be exchanged and result in the same molecule.
    For example, the three hydrogens of a methyl group will have the same symmetry class.

    Args:
        mol: Molecule to get symmetry classes for.

    Returns:
        NumPy array of integers, indicating the symmetry class of each atom.
    """
    oechem.OEPerceiveSymmetry(mol)

    symmetry_classes_unordered = [atom.GetSymmetryClass() for atom in mol.GetAtoms()]

    # Get mapping between atom indices and indices in symlist (a result of ignoring hydrogens for symlist but not for atom indexing)
    atom_idx_to_unordered_idx = {
        aidx: smidx
        for smidx, aidx in enumerate([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1])
    }

    # patch symmetry classes belonging to uncatched groups due to charge delocalization
    for smarts_group, indices in zip(SYM_SMARTS, SYM_IDXS, strict=True):
        for smarts in smarts_group:
            searcher = oechem.OESubSearch(smarts)

            # Get matches
            match_iterator = searcher.Match(mol)

            # Iterate through matches (there may be several e.g. carboxylic acids) and matched atoms in each match
            for match in match_iterator:
                atoms_to_symmetrize = []
                for matchpair in match.GetAtoms():
                    pattern_idx = matchpair.pattern.GetIdx()

                    # Collect indices in oemol of atoms corresponding to the symmetric ones for each SMARTS
                    if pattern_idx in indices:
                        atoms_to_symmetrize.append(matchpair.target.GetIdx())

                for atom in mol.GetAtoms():
                    # Only consider non-hydrogens
                    if atom.GetAtomicNum() != 1:
                        atom_idx = atom.GetIdx()
                        symlist_idx = atom_idx_to_unordered_idx[atom_idx]
                        if atom_idx in atoms_to_symmetrize:
                            # Get symmetry token of first atom in symmetry group
                            root_symmetry_atom_idx = atoms_to_symmetrize[0]
                            # Get symlist index of root atom
                            root_symlist_idx = atom_idx_to_unordered_idx[root_symmetry_atom_idx]
                            # Set symmetry token of this atom to that of root atom
                            symmetry_classes_unordered[symlist_idx] = symmetry_classes_unordered[
                                root_symlist_idx
                            ]

    return order_symmetry_classes(np.array(symmetry_classes_unordered))
