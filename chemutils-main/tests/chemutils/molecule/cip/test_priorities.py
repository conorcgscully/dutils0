import itertools

import pytest
from openeye import oechem

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.cip.priority import neighbours_ordered_by_cip_priority

FOUR_LIGAND_CASES = [
    ("C", ["H", "C", "Cl", "Br"]),
    ("C", ["H", "C", "CC", "O"]),
    ("C", ["C", "CC", "CBr", "Cl"]),
    ("C", ["H", "[2H]", "C", "O"]),
    ("C", ["H", "C", "[13CH3]", "CC"]),
    # IUPAC Blue Book, P-92.2.1.1.1, Example 1
    ("C", ["H", "F", "Cl", "Br"]),
    # IUPAC Blue Book, P-92.2.1.1.1, Example 2
    ("C", ["H", "C", "OC", "SC"]),
    # IUPAC Blue Book, P-92.2.1.1.1, Example 3
    ("C", ["OC", "[SiH3]", "SC", "[GeH3]"]),
    # IUPAC Blue Book, P-92.2.1.1.2, Example 1
    ("C", ["H", "C", "CO", "Cl"]),
    # IUPAC Blue Book, P-92.2.1.1.2, Example 2
    ("C", ["H", "CO", "CCl", "Cl"]),
    # IUPAC Blue Book, P-92.2.1.1.2, Example 3
    ("C", ["H", "CO", "C(C)Cl", "Cl"]),
    # IUPAC Blue Book, P-92.2.1.1.3, Example 1
    ("C", ["H", "C(CCI)C(F)CF", "C(CCBr)C(F)CCl", "O"]),
    # IUPAC Blue Book, P-92.2.1.1.3, Example 2
    ("C", ["H", "C(CCBr)C(F)CF", "C(CCI)C(F)CF", "O"]),
    # IUPAC Blue Book, P-92.2.1.2, Example 1
    ("C", ["H", "CC", "C=C", "O"]),
    # IUPAC Blue Book, P-92.2.1.2, Example 2
    ("C", ["H", "C=C", "C=O", "O"]),
    # IUPAC Blue Book, P-92.2.1.2, Example 3
    ("C", ["H", "C=C", "C#N", "O"]),
    # IUPAC Blue Book, P-92.2.1.2, Example 4
    ("C", ["H", "CO", "C=O", "O"]),
    # IUPAC Blue Book, P-92.2.1.3, Example 1
    ("C", ["H", "C(C)C", "C1CC1", "O"]),
    # IUPAC Blue Book, P-92.2.1.3, Example 2
    ("C", ["H", "C1CC1", "C2CCC2", "O"]),
    # IUPAC Blue Book, P-92.2.1.3, Example 3
    ("C", ["H", "C1CC1", "C=O", "O"]),
    # IUPAC Blue Book, P-92.2.1.3, Example 4
    ("C", ["H", "CC", "C1CC1", "Cl"]),
    # IUPAC Blue Book, P-92.2.2 Sequence Subrule 1b: Priority due to duplicate atoms
    # IUPAC Blue Book, P-92.2.2, Example 1
    ("C", ["H", "C(CCC1CC1)(CCC2CC2)CCC3CC3", "C45CCC(CC4)CC5", "O"]),
    # IUPAC Blue Book, P-92.2.2.3, Example 1
    ("C", ["H", "[2H]", "C", "O"]),
    # IUPAC Blue Book, P-92.2.2.3, Example 2
    ("C", ["H", "C[125I]", "CI", "O"]),
    # IUPAC Blue Book, P-92.4.2.1, Example 1
    ("C", ["H", "CCC/C=C/C", r"CCC/C=C\C", "O"]),
    # IUPAC Blue Book, P-92.4.2.1, Example 2
    ("C", ["H", r"CCC/C(Cl)=C\C", "CCC/C(Cl)=C/C", "O"]),
    # IUPAC Blue Book, P-92.4.2.2, Example 2
    ("C", ["H", r"/C=C(\C/C=C\C)C/C=C\C", r"/C=C(/C/C=C\C)C/C=C/C", "O"]),
    # IUPAC Blue Book, P-93.2.2, Example 1
    ("[N+]", ["C", "CC=C", "Cc1ccccc1", "c2ccccc2"]),
    # IUPAC Blue Book, P-93.2.3, Example 1
    ("[N+]", ["C", "CC", "c1ccccc1", "[O-]"]),
    # IUPAC Blue Book, P-93.2.3, Example 2
    ("P", ["C", "CCC", "c1ccccc1", "=O"]),
    # IUPAC Blue Book, P-93.2.4, Example 1
    ("P", ["[O-]", "OC", "Oc1ccccc1", "=S"]),
    # IUPAC Blue Book, P-93.2.4, Example 2
    ("P", ["=O", "OC", "Oc1ccccc1", "[S-]"]),
    # IUPAC Blue Book, P-93.2.4, Example 3
    ("P", ["C", "=O", "OC(C)C", "[S-]"]),
    # IUPAC Blue Book, P-93.2.4, Example 4
    ("P", ["C", "[O-]", "OC(C)C", "=S"]),
    # IUPAC Blue Book, P-93.2.4, Example 5
    ("P", ["C", "c1ccccc1", "=O", "OC"]),
    # IUPAC Blue Book, P-93.2.4, Example 6
    ("P", ["=O", "[17O-]", "[18O-]", "OC"]),
    # IUPAC Blue Book, P-93.2.5, Example 1
    ("S", ["c1ccccc1", "=O", "OC", "=S"]),
    # IUPAC Blue Book, P-93.2.5, Example 2
    ("S", ["Cc1ccccc1", "c2cc(C)ccc2", "=O", "=[18O]"]),
    # IUPAC Blue Book, P-93.2.5, Example 3
    ("S", ["[O-]", "=[17O]", "=[18O]", "Oc1ccccc1"]),
    # IUPAC Blue Book, P-93.2.6, Example 1
    ("[Si]", ["H", "C", "CCC", "O"]),
]

PERMUTATIONS4 = itertools.permutations([0, 1, 2, 3])


@pytest.mark.parametrize("permutation", PERMUTATIONS4)
@pytest.mark.parametrize(["root_smiles", "ligand_smiles_list"], FOUR_LIGAND_CASES)
def test_four_ligand(root_smiles, ligand_smiles_list, permutation):
    # Randomly ordered ligands
    a = ligand_smiles_list[permutation[0]]
    b = ligand_smiles_list[permutation[1]]
    c = ligand_smiles_list[permutation[2]]
    d = ligand_smiles_list[permutation[3]]

    oemol = oemol_from_smiles(f"{root_smiles}({a})({b})({c}){d}")

    root_atom = next(iter(oemol.GetAtoms()))

    neighbour_atoms = [bond.GetNbr(root_atom) for bond in root_atom.GetBonds()]
    correct_order = [
        neighbour_atoms[permutation.index(0)],
        neighbour_atoms[permutation.index(1)],
        neighbour_atoms[permutation.index(2)],
        neighbour_atoms[permutation.index(3)],
    ]
    result = neighbours_ordered_by_cip_priority(root_atom)
    assert result == correct_order


# Stereocenter is labelledID 1, neighbours are labelled in increasing priority ID 2-5
CYCLIC_CASES = [
    # IUPAC Blue Book, P-92.2.1.3, Example 5
    ("[C:1]1([H:2])([CH:3])[CH:5](C)CCC[CH2:4]1"),
    ("[CH:5]1(C)[C:1]([H:2])([CH3:3])[CH2:4]CCC1"),
    # IUPAC Blue Book, P-92.2.2.2, Example 1
    ("[C:1]([H:2])([O:5])([C:4]1(CCC2CC1)CC2)[C:3](CCC3CC3)(CCC4CC4)CCC5CC5"),
    # IUPAC Blue Book, P-92.2.2.2, Example 2
    ("[C:1]1([H:2])([CH2:3]C=C[CH:5]12)[CH2:4]2"),
    # IUPAC Blue Book, P-93.2.2, Example 2
    ("[P+:1]1([CH2:2]CC2=CC=CC=C2[CH2:3]1)([c:4]3ccccc3)[c:5]4ccc(O)cc4"),
]


@pytest.mark.parametrize("smiles", CYCLIC_CASES)
def test_cyclic(smiles):
    oemol = oemol_from_smiles(smiles)

    root_atom = oemol.GetAtom(oechem.OEHasMapIdx(1))

    ordered_neighbours = neighbours_ordered_by_cip_priority(root_atom)
    expected = [
        oemol.GetAtom(oechem.OEHasMapIdx(2)),
        oemol.GetAtom(oechem.OEHasMapIdx(3)),
        oemol.GetAtom(oechem.OEHasMapIdx(4)),
        oemol.GetAtom(oechem.OEHasMapIdx(5)),
    ]
    assert ordered_neighbours == expected
