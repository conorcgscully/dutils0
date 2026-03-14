# Tests that ensure the matches returned are not going to segfault
# The segfaults are caused by the deletion of the original query
# molecule after the iterators have returned.

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.graph import (
    iterate_automorphism_matches,
    iterate_homomorphism_matches,
    iterate_maximum_common_substructure_matches,
    iterate_substructure_matches,
)


def test_segfault_iterate_substructure_matches():
    pattern = oemol_from_smiles("CCC")
    target = oemol_from_smiles("CCC")

    # Cast to list to force generator to finish
    matches = list(iterate_substructure_matches(pattern=pattern, target=target))

    for match in matches:
        for pattern_atom in match.GetPatternAtoms():
            assert isinstance(
                str(pattern_atom), str
            )  # Sefaults without the fix that remaps matches
            assert pattern_atom in pattern.GetAtoms()
        for pattern_bond in match.GetPatternBonds():
            assert pattern_bond in pattern.GetBonds()
        for target_atom in match.GetTargetAtoms():
            assert target_atom in target.GetAtoms()
        for target_bond in match.GetTargetBonds():
            assert target_bond in target.GetBonds()


def test_segfault_iterate_automorphism_matches():
    oemol = oemol_from_smiles("CCC")

    # Cast to list to force generator to finish
    matches = list(iterate_automorphism_matches(oemol))

    for match in matches:
        for pattern_atom in match.GetPatternAtoms():
            assert isinstance(
                str(pattern_atom), str
            )  # Sefaults without the fix that remaps matches
            assert pattern_atom in oemol.GetAtoms()
        for pattern_bond in match.GetPatternBonds():
            assert pattern_bond in oemol.GetBonds()
        for target_atom in match.GetTargetAtoms():
            assert target_atom in oemol.GetAtoms()
        for target_bond in match.GetTargetBonds():
            assert target_bond in oemol.GetBonds()


def test_segfault_iterate_homomorphism_matches():
    pattern = oemol_from_smiles("CCC")
    target = oemol_from_smiles("CCC")

    # Cast to list to force generator to finish
    matches = list(iterate_homomorphism_matches(pattern=pattern, target=target))

    for match in matches:
        for pattern_atom in match.GetPatternAtoms():
            assert isinstance(
                str(pattern_atom), str
            )  # Sefaults without the fix that remaps matches
            assert pattern_atom in pattern.GetAtoms()
        for pattern_bond in match.GetPatternBonds():
            assert pattern_bond in pattern.GetBonds()
        for target_atom in match.GetTargetAtoms():
            assert target_atom in target.GetAtoms()
        for target_bond in match.GetTargetBonds():
            assert target_bond in target.GetBonds()


def test_segfault_iterate_maximum_common_substructure_matches():
    pattern = oemol_from_smiles("CCC")
    target = oemol_from_smiles("CCC")

    # Cast to list to force generator to finish
    matches = list(iterate_maximum_common_substructure_matches(pattern=pattern, target=target))

    for match in matches:
        for pattern_atom in match.GetPatternAtoms():
            assert isinstance(
                str(pattern_atom), str
            )  # Sefaults without the fix that remaps matches
            assert pattern_atom in pattern.GetAtoms()
        for pattern_bond in match.GetPatternBonds():
            assert pattern_bond in pattern.GetBonds()
        for target_atom in match.GetTargetAtoms():
            assert target_atom in target.GetAtoms()
        for target_bond in match.GetTargetBonds():
            assert target_bond in target.GetBonds()
