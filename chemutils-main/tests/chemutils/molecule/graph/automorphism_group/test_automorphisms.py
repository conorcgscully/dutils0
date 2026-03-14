from typing import NamedTuple

import numpy as np
import pytest
from openeye import oechem

from chemutils.molecule import (
    get_automorphism_group,
    get_num_automorphisms,
    iterate_automorphism_matches,
    oemol_from_smiles,
)


class TestCase(NamedTuple):
    smiles: str
    num_automorphisms: int
    num_automorphism_subgroups: int


CASES = [
    # Benzene
    TestCase("c1ccccc1", 12, 1),
    # Cyclobutane
    TestCase("C1CCC1", 8, 1),
    # Cubane
    TestCase("C12C3C4C1C5C2C3C45", 48, 1),
    # Buckyball
    TestCase(
        "c12c3c4c5c2c2c6c7c1c1c8c3c3c9c4c4c%10c5c5c2c2c6c6c%11c7c1c1c7c8c3c3c8c9c4c4c9c%10c5c5c2c2c6c6c%11c1c1c7c3c3c8c4c4c9c5c2c2c6c1c3c42",
        120,
        1,
    ),
    # Ensure aromatic and non-aromatic rings are not exchanged
    TestCase("c1ccccc1-C1CCCCC1", 4, 2),
    # Difficult examples from the PDB
    TestCase(
        "C1(C(C(C(C(C1OS(=O)(=O)O)OS(=O)(=O)O)OS(=O)(=O)O)OS(=O)(=O)O)OS(=O)(=O)O)OS(=O)(=O)O",
        559872,
        7,
    ),
    TestCase(
        "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@H]2[C@@H]([C@H]([C@@H](O[C@H]2C(=O)O)O[C@@H]3[C@H](O[C@@H]([C@@H]([C@H]3O)NS(=O)(=O)O)O[C@H]4[C@@H]([C@H]([C@@H](O[C@H]4C(=O)O)O)OS(=O)(=O)O)O)COS(=O)(=O)O)OS(=O)(=O)O)O)NS(=O)(=O)O)O)O[C@@H]5[C@@H]([C@H]([C@@H]([C@@H](O5)C(=O)O)O)O)OS(=O)(=O)O)OS(=O)(=O)O",
        2239488,
        10,
    ),
    TestCase(
        "CO[C@@H]1[C@H](O[C@@H]([C@@H]([C@H]1OC)OS(=O)(=O)O)O[C@H]2[C@@H]([C@H]([C@@H](O[C@@H]2C(=O)O)O[C@@H]3[C@H](O[C@@H]([C@@H]([C@H]3OS(=O)(=O)O)OS(=O)(=O)O)O[C@H]4[C@@H]([C@H]([C@@H](O[C@H]4C(=O)O)O[C@@H]5[C@H](O[C@@H]([C@@H]([C@H]5OS(=O)(=O)O)OS(=O)(=O)O)OC)COS(=O)(=O)O)OS(=O)(=O)O)OC)COS(=O)(=O)O)OC)OC)COS(=O)(=O)O",
        40310784,
        11,
    ),
]


@pytest.mark.parametrize(["smiles", "num_automorphisms", "num_automorphism_subgroups"], CASES)
def test_num_automorphisms(smiles, num_automorphisms, num_automorphism_subgroups):
    oemol = oemol_from_smiles(smiles)
    num_automorphisms_direct = get_num_automorphisms(oemol)
    assert num_automorphisms_direct == num_automorphisms

    group = get_automorphism_group(oemol)
    assert num_automorphisms == group.num_automorphisms
    assert num_automorphism_subgroups == group.num_subgroups

    # Only test explicit lists if there aren't too many
    if num_automorphisms <= 10_000:
        automorphisms = list(group.iterate_automorphisms())

        assert num_automorphisms == len(automorphisms)

        chemutils_automorphisms = [
            np.array([atom.GetIdx() for atom in match.GetTargetAtoms()])
            for match in iterate_automorphism_matches(oemol)
        ]

        for automorphism in automorphisms:
            assert automorphism.shape == (oemol.NumAtoms(),)

        assert {tuple(automorphism) for automorphism in automorphisms} == {
            tuple(automorphism) for automorphism in chemutils_automorphisms
        }

        # OpenEye can only calculate up to 1024 automorphisms
        if num_automorphisms <= 1024:
            openeye_automorphisms = [
                np.array([atom.GetIdx() for atom in match.GetTargetAtoms()])
                for match in oechem.OEGetAutomorphs(oemol)
            ]

            assert num_automorphisms == len(openeye_automorphisms)

            assert {tuple(automorphism) for automorphism in automorphisms} == {
                tuple(automorphism) for automorphism in openeye_automorphisms
            }

        assert num_automorphisms == len(chemutils_automorphisms)
