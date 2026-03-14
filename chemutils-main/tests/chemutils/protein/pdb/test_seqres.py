from pathlib import Path

import pytest

from chemutils import fs
from chemutils.protein.pdb.modify import filter_lines_in_pdb, insert_lines_into_pdb
from chemutils.protein.pdb.seqres import (
    generate_seqres,
    parse_seqres_lines,
    read_seqres_from_pdb,
)

CASES = fs.read_yaml(Path(__file__).parent / "seqres.yaml")


@pytest.mark.parametrize(
    ["seqres", "sequences"], [(case["seqres"], case["sequences"]) for case in CASES]
)
def test_read_seqres_from_pdb(seqres, sequences):
    assert read_seqres_from_pdb(seqres) == sequences


@pytest.mark.parametrize(
    ["seqres", "sequences"], [(case["seqres"], case["sequences"]) for case in CASES]
)
def test_parse_seqres_lines(seqres, sequences):
    assert parse_seqres_lines(seqres.splitlines()) == sequences


@pytest.mark.parametrize(
    ["seqres", "sequences"], [(case["seqres"], case["sequences"]) for case in CASES]
)
def test_generate_seqres(seqres, sequences):
    assert generate_seqres(sequence_residue_names=sequences) == seqres.splitlines()


@pytest.mark.parametrize(
    ["seqres", "sequences"], [(case["seqres"], case["sequences"]) for case in CASES]
)
def test_add_seqres_to_pdb(seqres, sequences):
    pdb_before = "REMARK   1\n"
    pdb_existing = "SEQRES I_WAS_HERE_BEFORE\n"
    pdb_after = "ATOM  "

    pdb_str = pdb_before + pdb_existing + pdb_after

    pdb_str = filter_lines_in_pdb(pdb=pdb_str, filter=lambda line: not line.startswith("SEQRES"))

    pdb_str = insert_lines_into_pdb(
        pdb=pdb_str, lines=generate_seqres(sequence_residue_names=sequences)
    )
    assert pdb_str == pdb_before + seqres + pdb_after
