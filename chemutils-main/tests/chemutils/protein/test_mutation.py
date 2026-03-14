import pytest

from chemutils.protein.mutation import (
    DeletionMutation,
    InsertionMutation,
    MissenseMutation,
    parse_mutation_nomenclature,
)


@pytest.mark.parametrize(
    ["nomenclature", "mutations"],
    [
        ("A1G", [MissenseMutation(resname="ALA", resnum=1, replacement_resname="GLY")]),
        ("C12K", [MissenseMutation(resname="CYS", resnum=12, replacement_resname="LYS")]),
        (
            "D18P Y221H",
            [
                MissenseMutation(resname="ASP", resnum=18, replacement_resname="PRO"),
                MissenseMutation(resname="TYR", resnum=221, replacement_resname="HIS"),
            ],
        ),
        (
            "_A12insG",
            [
                InsertionMutation(
                    preceding_resname=None,
                    preceding_resnum=None,
                    subsequent_resname="ALA",
                    subsequent_resnum=12,
                    inserted_resnames=["GLY"],
                )
            ],
        ),
        (
            "C18_P19insHE",
            [
                InsertionMutation(
                    preceding_resname="CYS",
                    preceding_resnum=18,
                    subsequent_resname="PRO",
                    subsequent_resnum=19,
                    inserted_resnames=["HIS", "GLU"],
                )
            ],
        ),
        (
            "R21_insQQE",
            [
                InsertionMutation(
                    preceding_resname="ARG",
                    preceding_resnum=21,
                    subsequent_resname=None,
                    subsequent_resnum=None,
                    inserted_resnames=["GLN", "GLN", "GLU"],
                )
            ],
        ),
        (
            "A24del",
            [
                DeletionMutation(
                    from_resname="ALA",
                    from_resnum=24,
                    to_resname="ALA",
                    to_resnum=24,
                )
            ],
        ),
        (
            "Y205_W220del",
            [
                DeletionMutation(
                    from_resname="TYR",
                    from_resnum=205,
                    to_resname="TRP",
                    to_resnum=220,
                )
            ],
        ),
    ],
)
def test_parse_mutation_nomenclature(nomenclature, mutations):
    results = parse_mutation_nomenclature(nomenclature)
    assert results == mutations
    assert " ".join(result.nomenclature for result in results) == nomenclature
