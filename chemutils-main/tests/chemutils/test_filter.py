import pathlib
import tempfile

import pytest

from chemutils.filter import (
    FilterResult,
    charmtx_openeye_filter,
    openeye_filter,
    rdkit_filter,
)
from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.rdmol_from_oemol import rdmol_from_oemol

from .props.drug_smiles import DRUG_SMILES


@pytest.mark.parametrize("name", ["BRENK", "NIH", "PAINS", "PAINS_A", "PAINS_B", "PAINS_C", "ZINC"])
def test_rd_filter_exists(name):
    _ = rdkit_filter(name)


@pytest.mark.parametrize(
    "name",
    [
        "BLOCKBUSTER",
        "DRUG",
        "FRAGMENT",
        "LEAD",
        "PAINS",
    ],
)
def test_oe_filter_exists(name):
    _ = openeye_filter(name)


def test_rdkit_filter_invalid():
    with pytest.raises(ValueError):
        _ = rdkit_filter("INVALID")


def test_oe_filter_invalid():
    with pytest.raises(ValueError):
        _ = openeye_filter("INVALID")


@pytest.mark.parametrize("smiles", [*DRUG_SMILES.values()])
def test_oe_filter_passed(smiles):
    oemol = oemol_from_smiles(smiles)
    filter = openeye_filter("BLOCKBUSTER")
    assert filter.passed(oemol) == filter(oemol).passed


@pytest.mark.parametrize("smiles", [*DRUG_SMILES.values()])
def test_rd_filter_passed(smiles):
    oemol = oemol_from_smiles(smiles)
    rdmol = rdmol_from_oemol(oemol)
    filter = rdkit_filter("BRENK")
    assert filter.passed(rdmol) == filter(rdmol).passed


@pytest.mark.parametrize(
    ["smiles", "rules"],
    [
        # Based on examples from Figure 1, Oprea (2000)
        (
            "CCS(Cl)(=O)=O",
            {
                "atom count",
                "sulfonyl_halide",
                "molecular weight",
                "heteroatom to carbon ratio",
                "carbons",
            },
        ),
        ("CC(=O)Cl", {"atom count", "acid_halide", "molecular weight", "carbons"}),
        ("CCCCl", {"atom count", "molecular weight", "carbons", "heteroatoms", "alkyl_halide"}),
        ("Clc1ncccn1", {"atom count", "halopyrimidine", "molecular weight", "carbons"}),
        ("CC(=O)OC(=O)C", {"molecular weight", "carbons", "atom count", "anhydride"}),
        ("CC(=O)C(=O)C", {"atom count", "molecular weight", "carbons"}),
        (
            "CCC(=O)C(Cl)(Cl)Cl",
            {
                "atom count",
                "molecular weight halide fraction",
                "alphahalo_ketone",
                "molecular weight",
                "carbons",
                "perhalo_ketone",
                "alkyl_halide",
            },
        ),
        ("CCC(=O)C", {"carbons", "atom count", "heteroatoms", "molecular weight"}),
        ("CC1OC1C", {"carbons", "molecular weight", "atom count", "heteroatoms", "epoxide"}),
        ("CC1N(C)C1C", {"atom count", "molecular weight", "carbons", "heteroatoms", "aziridine"}),
        ("CCC(=O)OC", {"atom count", "molecular weight", "carbons"}),
        ("CCC(=O)SC", {"atom count", "molecular weight", "carbons"}),  # Seems OpenEye misses this
        (
            "CCS(=O)(=O)OC",
            {
                "atom count",
                "molecular weight",
                "heteroatom to carbon ratio",
                "carbons",
                "sulfonic_ester",
            },
        ),
        (
            "CCP(=O)(O)OC",
            {
                "atom count",
                "molecular weight",
                "heteroatom to carbon ratio",
                "phosphonic_acid",
                "carbons",
                "phosphonic_ester",
            },
        ),
        ("CCC(C)=NC", {"atom count", "heteroatoms", "molecular weight", "carbons"}),
        ("CC=O", {"atom count", "molecular weight", "aldehyde", "carbons", "heteroatoms"}),
        (
            "CC(=O)C=CC",
            {"atom count", "molecular weight", "carbons", "heteroatoms", "michael_acceptor"},
        ),
        (
            "CC(=O)CC(Cl)C",
            {"atom count", "molecular weight", "beta_halo_carbonyl", "carbons", "alkyl_halide"},
        ),
        ("COOC", {"atom count", "molecular weight", "carbons"}),  # Seems OpenEye misses this
        ("CNOC", {"atom count", "molecular weight", "carbons"}),  # Seems OpenEye misses this
        ("CSOC", {"atom count", "molecular weight", "carbons"}),  # Seems OpenEye misses this
        ("CNNC", {"atom count", "hydrazine", "molecular weight", "carbons"}),
        ("CNSC", {"atom count", "acyclic_NS", "molecular weight", "carbons"}),
        ("CSSC", {"atom count", "molecular weight", "carbons", "disulfide"}),
    ],
)
def test_drug_filter(smiles, rules):
    oemol = oemol_from_smiles(smiles)
    oe_filter = openeye_filter("DRUG")
    assert oe_filter(oemol) == FilterResult(passed=False, failed_rules=rules)


@pytest.mark.parametrize(
    ["smiles", "openeye_rules", "rdkit_rules"],
    [
        ("O=CCl", {"pains_pre_acid_halide", "pains_pre_aldehyde"}, set()),
        (
            "Oc1ccccc1C=NN",
            {"pains_pre_hydrazine", "pains_a_hzone_phenol_A"},
            {"hzone_phenol_A(479)"},
        ),
    ],
)
def test_pains_filter(smiles, openeye_rules, rdkit_rules):
    oemol = oemol_from_smiles(smiles)
    rdmol = rdmol_from_oemol(oemol)
    oe_filter = openeye_filter("PAINS")
    rd_filter = rdkit_filter("PAINS")
    assert oe_filter(oemol) == FilterResult(
        passed=openeye_rules == set(), failed_rules=openeye_rules
    )
    assert rd_filter(rdmol) == FilterResult(passed=rdkit_rules == set(), failed_rules=rdkit_rules)


@pytest.mark.parametrize(
    ["drug_name", "openeye_rules", "rdkit_rules"],
    [
        ("azidothymidine", {"pains_a_azo_A"}, {"azo_A(324)"}),
        ("ciprofloxacin", {"pains_a_anil_di_alk_B"}, set()),
        ("cimetidine", {"pains_c_imidazole_B"}, {"imidazole_B(2)"}),
        ("olsalazine", {"pains_a_azo_A"}, {"azo_A(324)"}),
        ("sulfasalazine", {"pains_a_azo_A"}, {"azo_A(324)"}),
    ],
)
def test_pains_filter_drugs(drug_name, openeye_rules, rdkit_rules):
    smiles = DRUG_SMILES[drug_name]
    oemol = oemol_from_smiles(smiles)
    rdmol = rdmol_from_oemol(oemol)
    oe_filter = openeye_filter("PAINS")
    rd_filter = rdkit_filter("PAINS")
    # Ignore pains_pre part for OpenEye, as it is not part of PAINS really
    assert {
        rule for rule in oe_filter(oemol).failed_rules if not rule.startswith("pains_pre")
    } == openeye_rules
    assert rd_filter(rdmol).failed_rules == rdkit_rules


@pytest.mark.parametrize("filter_name", ["PAINS", "BLOCKBUSTER", "DRUG", "LEAD"])
def test_oefilter_passall(filter_name):
    smiles = "O=C(N)Cc1ccc(cc1)OCC(O)CNC(C)C"
    oemol = oemol_from_smiles(smiles)
    oe_filter = openeye_filter(filter_name)
    assert oe_filter(oemol) == FilterResult(passed=True, failed_rules=set())


@pytest.mark.parametrize(["smiles", "valid"], [("C", False), ("CC", True), ("CCC", False)])
def test_charmtx_oefilter(monkeypatch, smiles, valid):
    with tempfile.NamedTemporaryFile() as f:
        input_file_string = b"""
MIN_CARBONS      2      "Minimum number of carbons"
MAX_CARBONS      2      "Maximum number of carbons"
        """
        f.write(input_file_string)
        f.flush()
        FIXTURE = {
            "test": pathlib.Path(f.name),
        }
        monkeypatch.setattr("chemutils.filter.CHARMTX_OE_FILTER_NAME_TO_FILE", FIXTURE)
        oemol = oemol_from_smiles(smiles)
        oe_filter = charmtx_openeye_filter("test")
        assert oe_filter.passed(oemol) == valid


def test_charmtx_oefilter_invalid_filter_name_fails():
    with pytest.raises(ValueError):
        charmtx_openeye_filter("CHARM_20221101_invalid")


def test_charmtx_oefilter_filter_file_not_found_fails(monkeypatch):
    FIXTURE = {
        "test": pathlib.Path("does_not_exist.txt"),
    }
    monkeypatch.setattr("chemutils.filter.CHARMTX_OE_FILTER_NAME_TO_FILE", FIXTURE)
    with pytest.raises(ValueError, match="Could not open"):
        charmtx_openeye_filter("test")


def test_oefilter_repeated_use():
    good_smiles = "O=C(N)Cc1ccc(cc1)OCC(O)CNC(C)C"
    bad_smiles = "O=C1NC(C(C)=CN1[C@@H]2O[C@H](CO)[C@@H](N=[N+]=[N-])C2)=O"
    fail_reasons = {"pains_a_azo_A", "pains_pre_azido", "pains_pre_azo"}
    oe_filter = openeye_filter("PAINS")
    assert oe_filter(oemol_from_smiles(good_smiles)) == FilterResult(
        passed=True, failed_rules=set()
    )
    assert oe_filter(oemol_from_smiles(bad_smiles)) == FilterResult(
        passed=False, failed_rules=fail_reasons
    )
    assert oe_filter(oemol_from_smiles(good_smiles)) == FilterResult(
        passed=True, failed_rules=set()
    )
    assert oe_filter(oemol_from_smiles(bad_smiles)) == FilterResult(
        passed=False, failed_rules=fail_reasons
    )
