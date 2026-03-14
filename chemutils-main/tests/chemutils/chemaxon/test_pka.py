from dataclasses import dataclass

import pytest

from chemutils.chemaxon.pka import (
    get_fraction_microspecies_uncharged,
    get_macro_pka,
    get_macro_pka_conjugate_acid,
    get_microspecies_distribution,
    get_num_acidic_atoms,
    get_num_basic_atoms,
)
from chemutils.molecule import oemol_from_smiles


@dataclass
class pKaClass:
    smiles: str
    microspecies: set[str]
    distr_neutral: dict[str, float]
    distr_acid: dict[str, float]
    distr_alkaline: dict[str, float]
    pKa: float | None
    pKa_conjugate_acid: float | None
    fraction_uncharged: float


CASES = [
    # Water
    pKaClass(
        smiles="O",
        microspecies={"O", "[OH-]"},
        distr_neutral={"O": 1.0},
        distr_acid={"O": 1.0},
        distr_alkaline={"O": 0.9998, "[OH-]": 0.0002},
        pKa=15.70,
        pKa_conjugate_acid=-1.80,
        fraction_uncharged=1.0,
    ),
    # Nitrogen
    pKaClass(
        smiles="N",
        microspecies={"N", "[NH4+]"},
        distr_neutral={"N": 0.0334, "[NH4+]": 0.9666},
        distr_acid={"[NH4+]": 1.0},
        distr_alkaline={"N": 0.9993, "[NH4+]": 0.0007},
        pKa=None,
        pKa_conjugate_acid=8.86,
        fraction_uncharged=0.0334,
    ),
    pKaClass(
        smiles="C(=O)O",
        microspecies={"C(=O)O", "C(=O)[O-]"},
        distr_neutral={"C(=O)O": 0.0007, "C(=O)[O-]": 0.9993},
        distr_acid={"C(=O)O": 0.9946, "C(=O)[O-]": 0.0054},
        distr_alkaline={"C(=O)[O-]": 1.0},
        pKa=4.27,
        pKa_conjugate_acid=None,
        fraction_uncharged=0.0007,
    ),
    pKaClass(
        smiles="OS(O)(=O)=O",
        microspecies={"OS(O)(=O)=O", "OS(=O)(=O)[O-]", "[O-]S(=O)(=O)[O-]"},
        distr_neutral={"[O-]S(=O)(=O)[O-]": 1.0},
        distr_acid={"OS(=O)(=O)[O-]": 0.4407, "[O-]S(=O)(=O)[O-]": 0.5593},
        distr_alkaline={"[O-]S(=O)(=O)[O-]": 1.0},
        pKa=-3.03,
        pKa_conjugate_acid=None,
        fraction_uncharged=0.0,
    ),
    pKaClass(
        smiles="C[C@H](N)C(O)=O",
        microspecies={
            "C[C@H](N)C(O)=O",
            "C[C@@H](C(=O)O)[NH3+]",
            "C[C@@H](C(=O)[O-])N",
            "C[C@@H](C(=O)[O-])[NH3+]",
        },
        distr_neutral={
            "C[C@@H](C(=O)[O-])N": 0.0083,
            "C[C@@H](C(=O)[O-])[NH3+]": 0.9917,
        },
        distr_acid={
            "C[C@@H](C(=O)O)[NH3+]": 0.7490,
            "C[C@@H](C(=O)[O-])[NH3+]": 0.2510,
        },
        distr_alkaline={
            "C[C@@H](C(=O)[O-])N": 0.9970,
            "C[C@@H](C(=O)[O-])[NH3+]": 0.0029,
        },
        pKa=2.47,
        pKa_conjugate_acid=9.477,
        fraction_uncharged=0.0,
    ),
    pKaClass(
        smiles="N[C@@H](CC(O)=O)C(O)=O",
        microspecies={
            "N[C@@H](CC(O)=O)C(O)=O",
            "C([C@@H](C(=O)[O-])N)C(=O)O",
            "C([C@@H](C(=O)[O-])N)C(=O)[O-]",
            "C([C@@H](C(=O)[O-])[NH3+])C(=O)O",
            "C([C@@H](C(=O)[O-])[NH3+])C(=O)[O-]",
            "C([C@@H](C(=O)O)[NH3+])C(=O)O",
            "C([C@@H](C(=O)O)[NH3+])C(=O)[O-]",
        },
        distr_neutral={
            "C([C@@H](C(=O)[O-])N)C(=O)O": 0.0001,
            "C([C@@H](C(=O)[O-])N)C(=O)[O-]": 0.0060,
            "C([C@@H](C(=O)[O-])[NH3+])C(=O)O": 0.0051,
            "C([C@@H](C(=O)[O-])[NH3+])C(=O)[O-]": 0.9888,
        },
        distr_acid={
            "C([C@@H](C(=O)O)[NH3+])C(=O)O": 0.3350,
            "C([C@@H](C(=O)[O-])[NH3+])C(=O)O": 0.6641,
            "C([C@@H](C(=O)O)[NH3+])C(=O)[O-]": 0.0004,
            "C([C@@H](C(=O)[O-])[NH3+])C(=O)[O-]": 0.0005,
        },
        distr_alkaline={
            "C([C@@H](C(=O)[O-])N)C(=O)[O-]": 0.9959,
            "C([C@@H](C(=O)[O-])[NH3+])C(=O)[O-]": 0.0041,
        },
        pKa=1.70,
        pKa_conjugate_acid=9.61,
        fraction_uncharged=0.0,
    ),
]


@pytest.mark.parametrize(
    ["smiles", "pH", "expected"],
    [(smiles, 7.4, case.distr_neutral) for case in CASES for smiles in sorted(case.microspecies)]
    + [(smiles, 2, case.distr_acid) for case in CASES for smiles in sorted(case.microspecies)]
    + [(smiles, 12, case.distr_alkaline) for case in CASES for smiles in sorted(case.microspecies)],
)
def test_get_microspecies_distribution(smiles, pH, expected):
    distr = get_microspecies_distribution(smiles, pH=pH)
    distr = {key: value for key, value in distr.items() if value >= 0.00005}
    assert distr == pytest.approx(expected, abs=0.0001)


@pytest.mark.parametrize(
    ["smiles", "fraction"],
    [(smiles, case.fraction_uncharged) for case in CASES for smiles in sorted(case.microspecies)],
)
def test_get_fraction_microspecies_uncharged(smiles, fraction):
    assert get_fraction_microspecies_uncharged(smiles) == pytest.approx(fraction, abs=0.0001)


def test_get_fraction_microspecies_uncharged_dearomatize():
    # Previously was failing, if we don't dearomatize the cxmol before conversion
    assert get_fraction_microspecies_uncharged(
        oemol_from_smiles("CB1C=Cc2cccc(OCC(O)CN(C)C)c2N1")
    ) == pytest.approx(0.04787, abs=0.00001)


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [(smiles, case.pKa) for case in CASES for smiles in sorted(case.microspecies)],
)
def test_get_macro_pka(smiles, expected):
    assert get_macro_pka(smiles) == pytest.approx(expected, abs=0.01)


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [(smiles, case.pKa_conjugate_acid) for case in CASES for smiles in sorted(case.microspecies)],
)
def test_get_macro_pka_conjugate_acid(smiles, expected):
    assert get_macro_pka_conjugate_acid(smiles) == pytest.approx(expected, abs=0.01)


@pytest.mark.parametrize(
    ["smiles", "num_acidic", "num_basic"],
    [
        # Alanine
        ("C[C@@H](C(=O)O)N", 1, 1),
        # Arginine
        ("C(C[C@@H](C(=O)O)N)CNC(=N)N", 1, 2),
        # Asparagine
        ("C(N)C[C@H](N)C(=O)O", 1, 2),
        # Aspartate
        ("C([C@@H](C(=O)O)N)C(=O)O", 2, 1),
        # Cysteine
        ("C([C@@H](C(=O)O)N)S", 2, 1),
        # Glutamine
        ("O=C(N)CCC(N)C(=O)O", 1, 1),
        # Glutamtic Acid
        ("C(CC(=O)O)[C@@H](C(=O)O)N", 2, 1),
        # Glycine
        ("C(C(=O)O)N", 1, 1),
        # Histidine
        ("O=C([C@H](CC1=CNC=N1)N)O", 1, 1),
        # Isoleucine
        ("CC[C@H](C)[C@@H](C(=O)O)N", 1, 1),
        # Leucine
        ("CC(C)C[C@@H](C(=O)O)N", 1, 1),
        # Lysine
        ("C(CCN)C[C@@H](C(=O)O)N", 1, 2),
        # Methionine
        ("CSCC[C@H](N)C(=O)O", 1, 1),
        # Phenylalanine
        ("c1ccc(cc1)C[C@@H](C(=O)O)N", 1, 1),
        # Proline
        ("C1C[C@H](NC1)C(=O)O", 1, 1),
        # Serine
        ("C([C@@H](C(=O)O)N)O", 1, 0),
        # Threonine
        ("C[C@H]([C@@H](C(=O)O)N)O", 1, 1),
        # Tryptophan
        ("c1[nH]c2ccccc2c1C[C@H](N)C(=O)O", 1, 1),
        # Tyrosine
        ("N[C@@H](Cc1ccc(O)cc1)C(O)=O", 1, 1),
        # Valine
        ("CC(C)[C@@H](C(=O)O)N", 1, 1),
        # Acetic acid
        ("CC(=O)O", 1, 0),
        # Oxalic acid
        ("OC(=O)C(=O)O", 2, 0),
        # Phosphoric acid
        ("OP(=O)(O)O", 3, 0),
        # Ammonia
        ("N", 0, 1),
        # Aniline
        ("c1ccccc1N", 0, 0),
        # Pyrrole
        ("[nH]1cccc1", 0, 0),
        # TREN
        ("NCCN(CCN)CCN", 0, 4),
        # Urea
        ("C(=O)(N)N", 0, 0),
        # Examples from Hannah
        ("NCCCN1CCCCC1", 0, 2),
        ("CCCCN2CCCCC2", 0, 1),
        ("NCCC(N3CCCCC3)=O", 0, 1),
        ("OC1=C(F)C(F)=C(F)C(F)=C1F", 1, 0),
        ("O=C(OC1)NC1=O", 1, 0),
        ("CS(NC(C)=O)(=O)=O", 1, 0),
    ],
)
def test_num_acidic_basic_atoms(smiles, num_acidic, num_basic):
    assert get_num_acidic_atoms(smiles) == num_acidic
    assert get_num_basic_atoms(smiles) == num_basic
