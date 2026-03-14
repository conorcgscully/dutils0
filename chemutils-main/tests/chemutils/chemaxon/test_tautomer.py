from dataclasses import dataclass

import pytest

from chemutils.chemaxon.tautomer import (
    get_canonical_tautomer,
    get_major_protomer,
    get_major_tautomer,
    get_major_tautomer_and_protomer,
    get_protomer_distribution,
    get_tautomer_and_protomer_distribution,
    get_tautomer_distribution,
)
from chemutils.molecule import (
    assert_mol_equal,
    make_hydrogens_implicit,
    oemol_from_smiles,
    smiles_from_oemol,
)
from chemutils.molecule.equal import (
    AllAtomProperties,
    AllBondProperties,
    AtomProperty,
    BondProperty,
)

AMINE = "CN"
AMINE_CATION = "C[NH3+]"

CARBOXYYLIC_ACID = "CC(=O)O"
CARBOXYYLIC_ACID_ANION = "CC(=O)[O-]"

ENOL = "CC=C(C)O"
ENOL_ANION = "CC=C(C)[O-]"
KETO = "CCC(=O)C"

DIENOL = "CC(=CC(=O)C)O"
DIKETO = "CC(=O)CC(=O)C"
DIENOL_OXYGEN_ANION = "CC(=CC(=O)C)[O-]"
DIKETO_CARBON_ANION = "CC(=O)[CH-]C(=O)C"

AMIDE = "CC(=O)NC"
IMIDIC_ACID = "CC(=NC)O"
IMIDIC_ACID_ANION = "CC(=NC)[O-]"
IMIDIC_ACID_CATION = "CC(=[NH+]C)O"
IMIDIC_ACID_ZWITTERION = "CC(=[NH+]C)[O-]"

LACTAM = "C1CC(=O)NC1"
LACTIM = "C1CC(N=C1)O"

TWOPYRIDONE = "c1cc[nH]c(=O)c1"
TWOPYRIDONE_ANION = "c1cc[n-]c(=O)c1"
TWOHYDROXYPYRIDINE = "c1ccnc(c1)O"
TWOHYDROXYPYRIDINE_ANION = "c1ccnc(c1)[O-]"

IMINE = "CC=N"
IMINE_CATION = "CC=[NH2+]"
ENAMINE = "C=CN"
ENAMINE_CATION = "C=C[NH3+]"

CYANAMIDE = "CCN(C)C#N"
CARBODIIMIDE = "CC[N](=C=N)C"
CARBODIIMIDE_CATION = "C(=[NH2+])=NC"

ALANINE = "C[C@@H](C(=O)O)N"
ALANINE_ANION = "C[C@@H](C(=O)[O-])N"
ALANINE_ZWITTERION = "C[C@@H](C(=O)[O-])[NH3+]"

AZIDE = "CN=[N+]=[N-]"
AZIDE_TRIPLE = "C[N-][N+]#N"
AZIDE_CATION = "CN=[N+]=[NH2+]"
AZIDE_TRIPLE_CATION = "CN[N+]#N"

TRIAZOLE_1H = "c1[nH]ncn1"
TRIAZOLE_1H_ANION = "c1[n-]ncn1"
TRIAZOLE_4H = "c1[nH]cnn1"
TRIAZOLE_4H_ANION = "c1[n-]cnn1"

GUANIDINE1 = "CCCNC(=NC)NCC"
GUANIDINE1_CATION = "CCCNC(=[NH+]C)NCC"
GUANIDINE2 = "CCCNC(=NCC)NC"
GUANIDINE2_CATION = "CCCNC(=[NH+]CC)NC"
GUANIDINE3 = "CCCN=C(NC)NCC"
GUANIDINE3_CATION = "CCC[NH+]=C(NC)NCC"

NITROSO = "CCN=O"
OXIME = "CC=NO"
OXIME_CATION = "CC=[NH+]O"
OXIME_ANION = "CC=N[O-]"


@dataclass
class TautomerCase:
    smiles: str
    canonical_tautomer: str
    major_tautomer: str
    tautomer_distribution: dict[str, str]
    normal_canonical_tautomer: str | None = None
    major_tautomer_protomer: str | None = None
    major_protomer: str | None = None
    protomer_distribution: dict[str, str] | None = None
    tautomer_protomer_distribution: dict[str, str] | None = None


CASES = [
    # Simple Amine (no tautomerisation, cation form preferred)
    TautomerCase(
        smiles=AMINE,
        canonical_tautomer=AMINE,
        major_tautomer=AMINE,
        major_tautomer_protomer=AMINE_CATION,
        major_protomer=AMINE_CATION,
        tautomer_distribution={AMINE: "1.0000"},
        protomer_distribution={AMINE_CATION: "0.9979", AMINE: "0.0021"},
        tautomer_protomer_distribution={AMINE_CATION: "0.9979", AMINE: "0.0021"},
    ),
    # Simple carboxylic acid (no tautomerisation, anion is preffered)
    TautomerCase(
        smiles=CARBOXYYLIC_ACID,
        canonical_tautomer=CARBOXYYLIC_ACID,
        major_tautomer=CARBOXYYLIC_ACID,
        major_tautomer_protomer=CARBOXYYLIC_ACID_ANION,
        major_protomer=CARBOXYYLIC_ACID_ANION,
        tautomer_distribution={CARBOXYYLIC_ACID: "1.0000"},
        protomer_distribution={CARBOXYYLIC_ACID_ANION: "0.9986", CARBOXYYLIC_ACID: "0.0014"},
        tautomer_protomer_distribution={
            CARBOXYYLIC_ACID_ANION: "0.9986",
            CARBOXYYLIC_ACID: "0.0014",
        },
    ),
    # Keto-enol tautomerisation
    # Chemaxon chooses keto form
    TautomerCase(
        smiles=ENOL,
        canonical_tautomer=KETO,
        major_tautomer=KETO,
        tautomer_distribution={KETO: "1.0000"},
        protomer_distribution={ENOL: "0.9999", ENOL_ANION: "0.0001"},
    ),
    TautomerCase(
        smiles=KETO,
        canonical_tautomer=KETO,
        major_tautomer=KETO,
        tautomer_distribution={KETO: "1.0000"},
    ),
    # Acetylacetone, example of enol being prefered
    TautomerCase(
        smiles=DIENOL,
        canonical_tautomer=DIENOL,
        major_tautomer=DIENOL,
        tautomer_distribution={DIENOL: "0.9543", DIKETO: "0.0457"},
        protomer_distribution={DIENOL: "0.9968", DIENOL_OXYGEN_ANION: "0.0032"},
        tautomer_protomer_distribution={
            DIENOL: "0.9486",
            DIKETO: "0.0454",
            DIENOL_OXYGEN_ANION: "0.0031",
            DIKETO_CARBON_ANION: "0.0029",
        },
    ),
    TautomerCase(
        smiles=DIKETO,
        canonical_tautomer=DIENOL,
        major_tautomer=DIENOL,
        tautomer_distribution={DIENOL: "0.9543", DIKETO: "0.0457"},
        protomer_distribution={DIKETO: "0.9392", DIKETO_CARBON_ANION: "0.0608"},
        tautomer_protomer_distribution={
            DIENOL: "0.9486",
            DIKETO: "0.0454",
            DIENOL_OXYGEN_ANION: "0.0031",
            DIKETO_CARBON_ANION: "0.0029",
        },
    ),
    # Lactam <=> Lactim
    # Chemaxon does not interconvert lactams and lactims
    TautomerCase(
        smiles=LACTAM,
        canonical_tautomer=LACTAM,
        major_tautomer=LACTAM,
        tautomer_distribution={LACTAM: "1.0000"},
    ),
    TautomerCase(
        smiles=LACTIM,
        canonical_tautomer=LACTIM,
        major_tautomer=LACTIM,
        tautomer_distribution={LACTIM: "1.0000"},
    ),
    # Amide <=> Imidic Acid
    # Chemaxon prefers amides
    TautomerCase(
        smiles=AMIDE,
        canonical_tautomer=AMIDE,
        major_tautomer=AMIDE,
        tautomer_distribution={AMIDE: "1.0000"},
    ),
    TautomerCase(
        smiles=IMIDIC_ACID,
        canonical_tautomer=AMIDE,
        major_tautomer=AMIDE,
        major_protomer=IMIDIC_ACID_ANION,
        tautomer_distribution={AMIDE: "1.0000"},
        protomer_distribution={
            IMIDIC_ACID_ANION: "0.5300",
            IMIDIC_ACID_ZWITTERION: "0.3426",
            IMIDIC_ACID: "0.1207",
            IMIDIC_ACID_CATION: "0.0068",
        },
    ),
    # 2-hydroxylpyridine and 2-pyridone
    # Chemaxon prefers 2-pyridone
    TautomerCase(
        smiles=TWOHYDROXYPYRIDINE,
        canonical_tautomer=TWOPYRIDONE,
        major_tautomer=TWOPYRIDONE,
        tautomer_distribution={TWOPYRIDONE: "1.0000"},
        protomer_distribution={TWOHYDROXYPYRIDINE: "0.9950", TWOHYDROXYPYRIDINE_ANION: "0.0050"},
        tautomer_protomer_distribution={TWOPYRIDONE: "0.9999", TWOPYRIDONE_ANION: "0.0001"},
    ),
    TautomerCase(
        smiles=TWOPYRIDONE,
        canonical_tautomer=TWOPYRIDONE,
        major_tautomer=TWOPYRIDONE,
        tautomer_distribution={TWOPYRIDONE: "1.0000"},
        protomer_distribution={TWOPYRIDONE: "0.9999", TWOPYRIDONE_ANION: "0.0001"},
        tautomer_protomer_distribution={TWOPYRIDONE: "0.9999", TWOPYRIDONE_ANION: "0.0001"},
    ),
    # Imine <=> Enamine
    # Same canonical tautomer, different major tautomer
    TautomerCase(
        smiles=IMINE,
        canonical_tautomer=IMINE,
        major_tautomer=IMINE,
        major_protomer=IMINE_CATION,
        major_tautomer_protomer=IMINE_CATION,
        tautomer_distribution={IMINE: "1.0000"},
        protomer_distribution={IMINE_CATION: "0.9999", IMINE: "0.0001"},
        tautomer_protomer_distribution={IMINE_CATION: "0.9999", IMINE: "0.0001"},
    ),
    TautomerCase(
        smiles=ENAMINE,
        canonical_tautomer=IMINE,
        major_tautomer=ENAMINE,
        tautomer_distribution={ENAMINE: "1.0000"},
        protomer_distribution={ENAMINE: "0.6911", ENAMINE_CATION: "0.3089"},
        tautomer_protomer_distribution={ENAMINE: "0.6911", ENAMINE_CATION: "0.3089"},
    ),
    # Cyanamide <=> Carbodiimide
    # Chemaxon does not interconvert
    TautomerCase(
        smiles=CYANAMIDE,
        canonical_tautomer=CYANAMIDE,
        major_tautomer=CYANAMIDE,
        tautomer_distribution={CYANAMIDE: "1.0000"},
    ),
    TautomerCase(
        smiles=CARBODIIMIDE,
        canonical_tautomer=CARBODIIMIDE,
        normal_canonical_tautomer=CARBODIIMIDE,
        major_tautomer=CARBODIIMIDE,
        # Seemingly a bug in Chemaxon, no protomers generated
        major_protomer="error",
        tautomer_distribution={CARBODIIMIDE: "1.0000"},
        # Seemingly a bug in Chemaxon, no protomers generated
        protomer_distribution={},
        tautomer_protomer_distribution={CARBODIIMIDE: "1.0000"},
    ),
    # Alanine
    TautomerCase(
        smiles=ALANINE,
        canonical_tautomer=ALANINE,
        major_tautomer=ALANINE,
        major_protomer=ALANINE_ZWITTERION,
        major_tautomer_protomer=ALANINE_ZWITTERION,
        tautomer_distribution={ALANINE: "1.0000"},
        protomer_distribution={ALANINE_ZWITTERION: "0.9917", ALANINE_ANION: "0.0083"},
        tautomer_protomer_distribution={ALANINE_ZWITTERION: "0.9917", ALANINE_ANION: "0.0083"},
    ),
    # Azide
    TautomerCase(
        smiles=AZIDE,
        canonical_tautomer=AZIDE_TRIPLE,
        major_tautomer=AZIDE,
        major_protomer=AZIDE_CATION,
        major_tautomer_protomer=AZIDE,
        tautomer_distribution={AZIDE: "1.0000"},
        protomer_distribution={AZIDE_CATION: "1.0000"},
    ),
    TautomerCase(
        smiles=AZIDE_TRIPLE,
        canonical_tautomer=AZIDE_TRIPLE,
        major_tautomer=AZIDE,
        major_protomer=AZIDE_TRIPLE,
        major_tautomer_protomer=AZIDE,
        tautomer_distribution={AZIDE: "1.0000"},
        protomer_distribution={AZIDE_TRIPLE: "0.9987", AZIDE_TRIPLE_CATION: "0.0013"},
    ),
    # 1,2,4-triazole
    TautomerCase(
        smiles=TRIAZOLE_1H,
        canonical_tautomer=TRIAZOLE_1H,
        major_tautomer=TRIAZOLE_1H,
        tautomer_distribution={TRIAZOLE_1H: "0.7572", TRIAZOLE_4H: "0.2428"},
        protomer_distribution={TRIAZOLE_1H: "0.9968", TRIAZOLE_1H_ANION: "0.0032"},
        tautomer_protomer_distribution={
            TRIAZOLE_1H: "0.7552",
            TRIAZOLE_4H: "0.2422",
            TRIAZOLE_1H_ANION: "0.0024",
            TRIAZOLE_4H_ANION: "0.0002",
        },
    ),
    TautomerCase(
        smiles=TRIAZOLE_4H,
        canonical_tautomer=TRIAZOLE_1H,
        major_tautomer=TRIAZOLE_1H,
        tautomer_distribution={TRIAZOLE_1H: "0.7572", TRIAZOLE_4H: "0.2428"},
        protomer_distribution={TRIAZOLE_4H: "0.9990", TRIAZOLE_4H_ANION: "0.0010"},
        tautomer_protomer_distribution={
            TRIAZOLE_1H: "0.7552",
            TRIAZOLE_4H: "0.2422",
            TRIAZOLE_1H_ANION: "0.0024",
            TRIAZOLE_4H_ANION: "0.0002",
        },
    ),
    TautomerCase(
        smiles=GUANIDINE1,
        canonical_tautomer=GUANIDINE1,
        major_tautomer=GUANIDINE1,
        major_protomer=GUANIDINE1_CATION,
        major_tautomer_protomer=GUANIDINE1_CATION,
        tautomer_distribution={GUANIDINE1: "0.5368", GUANIDINE2: "0.2917", GUANIDINE3: "0.1715"},
        protomer_distribution={GUANIDINE1_CATION: "1.0000"},
        tautomer_protomer_distribution={
            GUANIDINE1_CATION: "0.5369",
            GUANIDINE2_CATION: "0.2916",
            GUANIDINE3_CATION: "0.1715",
        },
    ),
    TautomerCase(
        smiles=GUANIDINE2,
        canonical_tautomer=GUANIDINE1,
        normal_canonical_tautomer=GUANIDINE2,
        major_tautomer=GUANIDINE1,
        major_protomer=GUANIDINE2_CATION,
        major_tautomer_protomer=GUANIDINE1_CATION,
        tautomer_distribution={GUANIDINE1: "0.5368", GUANIDINE2: "0.2917", GUANIDINE3: "0.1715"},
        protomer_distribution={GUANIDINE2_CATION: "1.0000"},
        tautomer_protomer_distribution={
            GUANIDINE1_CATION: "0.5369",
            GUANIDINE2_CATION: "0.2916",
            GUANIDINE3_CATION: "0.1715",
        },
    ),
    TautomerCase(
        smiles=GUANIDINE3,
        canonical_tautomer=GUANIDINE1,
        normal_canonical_tautomer=GUANIDINE3,
        major_tautomer=GUANIDINE1,
        major_protomer=GUANIDINE3_CATION,
        major_tautomer_protomer=GUANIDINE1_CATION,
        tautomer_distribution={GUANIDINE1: "0.5368", GUANIDINE2: "0.2917", GUANIDINE3: "0.1715"},
        protomer_distribution={GUANIDINE3_CATION: "1.0000"},
        tautomer_protomer_distribution={
            GUANIDINE1_CATION: "0.5369",
            GUANIDINE2_CATION: "0.2916",
            GUANIDINE3_CATION: "0.1715",
        },
    ),
    # Nitroso <=> Oxime tautomers
    # Oxime is the canonical tautomer, but cannot interconvert using normal mode
    TautomerCase(
        smiles=NITROSO,
        canonical_tautomer=OXIME,
        major_tautomer=NITROSO,
        tautomer_distribution={NITROSO: "1.0000"},
    ),
    TautomerCase(
        smiles=OXIME,
        canonical_tautomer=OXIME,
        major_tautomer=OXIME,
        tautomer_distribution={OXIME: "1.0000"},
        protomer_distribution={OXIME: "0.9998", OXIME_ANION: "0.0001", OXIME_CATION: "0.0001"},
        tautomer_protomer_distribution={
            OXIME: "0.9998",
            OXIME_CATION: "0.0001",
            OXIME_ANION: "0.0001",
        },
    ),
]


@pytest.mark.parametrize(
    ["smiles", "canonical_smiles"], [(case.smiles, case.canonical_tautomer) for case in CASES]
)
def test_get_canonical_tautomer(smiles, canonical_smiles):
    assert get_canonical_tautomer(smiles) == canonical_smiles


def _clean_oemol_from_smiles(smiles):
    oemol = oemol_from_smiles(smiles)
    make_hydrogens_implicit(oemol, remove_isotopic_hydrogens=True)
    oemol.Sweep()
    return oemol


def _assert_tautomer_equal(oemol1, oemol2):
    assert_mol_equal(
        oemol1,
        oemol2,
        atom_properties=AllAtomProperties
        & ~AtomProperty.Degree
        & ~AtomProperty.Valence
        & ~AtomProperty.FormalCharge
        & ~AtomProperty.ImplicitHCount
        & ~AtomProperty.ExplicitHCount,
        bond_properties=AllBondProperties
        & ~BondProperty.Type
        & ~BondProperty.Order
        & ~BondProperty.IsRotor,
    )


@pytest.mark.parametrize(
    ["smiles", "canonical_smiles"], [(case.smiles, case.canonical_tautomer) for case in CASES]
)
def test_get_canonical_tautomer_oemol(smiles, canonical_smiles):
    oemol = _clean_oemol_from_smiles(smiles)
    canonical_tautomer_oemol = get_canonical_tautomer(oemol)
    assert smiles_from_oemol(canonical_tautomer_oemol) == canonical_smiles
    _assert_tautomer_equal(oemol, canonical_tautomer_oemol)


@pytest.mark.parametrize(
    ["smiles", "canonical_smiles"],
    [
        (
            case.smiles,
            case.normal_canonical_tautomer
            if case.normal_canonical_tautomer
            else case.canonical_tautomer,
        )
        for case in CASES
    ],
)
def test_get_normal_canonical_tautomer_oemol(smiles, canonical_smiles):
    oemol = _clean_oemol_from_smiles(smiles)
    canonical_tautomer_oemol = get_canonical_tautomer(oemol, normal=True)
    assert smiles_from_oemol(canonical_tautomer_oemol) == canonical_smiles
    _assert_tautomer_equal(oemol, canonical_tautomer_oemol)


@pytest.mark.parametrize(
    ["smiles", "major_smiles"], [(case.smiles, case.major_tautomer) for case in CASES]
)
def test_get_major_tautomer(smiles, major_smiles):
    assert get_major_tautomer(smiles) == major_smiles


@pytest.mark.parametrize(
    ["smiles", "major_smiles"], [(case.smiles, case.major_tautomer) for case in CASES]
)
def test_get_major_tautomer_oemol(smiles, major_smiles):
    oemol = _clean_oemol_from_smiles(smiles)
    major_oemol = get_major_tautomer(oemol)
    assert smiles_from_oemol(major_oemol) == major_smiles
    _assert_tautomer_equal(oemol, major_oemol)


@pytest.mark.parametrize(
    ["smiles", "major_smiles"],
    [(case.smiles, case.major_protomer if case.major_protomer else case.smiles) for case in CASES],
)
def test_get_major_protomer(smiles, major_smiles):
    if major_smiles == "error":
        return
    assert get_major_protomer(smiles) == major_smiles


@pytest.mark.parametrize(
    ["smiles", "major_smiles"],
    [(case.smiles, case.major_protomer if case.major_protomer else case.smiles) for case in CASES],
)
def test_get_major_protomer_oemol(smiles, major_smiles):
    if major_smiles == "error":
        return
    oemol = _clean_oemol_from_smiles(smiles)
    major_oemol = get_major_protomer(oemol)
    assert smiles_from_oemol(major_oemol) == major_smiles
    _assert_tautomer_equal(oemol, major_oemol)


@pytest.mark.parametrize(
    ["smiles", "major_tautomer"],
    [
        (
            case.smiles,
            case.major_tautomer_protomer if case.major_tautomer_protomer else case.major_tautomer,
        )
        for case in CASES
    ],
)
def test_get_major_tautomer_and_protomer(smiles, major_tautomer):
    assert get_major_tautomer_and_protomer(smiles) == major_tautomer


@pytest.mark.parametrize(
    ["smiles", "major_tautomer"],
    [
        (
            case.smiles,
            case.major_tautomer_protomer if case.major_tautomer_protomer else case.major_tautomer,
        )
        for case in CASES
    ],
)
def test_get_major_tautomer_and_protomer_oemol(smiles, major_tautomer):
    oemol = _clean_oemol_from_smiles(smiles)
    major_oemol = get_major_tautomer_and_protomer(oemol)
    assert smiles_from_oemol(major_oemol) == major_tautomer
    _assert_tautomer_equal(oemol, major_oemol)


@pytest.mark.parametrize(
    ["smiles", "distribution"], [(case.smiles, case.tautomer_distribution) for case in CASES]
)
def test_get_tautomer_distribution(smiles, distribution):
    assert {
        key: f"{value:.4f}" for key, value in get_tautomer_distribution(smiles).items()
    } == distribution


@pytest.mark.parametrize(
    ["smiles", "distribution"], [(case.smiles, case.tautomer_distribution) for case in CASES]
)
def test_get_tautomer_distribution_oemol(smiles, distribution):
    oemol = _clean_oemol_from_smiles(smiles)
    distr = get_tautomer_distribution(oemol)
    assert {smiles_from_oemol(key): f"{value:.4f}" for key, value in distr.items()} == distribution
    for tautomer in distr:
        _assert_tautomer_equal(oemol, tautomer)


@pytest.mark.parametrize(
    ["smiles", "distribution"],
    [
        (smiles, case.tautomer_distribution)
        for case in CASES
        for smiles in case.tautomer_distribution
    ],
)
def test_get_tautomer_distribution_independent(smiles, distribution):
    assert {
        key: f"{value:.4f}" for key, value in get_tautomer_distribution(smiles).items()
    } == distribution


@pytest.mark.parametrize(
    ["smiles", "distribution"],
    [
        (
            case.smiles,
            case.protomer_distribution
            if case.protomer_distribution is not None
            else case.tautomer_distribution,
        )
        for case in CASES
    ],
)
def test_get_protomer_distribution_oemol(smiles, distribution):
    oemol = _clean_oemol_from_smiles(smiles)
    distr = get_protomer_distribution(oemol)
    distr_with_smiles = {smiles_from_oemol(key): f"{value:.4f}" for key, value in distr.items()}
    assert list(distr_with_smiles.items()) == list(distribution.items())
    for tautomer in distr:
        _assert_tautomer_equal(oemol, tautomer)


@pytest.mark.parametrize(
    ["smiles", "distribution"],
    [
        (smiles, case.protomer_distribution)
        for case in CASES
        if case.protomer_distribution is not None
        for smiles in case.protomer_distribution
    ],
)
def test_get_protomer_distribution_independent(smiles, distribution):
    distr = {key: f"{value:.4f}" for key, value in get_protomer_distribution(smiles).items()}
    assert list(distr.items()) == list(distribution.items())


@pytest.mark.parametrize(
    ["smiles", "distribution"],
    [
        (
            case.smiles,
            case.tautomer_protomer_distribution
            if case.tautomer_protomer_distribution is not None
            else case.tautomer_distribution,
        )
        for case in CASES
    ],
)
def test_get_tautomer_and_protomer_distribution(smiles, distribution):
    assert {
        key: f"{value:.4f}" for key, value in get_tautomer_and_protomer_distribution(smiles).items()
    } == distribution


@pytest.mark.parametrize(
    ["smiles", "distribution"],
    [
        (
            case.smiles,
            case.tautomer_protomer_distribution
            if case.tautomer_protomer_distribution is not None
            else case.tautomer_distribution,
        )
        for case in CASES
    ],
)
def test_get_tautomer_and_protomer_distribution_oemol(smiles, distribution):
    oemol = _clean_oemol_from_smiles(smiles)
    distr = get_tautomer_and_protomer_distribution(oemol)
    assert {smiles_from_oemol(key): f"{value:.4f}" for key, value in distr.items()} == distribution
    for tautomer in distr:
        _assert_tautomer_equal(oemol, tautomer)


@pytest.mark.parametrize(
    ["smiles", "distribution"],
    [
        (smiles, case.tautomer_protomer_distribution)
        for case in CASES
        if case.tautomer_protomer_distribution is not None
        for smiles in case.tautomer_protomer_distribution
    ],
)
def test_get_tautomer_and_protomer_distribution_independent(smiles, distribution):
    assert {
        key: f"{value:.4f}" for key, value in get_tautomer_and_protomer_distribution(smiles).items()
    } == distribution
