import pytest

from chemutils.props.blood_brain_barrier import (
    _aromatic_ring_score,
    _get_mwhbn,
    blood_brain_barrier_score,
)


@pytest.mark.parametrize(
    ["num_aro_rings", "expected_aro_ring_score"],
    [
        (3, 0.69),
    ],
)
def test_aromatic_ring_score(num_aro_rings, expected_aro_ring_score):
    assert _aromatic_ring_score(num_aromatic_rings=num_aro_rings) == pytest.approx(
        expected_aro_ring_score, abs=0.01
    )


@pytest.mark.parametrize(
    [
        "num_aromatic_rings",
        "num_heavy_atoms",
        "molecular_weight",
        "hbond_acceptor_count",
        "hbond_donor_count",
        "tpsa",
        "basic_pka",
        "mwhbn_score",
        "bbb_score",
    ],
    [
        (2, 13, 180.17, 3, 1, 69.300003, 7.82, 0.3, 4.62),  # Aminophylline
        (3, 26, 381.37, 3, 1, 77.980003, None, 0.2, 4.5),  # Celecoxib
        (1, 13, 202.63, 3, 2, 49.689999, None, 0.35, 4.46),  # Chlorphenesin
        (1, 14, 230.09, 3, 2, 36.419998, 8.16, 0.33, 4.8),  # Clonidine
        (1, 22, 303.36, 3, 0, 55.84, 8.85, 0.17, 4.99),  # Cocaine
        (3, 32, 426.56, 3, 1, 55.560001, 10.11, 0.19, 4.52),  # Darifenacin
        (0, 10, 168.04, 1, 0, 9.2299995, None, 0.08, 4.23),  # Desflurane
        (0, 28, 392.47, 5, 3, 94.830002, None, 0.4, 3.03),  # Dexamethasone
        (2, 15, 200.29, 1, 1, 28.68, 7.07, 0.14, 5.3),  # Dexmedetomidine
        (2, 18, 250.2, 3, 2, 57.529999, 2.69, 0.32, 4.37),  # Diflunisal
        (0, 16, 296.52, 0, 0, 6.48, None, 0, 3.7),  # Disulfiram
        (1, 15, 213.19, 6, 5, 124.01, 8.72, 0.75, 2.24),  # Droxidopa
        (1, 12, 160.18, 2, 0, 32.669998, 7.83, 0.16, 5.05),  # Edaravone
        (2, 22, 312.48, 2, 0, 6.48, 9.6, 0.11, 5.57),  # Ethopropazine
        (2, 21, 287.36, 3, 2, 62.32, 4.73, 0.29, 4.54),  # Etodolac
        (2, 18, 242.27, 2, 1, 46.529999, 3.96, 0.19, 4.89),  # Fenoprofen
        (2, 30, 411.59, 3, 1, 49.77, 10.64, 0.2, 4.95),  # Fesoterodine
        (2, 30, 434.52, 3, 1, 26.709999, 8.03, 0.19, 5.35),  # Flupenthixol
        (0, 24, 338.44, 5, 1, 72.830002, None, 0.33, 3.9),  # Idebenone
        (1, 23, 431.29, 2, 0, 29.540001, 9.46, 0.1, 4.92),  # Ioflupane
        (2, 22, 411.2, 3, 0, 64.43, None, 0.15, 4.97),  # Iomazenil
        (2, 19, 254.29, 3, 1, 54.369999, 3.88, 0.25, 4.72),  # Ketoprofen
        (2, 19, 255.27, 3, 1, 59.299999, 3.84, 0.25, 4.65),  # Ketorolac
        (1, 16, 259.13, 3, 1, 33.619999, 9.27, 0.25, 5.24),  # Lofexidine
        (2, 18, 241.29, 3, 2, 49.330002, 3.89, 0.32, 4.5),  # Mefenamic acid
        (0, 27, 374.48, 5, 3, 94.830002, None, 0.41, 3.01),  # Methylprednisolone
        (0, 15, 219.28, 5, 4, 84.160004, 8.49, 0.61, 2.58),  # Miglustat
        (3, 28, 396.51, 5, 4, 100.27, 9.62, 0.45, 2.6),  # Mirabegron
        (0, 27, 382.54, 4, 2, 83.830002, None, 0.31, 3.76),  # Misoprostol
        (2, 17, 228.29, 2, 0, 26.299999, None, 0.13, 5.43),  # Nabumetone
        (2, 17, 230.26, 3, 1, 46.529999, 4.19, 0.26, 4.8),  # Naproxen
        (1, 16, 223.3, 1, 0, 29.540001, None, 0.07, 4.62),  # Neostigmine
        (3, 31, 484.39, 4, 0, 56.59, 8.12, 0.18, 4.57),  # Nicergoline
        (1, 30, 418.45, 6, 1, 117, None, 0.34, 3.48),  # Nimodipine
        (2, 19, 270.72, 2, 1, 41.459999, None, 0.18, 5.4),  # Nordazepam
        (1, 26, 357.49, 3, 1, 49.77, 8.77, 0.21, 5.01),  # Oxybutynin
        (0, 9, 132.16, 3, 0, 27.690001, None, 0.26, 4.42),  # Paraldehyde
        (1, 14, 189.21, 2, 0, 37.380001, None, 0.15, 5.06),  # Phensuximide
        (1, 13, 181.21, 1, 0, 33.419998, None, 0.07, 4.54),  # Pyridostigmine
        (2, 26, 352.43, 3, 1, 54.560001, 7.33, 0.21, 5.07),  # Rauwolfia Serpentina
        (4, 34, 477.57, 3, 0, 35.91, 8, 0.14, 4.09),  # Ritanserin
        (2, 27, 362.47, 2, 0, 32.779999, 8.88, 0.11, 5.04),  # Solifenacin
        (2, 25, 356.41, 3, 1, 54.369999, 4.09, 0.21, 4.72),  # Sulindac
        (1, 18, 290.35, 4, 1, 97.540001, None, 0.29, 4.24),  # Sulthiame
        (1, 29, 405.42, 6, 4, 170.59, 6.74, 0.5, 1.98),  # Taltirelin
        (2, 25, 375.55, 3, 1, 40.540001, 9.26, 0.21, 5.35),  # Tiagabine
        (2, 16, 253.71, 5, 2, 62.200001, 7.49, 0.44, 4.02),  # Tizanidine
        (2, 19, 257.29, 3, 1, 59.299999, 3.93, 0.25, 4.66),  # Tolmetin
        (2, 23, 311.47, 2, 1, 23.469999, 10.14, 0.17, 5.55),  # Tolterodine
        (2, 29, 392.52, 2, 1, 46.529999, None, 0.15, 5.04),  # Trospium
    ],
)
def test_blood_brain_barrier(
    num_aromatic_rings,
    num_heavy_atoms,
    molecular_weight,
    hbond_acceptor_count,
    hbond_donor_count,
    tpsa,
    basic_pka,
    mwhbn_score,
    bbb_score,
):
    assert _get_mwhbn(
        molecular_weight=molecular_weight,
        acceptor_count=hbond_acceptor_count,
        donor_count=hbond_donor_count,
    ) == pytest.approx(mwhbn_score, abs=0.01)
    assert blood_brain_barrier_score(
        num_aromatic_rings=num_aromatic_rings,
        num_heavy_atoms=num_heavy_atoms,
        molecular_weight=molecular_weight,
        hbond_acceptor_count=hbond_acceptor_count,
        hbond_donor_count=hbond_donor_count,
        tpsa=tpsa,
        basic_pka=basic_pka,
    ) == pytest.approx(bbb_score, abs=0.01)
