import math
from collections.abc import Callable

from .atomic import MolecularWeight, NumberHeavyAtoms
from .bonds_rings import NumberAromaticRings
from .hydrogen_bonds import NumberHydrogenAcceptors, NumberHydrogenDonors
from .pka import PKaConjugateAcid
from .property import MolecularProperty
from .tpsa import TwoDimensionalTPSA


def linear_polynomial(a: float, b: float) -> Callable[[float], float]:
    return lambda x: a * x + b


def quadratic_polynomial(
    a: float, b: float, c: float, d: float, e: float
) -> Callable[[float], float]:
    return lambda x: a * x**4 + b * x**3 + c * x**2 + d * x + e


def cubic_polynomial(a: float, b: float, c: float, d: float) -> Callable[[float], float]:
    return lambda x: a * x**3 + b * x**2 + c * x + d


def _aromatic_ring_score(num_aromatic_rings: int) -> float:
    if num_aromatic_rings == 0:
        return 0.336376
    if num_aromatic_rings == 1:
        return 0.816016
    if num_aromatic_rings == 2:
        return 1.0
    if num_aromatic_rings == 3:
        return 0.691115
    if num_aromatic_rings == 4:
        return 0.199399
    return 0.0


def _heavy_atom_score(num_heavy_atoms: int) -> float:
    if 5 < num_heavy_atoms <= 45:
        return (
            cubic_polynomial(a=0.0000443, b=-0.004556, c=0.12775, d=-0.463)(num_heavy_atoms)
            / 0.624231
        )
    return 0.0


def _mwhbn_score(mwhbn_number: float) -> float:
    if 0.05 <= mwhbn_number <= 0.45:
        return cubic_polynomial(a=26.733, b=-31.495, c=9.5202, d=-0.1358)(mwhbn_number) / 0.72258
    return 0.0


def _get_mwhbn(molecular_weight: float, acceptor_count: int, donor_count: int) -> float:
    return (acceptor_count + donor_count) / math.sqrt(molecular_weight)


def _tpsa_score(tpsa: float) -> float:
    if 0 < tpsa <= 120:
        return linear_polynomial(a=-0.0067, b=0.9598)(tpsa) / 0.9598
    return 0.0


def _pka_score(pka: float) -> float:
    if 3.0 < pka <= 11.0:
        return (
            quadratic_polynomial(a=0.00045068, b=-0.016331, c=0.18618, d=-0.71043, e=0.8579)(pka)
            / 0.597488
        )
    return 0.0


def blood_brain_barrier_score(
    *,
    num_aromatic_rings: int,
    num_heavy_atoms: int,
    molecular_weight: float,
    hbond_acceptor_count: int,
    hbond_donor_count: int,
    tpsa: float,
    basic_pka: float | None,
) -> float:
    """
    Calculate the Blood-Brain Barrier (BBB) score.

    This is based on the paper:
    Gupta, Mayuri, et al. "The blood-brain barrier (BBB) score." Journal of medicinal chemistry 62.21 (2019): 9824-9836.
    https://pubs.acs.org/doi/10.1021/acs.jmedchem.9b01220

    BBB adds up 5 separate scores of different weight, resulting in the total score between 0 (worst) and 6 (best).

    Args:
        num_aromatic_rings: number of aromatic rings in molecule
        num_heavy_atoms: heavy atom count for molecule
        molecular_weight: molecular weight
        hbond_acceptor_count: number of acceptor atoms in molecule
        hbond_donor_count: number of donor atoms in molecule
        tpsa: topological polar surface area
        basic_pka: most basic center acid dissociation constant

    Returns:
        Blood-Brain Barrier score.
    """
    mwhbn = _get_mwhbn(molecular_weight, hbond_acceptor_count, hbond_donor_count)
    bbb_score = (
        1.0 * _aromatic_ring_score(num_aromatic_rings)
        + 1.0 * _heavy_atom_score(num_heavy_atoms)
        + 1.5 * _mwhbn_score(mwhbn)
        + 2.0 * _tpsa_score(tpsa)
        + 0.5 * (_pka_score(basic_pka) if basic_pka is not None else 1)
    )
    return bbb_score


BloodBrainBarrierScore = MolecularProperty(
    "bbb_score",
    composite_func=blood_brain_barrier_score,
    composite_properties={
        "num_aromatic_rings": NumberAromaticRings,
        "num_heavy_atoms": NumberHeavyAtoms,
        "molecular_weight": MolecularWeight,
        "hbond_acceptor_count": NumberHydrogenAcceptors,
        "hbond_donor_count": NumberHydrogenDonors,
        "tpsa": TwoDimensionalTPSA,
        "basic_pka": PKaConjugateAcid,
    },
)
