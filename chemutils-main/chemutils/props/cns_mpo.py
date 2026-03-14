from .atomic import MolecularWeight
from .hydrogen_bonds import NumberHydrogenDonors
from .logd import LogD
from .logp import LogP
from .pka import PKaConjugateAcid
from .property import MolecularProperty
from .tpsa import TwoDimensionalTPSA


def _clip_score(score: float) -> float:
    """
    Constrain score to the range [0, 1]. Values below 0 are set 0, and values above 1 are set to 1.

    Args:
        score: Arbitrary float score

    Returns:
        Score constrained to the range [0, 1]
    """
    return 1.0 if score > 1.0 else 0.0 if score < 0.0 else score


def _logp_score(logp: float) -> float:
    """
    Calculate the partition coefficient (logP) score.

    From https://squonk.it/docs/cells/CNS%20MPO%20(CXN)/
        if value < 3 then 1
        ramp from 1 to 0 over range of 3 to 5
        if value >= 5 then 0.

    Args:
        logp: (calculated) partition coefficient

    Returns:
        logP score
    """
    return _clip_score((5.0 - logp) * 0.5)


def _logd_score(logd: float) -> float:
    """
    Calculate distribution coefficient (logD) score.

    From https://squonk.it/docs/cells/CNS%20MPO%20(CXN)/
        if value < 2 then 1
        ramp from 1 to 0 over range of 2 to 4
        if value >= 4 then 0.

    Args:
        logd: (calculated) distribution coefficient (at pH = 7.4)

    Returns:
        logD score
    """
    return _clip_score((4.0 - logd) * 0.5)


def _mw_score(mw: float) -> float:
    """
    Calculate the molecular weight score.

    From https://squonk.it/docs/cells/CNS%20MPO%20(CXN)/
        if value < 360 then 1
        ramp from 1 to 0 over range of 360 to 500
        if value >= 500 then 0.

    Args:
        mw: molecular weight

    Returns:
        molecular weight score
    """
    return _clip_score((500.0 - mw) / (500.0 - 360.0))


def _tpsa_score(tpsa: float) -> float:
    """
    Calculate the topological polar surface area (TPSA) score.

    From https://squonk.it/docs/cells/CNS%20MPO%20(CXN)/
        if value < 20 then 0
        ramp from 0 to 1 over range of 20 to 40
        if value between 40 and 90 then 1
        ramp from 1 to 0 over range of 90 to 120
        if value >= 120 then 0.

    Args:
        tpsa: topological polar surface area

    Returns:
        TPSA score
    """
    if tpsa < 65:
        return _clip_score((tpsa - 20.0) / (40.0 - 20.0))
    else:
        return _clip_score((120.0 - tpsa) / (120.0 - 90.0))


def _hbond_donor_count_score(hbond_donor_count: int) -> float:
    """
    Calculate the H-bond donor count score.

    From https://squonk.it/docs/cells/CNS%20MPO%20(CXN)/
        if value = 0 then 1
        if value = 1 then 0.75
        if value = 2 then 0.5
        if value = 3 then 0.25
        if value >= 4 then 0.

    Args:
        hbond_donor_count: number of hydrogen bond donors

    Returns:
        H-bond donor count score
    """
    return _clip_score((4.0 - hbond_donor_count) * 0.25)


def _basic_pka_score(basic_pka: float | None) -> float:
    """
    Calculate the most basic center acid dissociation constant (pKa) score.

    From https://squonk.it/docs/cells/CNS%20MPO%20(CXN)/
        if value < 8 then 1
        ramp from 1 to 0 over range of 8 to 10
        if value >= 10 then 0.

    Args:
        basic_pka: most basic center acid dissociation constant

    Returns:
        most basic center pKa score
    """
    if basic_pka is None:
        return 0
    return _clip_score((10.0 - basic_pka) * 0.5)


def cns_mpo(
    *,
    logp: float,
    logd: float,
    mw: float,
    tpsa: float,
    hbond_donor_count: int,
    basic_pka: float | None,
) -> float:
    """
    Calculate the Multi-Parameter Optimisation (MPO) score for drugs targeted at the Central Nervous System (CNS).

    Based on the paper:
    "Moving beyond Rules: The Development of a Central Nervous System Multiparameter Optimization (CNS MPO) Approach To Enable Alignment of Druglike Properties"
    by Travis T. Wager, Xinjun Hou, Patrick R. Verhoest, and Anabella Villalobos
    https://pubs.acs.org/doi/10.1021/cn100008c

    CNS MPO adds up 6 separate scores each bounded between 0 (worst) and 1 (best), resulting in the total score between 0 (worst) and 6 (best).

    Args:
        logp: (calculated) partition coefficient
        logd: (calculated) distribution coefficient (at pH = 7.4)
        mw: molecular weight
        tpsa: topological polar surface area
        hbond_donor_count: number of hydrogen bond donors
        basic_pka: most basic center acid dissociation constant

    Returns:
        CNS MPO score
    """
    cns_mpo_score = (
        _logp_score(logp)
        + _logd_score(logd)
        + _mw_score(mw)
        + _tpsa_score(tpsa)
        + _hbond_donor_count_score(hbond_donor_count)
        + _basic_pka_score(basic_pka)
    )
    return cns_mpo_score


CentralNervousSystemMPO = MolecularProperty(
    "cns_mpo",
    composite_func=cns_mpo,
    composite_properties={
        "logp": LogP,
        "logd": LogD,
        "mw": MolecularWeight,
        "tpsa": TwoDimensionalTPSA,
        "hbond_donor_count": NumberHydrogenDonors,
        "basic_pka": PKaConjugateAcid,
    },
)
