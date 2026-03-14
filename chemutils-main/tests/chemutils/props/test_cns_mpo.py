import pytest

from chemutils.props.cns_mpo import (
    _basic_pka_score,
    _clip_score,
    _hbond_donor_count_score,
    _logd_score,
    _logp_score,
    _mw_score,
    _tpsa_score,
    cns_mpo,
)


@pytest.mark.parametrize(
    ["unconstrained_score", "expected_constrained_score"],
    [
        (-0.1, 0.0),
        (0.0, 0.0),
        (0.12, 0.12),
        (0.98, 0.98),
        (1.0, 1.0),
        (1.1, 1.0),
    ],
)
def test_clip_score(unconstrained_score, expected_constrained_score):
    assert _clip_score(unconstrained_score) == expected_constrained_score


@pytest.mark.parametrize(
    ["logp", "expected_logp_score"],
    [
        (2.9, 1.0),
        (3.0, 1.0),
        (3.5, 0.75),
        (4.0, 0.5),
        (4.5, 0.25),
        (5.0, 0.0),
        (5.2, 0.0),
    ],
)
def test_logp_score(logp, expected_logp_score):
    assert _logp_score(logp) == expected_logp_score


@pytest.mark.parametrize(
    ["logd", "expected_logd_score"],
    [
        (1.9, 1.0),
        (2.0, 1.0),
        (2.5, 0.75),
        (3.0, 0.5),
        (3.5, 0.25),
        (4.0, 0.0),
        (4.1, 0.0),
    ],
)
def test_logd_score(logd, expected_logd_score):
    assert _logd_score(logd) == expected_logd_score


@pytest.mark.parametrize(
    ["mw", "expected_mw_score"],
    [
        (359.5, 1.0),
        (360.0, 1.0),
        (395.0, 0.75),
        (430.0, 0.5),
        (465.0, 0.25),
        (500.0, 0.0),
        (501.1, 0.0),
    ],
)
def test_mw_score(mw, expected_mw_score):
    assert _mw_score(mw) == expected_mw_score


@pytest.mark.parametrize(
    ["tpsa", "expected_tpsa_score"],
    [
        (19.9, 0.0),
        (20.0, 0.0),
        (25.0, 0.25),
        (30.0, 0.5),
        (35.0, 0.75),
        (40.0, 1.0),
        (65.0, 1.0),
        (90.0, 1.0),
        (97.5, 0.75),
        (105.0, 0.5),
        (112.5, 0.25),
        (120.0, 0.0),
        (120.1, 0.0),
    ],
)
def test_tpsa_score(tpsa, expected_tpsa_score):
    assert _tpsa_score(tpsa) == expected_tpsa_score


@pytest.mark.parametrize(
    ["hbond_donor_count", "expected_hbond_donor_count_score"],
    [
        (0, 1.0),
        (1, 0.75),
        (2, 0.5),
        (3, 0.25),
        (4, 0.0),
        (5, 0.0),
    ],
)
def test_hbond_donor_count_score(hbond_donor_count, expected_hbond_donor_count_score):
    assert _hbond_donor_count_score(hbond_donor_count) == expected_hbond_donor_count_score


@pytest.mark.parametrize(
    ["basic_pka", "expected_basic_pka_score"],
    [
        (7.9, 1.0),
        (8.0, 1.0),
        (8.5, 0.75),
        (9.0, 0.5),
        (9.5, 0.25),
        (10.0, 0.0),
        (10.1, 0.0),
    ],
)
def test_basic_pka_score(basic_pka, expected_basic_pka_score):
    assert _basic_pka_score(basic_pka) == expected_basic_pka_score


@pytest.mark.parametrize(
    ["logp", "logd", "mw", "tpsa", "hbond_donor_count", "basic_pka", "expected_cns_mpo_score"],
    [
        (4.0, 3.0, 430.0, 112.5, 2, 12.0, 2.25),
        (3.0, 3.5, 465.0, 120.0, 1, 9.0, 2.75),
        (6.0, 5.0, 510.0, 125.0, 4, 11.0, 0.0),
    ],
)
def test_cns_mpo(logp, logd, mw, tpsa, hbond_donor_count, basic_pka, expected_cns_mpo_score):
    assert (
        cns_mpo(
            logp=logp,
            logd=logd,
            mw=mw,
            tpsa=tpsa,
            hbond_donor_count=hbond_donor_count,
            basic_pka=basic_pka,
        )
        == expected_cns_mpo_score
    )
