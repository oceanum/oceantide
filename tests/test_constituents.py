"""Tests for the constituent tables and nodal corrections."""

import numpy as np
import pytest

from oceantide.constituents import OMEGA, PERIODS, SHALLOW, V0U
from oceantide.core.utils import nodal

# Somewhere in the middle of a nodal cycle, so f and u are not near unity.
TIME = 48622.0 + 12000.0


def test_tables_cover_the_same_constituents():
    """A constituent in one table and not another is how silent errors start."""
    assert set(OMEGA) == set(V0U) == set(PERIODS)


@pytest.mark.parametrize("con", sorted(c for c in OMEGA if c != "Z0"))
def test_omega_agrees_with_period(con):
    """Angular speed and period must describe the same constituent."""
    assert OMEGA[con] == pytest.approx(2 * np.pi / (PERIODS[con] * 3600.0), rel=1e-5)


def test_ms4_is_the_sum_of_m2_and_s2():
    """Regression: MS4 was 2.811149e-04, which is not M2 + S2."""
    assert OMEGA["MS4"] == pytest.approx(OMEGA["M2"] + OMEGA["S2"], rel=1e-5)


@pytest.mark.parametrize("con,parents", sorted(SHALLOW.items()))
def test_compound_speeds_are_the_sum_of_their_parents(con, parents):
    expected = sum(OMEGA[p] * k for p, k in parents.items())
    assert OMEGA[con] == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize("con,parents", sorted(SHALLOW.items()))
def test_compound_equilibrium_arguments_follow_their_parents(con, parents):
    expected = sum(V0U[p] * k for p, k in parents.items()) % (2 * np.pi)
    got = V0U[con] % (2 * np.pi)
    diff = abs((got - expected + np.pi) % (2 * np.pi) - np.pi)
    assert diff < 1e-3


def test_every_constituent_has_a_nodal_entry():
    """No constituent may fall through to f=1, u=0, V0=0 unnoticed."""
    cons = sorted(OMEGA)
    pu, pf, v0u = nodal(TIME, cons)
    # Z0 is the mean, and the solar constituents genuinely have no nodal
    # modulation; everything else must be corrected.
    unmodulated = {"Z0", "SA", "SSA", "S1", "S2", "T2", "P1"}
    for con, f, u in zip(cons, pf, pu):
        if con in unmodulated:
            assert f == pytest.approx(1.0), con
        else:
            assert f != pytest.approx(1.0, abs=1e-6), f"{con} has no nodal factor"


def test_unknown_constituent_warns_and_is_dropped():
    with pytest.warns(UserWarning, match="No equilibrium argument"):
        pu, pf, v0u = nodal(TIME, ["M2", "NOPE"])
    assert pf[1] == 0.0
    assert v0u[1] == 0.0


@pytest.mark.parametrize(
    "con,parents",
    [("M4", {"M2": 2}), ("MS4", {"M2": 1, "S2": 1}), ("2MS6", {"M2": 2, "S2": 1})],
)
def test_compound_nodal_factors_are_products_of_parents(con, parents):
    names = sorted({con} | set(parents))
    pu, pf, _ = nodal(TIME, names)
    f = dict(zip(names, pf))
    u = dict(zip(names, pu))
    assert f[con] == pytest.approx(
        np.prod([f[p] ** k for p, k in parents.items()]), rel=1e-9
    )
    assert u[con] == pytest.approx(sum(u[p] * k for p, k in parents.items()), rel=1e-9)


def test_m2_group_shares_the_m2_nodal_factor():
    names = ["M2", "N2", "2N2", "MU2", "NU2"]
    _, pf, _ = nodal(TIME, names)
    assert np.allclose(pf, pf[0])


def test_mf_and_mm_are_modulated():
    """Mf swings by roughly 40% over the nodal cycle; Mm by about 13%."""
    fs = {c: [] for c in ("MF", "MM")}
    for day in np.linspace(0, 6798, 60):  # one 18.6 year cycle
        _, pf, _ = nodal(48622.0 + day, ["MF", "MM"])
        fs["MF"].append(pf[0])
        fs["MM"].append(pf[1])
    assert np.ptp(fs["MF"]) > 0.6
    assert np.ptp(fs["MM"]) > 0.2


def test_existing_constituents_keep_their_published_arguments():
    """The nine originally tabulated values are unchanged by this addition."""
    original = {
        "Z0": 0.0,
        "Q1": 5.877717569,
        "O1": 1.558553872,
        "P1": 6.110181633,
        "K1": 0.173003674,
        "N2": 6.050721243,
        "M2": 1.731557546,
        "S2": 0.000000000,
        "K2": 3.487600001,
    }
    for con, value in original.items():
        assert V0U[con] == pytest.approx(value)
