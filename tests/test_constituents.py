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


# Doodson numbers (tau, s, h, p, ps) of each constituent, and the mean motions
# of the astronomical arguments in degrees per day from the Explanatory
# Supplement. Together these fix every angular speed from first principles,
# independently of the OTIS table OMEGA is transcribed from. N' is omitted: it
# does not enter the frequencies, only the nodal corrections f and u.
DOODSON = {
    "SA": (0, 0, 1, 0, -1),
    "SSA": (0, 0, 2, 0, 0),
    "MM": (0, 1, 0, -1, 0),
    "MSF": (0, 2, -2, 0, 0),
    "MF": (0, 2, 0, 0, 0),
    "Q1": (1, -2, 0, 1, 0),
    "O1": (1, -1, 0, 0, 0),
    "P1": (1, 1, -2, 0, 0),
    "S1": (1, 1, -1, 0, 1),
    "K1": (1, 1, 0, 0, 0),
    "2N2": (2, -2, 0, 2, 0),
    "MU2": (2, -2, 2, 0, 0),
    "N2": (2, -1, 0, 1, 0),
    "NU2": (2, -1, 2, -1, 0),
    "M2": (2, 0, 0, 0, 0),
    "T2": (2, 2, -3, 0, 1),
    "S2": (2, 2, -2, 0, 0),
    "K2": (2, 2, 0, 0, 0),
    "MN4": (4, -1, 0, 1, 0),
    "M4": (4, 0, 0, 0, 0),
    "MS4": (4, 2, -2, 0, 0),
    "2MS6": (6, 2, -2, 0, 0),
}

# Mean motions in degrees per day: moon, sun, lunar perigee, solar perigee.
RATE_S, RATE_H, RATE_P, RATE_PS = 13.17639648, 0.98564736, 0.11140353, 0.0000470684
RATE_TAU = 360.0 + RATE_H - RATE_S  # one mean lunar day


@pytest.mark.parametrize("con", sorted(DOODSON))
def test_omega_matches_the_doodson_expansion(con):
    """Angular speeds must be the combination their Doodson numbers declare.

    PERIODS is derived from OMEGA and the compound speeds are derived from
    their parents, so neither can be checked against the other. This checks the
    whole table against the astronomy instead.
    """
    tau, s, h, p, ps = DOODSON[con]
    degrees_per_day = (
        tau * RATE_TAU + s * RATE_S + h * RATE_H + p * RATE_P + ps * RATE_PS
    )
    expected = np.deg2rad(degrees_per_day) / 86400.0
    assert OMEGA[con] == pytest.approx(expected, rel=5e-6)
    assert PERIODS[con] == pytest.approx(360.0 / degrees_per_day * 24.0, rel=5e-6)


def test_ms4_is_the_sum_of_m2_and_s2():
    """Regression: MS4 was 2.811149e-04, which is not M2 + S2."""
    assert OMEGA["MS4"] == OMEGA["M2"] + OMEGA["S2"]


@pytest.mark.parametrize("con,parents", sorted(SHALLOW.items()))
def test_compound_speeds_are_the_sum_of_their_parents(con, parents):
    """Exactly, not approximately: these are derived, not transcribed.

    M4 was tabulated as 2.810379e-04 against 2*OMEGA["M2"] of 2.810378e-04.
    That 1e-10 rad/s is 6 degrees of M4 phase 34 years from the 1992 epoch.
    """
    assert OMEGA[con] == sum(OMEGA[p] * k for p, k in parents.items())


@pytest.mark.parametrize("con,parents", sorted(SHALLOW.items()))
def test_compound_equilibrium_arguments_follow_their_parents(con, parents):
    """Exactly, as for the speeds: derived from the parents, not transcribed.

    2MS6 had been given M4's argument verbatim, which happens to be within
    0.03 degrees of the right answer and so read as plausible.
    """
    expected = sum(V0U[p] * k for p, k in parents.items()) % (2 * np.pi)
    assert V0U[con] == expected


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


@pytest.mark.parametrize("con,parents", sorted(SHALLOW.items()))
def test_compound_nodal_factors_are_products_of_parents(con, parents):
    """f multiplies the magnitudes; only u carries the sign of the exponent."""
    names = sorted({con} | set(parents))
    pu, pf, _ = nodal(TIME, names)
    f = dict(zip(names, pf))
    u = dict(zip(names, pu))
    assert f[con] == pytest.approx(
        np.prod([f[p] ** abs(k) for p, k in parents.items()]), rel=1e-9
    )
    assert u[con] == pytest.approx(sum(u[p] * k for p, k in parents.items()), rel=1e-9)


def test_msf_nodal_factor_is_not_inverted():
    """Regression: MSF = S2 - M2 gave f = 1/f(M2) instead of f(M2).

    f(S2) is unity, so MSF's factor is exactly M2's and its angle exactly the
    negative of M2's. The signed form inverted the factor, worth up to 7.8%
    across the nodal cycle.
    """
    worst = 0.0
    for mjd in CYCLE:
        pu, pf, _ = nodal(mjd, ["M2", "MSF"])
        assert pf[1] == pytest.approx(pf[0], rel=1e-12)
        assert pu[1] == pytest.approx(-pu[0], rel=1e-12, abs=1e-15)
        worst = max(worst, abs(1.0 / pf[0] - pf[0]) / pf[0])
    # The spread the inverted form used to introduce, for the record.
    assert 100 * worst == pytest.approx(7.8, abs=0.2)


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


# Complex nodal modulation terms from the OTIS ``nodal.f`` coefficients, as
# (real, imag) polynomials in (1, cosN, cos2N) and (sinN, sin2N). f is the
# modulus of this term and u its argument; spelling them out here checks the
# library's numbers against the published ones rather than against itself.
NODAL_TERMS = {
    "M2": ((1.0, -0.03731, 0.00052), (-0.03731, 0.00052)),
    "N2": ((1.0, -0.03731, 0.00052), (-0.03731, 0.00052)),
    "2N2": ((1.0, -0.03731, 0.00052), (-0.03731, 0.00052)),
    "MU2": ((1.0, -0.03731, 0.00052), (-0.03731, 0.00052)),
    "NU2": ((1.0, -0.03731, 0.00052), (-0.03731, 0.00052)),
    "K1": ((1.0, 0.1158, -0.0029), (-0.1554, 0.0029)),
    "K2": ((1.0, 0.2852, 0.0324), (-0.3108, -0.0324)),
}

# A full nodal cycle, sampled off the extremes where sinN vanishes and the
# imaginary term (the one K1 got wrong) drops out of f.
CYCLE = 48622.0 + np.linspace(0.0, 6798.0, 97)


def _expected_polar(con, mjd):
    """f and u in degrees from the published coefficients, independently."""
    from oceantide.core.utils import astrol

    (a0, a1, a2), (b1, b2) = NODAL_TERMS[con]
    _, _, _, N = astrol(mjd)
    n = np.deg2rad(N)
    real = a0 + a1 * np.cos(n) + a2 * np.cos(2 * n)
    imag = b1 * np.sin(n) + b2 * np.sin(2 * n)
    return np.abs(real + 1j * imag), np.rad2deg(np.angle(real + 1j * imag))


@pytest.mark.parametrize("con", sorted(NODAL_TERMS))
def test_nodal_factor_and_angle_match_the_published_term(con):
    """f and u must be the modulus and argument of the same complex term."""
    for mjd in CYCLE:
        pu, pf, _ = nodal(mjd, [con])
        f_exp, u_exp = _expected_polar(con, mjd)
        assert pf[0] == pytest.approx(f_exp, rel=1e-12), con
        assert np.rad2deg(pu[0]) == pytest.approx(u_exp, rel=1e-12, abs=1e-12), con


def test_k1_nodal_factor_is_not_the_dropped_digit_form():
    """Regression: f used 0.01554 where u (and every reference) uses 0.1554.

    The two forms agree wherever sinN vanishes, so this compares them across
    the cycle rather than at a single time.
    """
    from oceantide.core.utils import astrol

    _, _, _, N = astrol(CYCLE)
    n = np.deg2rad(N)
    real = 1.0 + 0.1158 * np.cos(n) - 0.0029 * np.cos(2 * n)
    typo = np.hypot(real, 0.01554 * np.sin(n) - 0.0029 * np.sin(2 * n))
    got = np.array([nodal(mjd, ["K1"])[1][0] for mjd in CYCLE])

    assert np.abs(got - typo).max() > 0.01, "K1 still carries the dropped digit"
    # The error the dropped digit used to introduce, for the record.
    assert (100 * np.abs(got - typo) / got).max() == pytest.approx(1.19, abs=0.05)
