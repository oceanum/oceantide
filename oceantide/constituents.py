"""Constituents.

Each constituent needs three entries to be predicted correctly:

``PERIODS``
    Period in hours, used for reporting. Derived from ``OMEGA``.
``OMEGA``
    Angular speed in rad/s, the rate at which the argument advances.
``V0U``
    The equilibrium argument V0 in radians at 1992-01-01T00:00 UTC, which is the
    zero point of the ``tsec`` axis :meth:`oceantide.tide.Tide.predict` builds.

A constituent listed in ``OMEGA`` but missing from ``V0U`` is not simply
unsupported -- it is predicted with a zero equilibrium argument, i.e. an
arbitrary phase reference, which is worse than dropping it. ``nodal()`` warns
when that happens; keep the tables in step.

Only the astronomically primary constituents are tabulated. Compound (shallow
water) constituents and the periods are derived below, because a table that
repeats a quantity it could compute is a table that will eventually disagree
with itself -- which is how MS4 came to carry a frequency 1.7% off its own
period, and M4 one that was not quite twice M2's.
"""

import math


# Compound (shallow water) constituents, as the parent constituents they are
# built from and the exponent of each:
#
#   M4   = 2*M2       MN4  = M2 + N2     MS4  = M2 + S2
#   2MS6 = 2*M2 + S2  MSF  = S2 - M2
#
# The angular speed and the equilibrium argument are both this same combination
# of the parents', and are derived as such below. The nodal amplitude factor is
# the product of the parents' factors raised to the *magnitude* of these
# exponents, and the nodal phase correction the signed weighted sum, which is
# how Schureman and Foreman define them (see ``oceantide.core.utils.nodal``).
SHALLOW = {
    "M4": {"M2": 2},
    "MN4": {"M2": 1, "N2": 1},
    "MS4": {"M2": 1, "S2": 1},
    "2MS6": {"M2": 2, "S2": 1},
    "MSF": {"S2": 1, "M2": -1},
}

# Angular speed in rad/s. These agree to every digit given with the OTIS
# `constit.h` table that TPXO and the other OTIS-format atlases are fit
# against, and with pyTMD's copy of it.
OMEGA = {
    "Z0": 0.0,
    "SA": 1.990969e-07,
    "SSA": 3.982128e-07,
    "MM": 0.026392e-04,
    "MF": 0.053234e-04,
    "Q1": 6.495854e-05,
    "O1": 6.759774e-05,
    "P1": 7.252295e-05,
    "S1": 7.2722083e-05,
    "K1": 7.292117e-05,
    "2N2": 1.352405e-04,
    "MU2": 1.355937e-04,
    "N2": 1.378797e-04,
    "NU2": 1.382329e-04,
    "M2": 1.405189e-04,
    "T2": 1.452450e-04,
    "S2": 1.454441e-04,
    "K2": 0.0001458423,
}

# V0 at 1992-01-01T00:00 UTC, radians. The nine originally present are kept at
# their published values; the rest are computed from the Doodson numbers at the
# same epoch and agree with the originals to better than 0.03 degrees.
V0U = {
    "Z0": 0.0,
    "SA": 6.232722191,
    "SSA": 3.487570030,
    "MM": 1.963863349,
    "MF": 1.755713874,
    "Q1": 5.877717569,
    "O1": 1.558553872,
    "P1": 6.110181633,
    "S1": 0.223451805,
    "K1": 0.173003674,
    "2N2": 4.087314766,
    "MU2": 3.463712312,
    "N2": 6.050721243,
    "NU2": 5.427575661,
    "M2": 1.731557546,
    "T2": 0.050463116,
    "S2": 0.000000000,
    "K2": 3.487600001,
}


def _combine(table: dict, parents: dict) -> float:
    """Sum a parent quantity weighted by the compound exponents."""
    return sum(table[parent] * power for parent, power in parents.items())


for _con, _parents in SHALLOW.items():
    OMEGA[_con] = _combine(OMEGA, _parents)
    V0U[_con] = _combine(V0U, _parents) % (2 * math.pi)
del _con, _parents

# Frequency order, so the tables read as a spectrum.
OMEGA = dict(sorted(OMEGA.items(), key=lambda item: item[1]))
V0U = {con: V0U[con] for con in OMEGA}

# Period in hours. Derived, so it cannot drift away from the angular speed.
PERIODS = {
    con: (2 * math.pi / omega / 3600.0) if omega else 0.0
    for con, omega in OMEGA.items()
}
