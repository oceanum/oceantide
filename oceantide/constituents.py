"""Constituents.

Each constituent needs three entries to be predicted correctly:

``PERIODS``
    Period in hours, used for reporting.
``OMEGA``
    Angular speed in rad/s, the rate at which the argument advances.
``V0U``
    The equilibrium argument V0 in radians at 1992-01-01T00:00 UTC, which is the
    zero point of the ``tsec`` axis :meth:`oceantide.tide.Tide.predict` builds.

A constituent listed in ``OMEGA`` but missing from ``V0U`` is not simply
unsupported -- it is predicted with a zero equilibrium argument, i.e. an
arbitrary phase reference, which is worse than dropping it. ``nodal()`` warns
when that happens; keep the three tables in step.
"""

PERIODS = {
    "Z0": 0.0,
    "SA": 8766.23178267,
    "SSA": 4382.906453154,
    "MM": 27.55463 * 24,
    "MSF": 354.367061882,
    "MF": 13.66079 * 24,
    "Q1": 26.8683567047119,
    "O1": 25.8193397521973,
    "P1": 24.0658893585205,
    "S1": 24.0,
    "K1": 23.9344692230225,
    "2N2": 12.90537297,
    "MU2": 12.871757597,
    "N2": 12.6583499908447,
    "NU2": 12.626004399,
    "M2": 12.420599937439,
    "T2": 12.016449193,
    "S2": 12.0,
    "K2": 11.9672346115112,
    "MN4": 6.269173724,
    "M4": 6.210300601,
    "MS4": 6.103339275,
    "2MS6": 4.092387536,
}

OMEGA = {
    "Z0": 0,
    "SA": 1.990969e-07,
    "SSA": 3.982128e-07,
    "MM": 0.026392e-04,
    "MSF": 4.925202e-06,
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
    "MN4": 2.783987e-04,
    "M4": 2.810379e-04,
    # Was 2.811149e-04, which is neither M2+S2 nor any other combination and
    # disagreed with this constituent's own PERIODS entry by 1.7%.
    "MS4": 2.859630e-04,
    "2MS6": 4.264819e-04,
}

# V0 at 1992-01-01T00:00 UTC, radians. The nine originally present are kept at
# their published values; the rest are computed from the Doodson numbers at the
# same epoch and agree with the originals to better than 0.03 degrees.
V0U = {
    "Z0": 0.0,
    "SA": 6.232722191,
    "SSA": 3.487570030,
    "MM": 1.963863349,
    "MSF": 4.551329151,
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
    "MN4": 1.499848963,
    "M4": 3.463712312,
    "MS4": 1.731856156,
    "2MS6": 3.463712312,
}

# Compound (shallow water) constituents, as the parent constituents they are
# built from and the exponent of each. The nodal amplitude factor is the product
# of the parents' factors raised to the *magnitude* of these exponents, and the
# phase correction is the signed weighted sum, which is how Schureman and
# Foreman define them (see the note in ``oceantide.core.utils.nodal``).
#
#   M4   = 2*M2       MN4  = M2 + N2     MS4  = M2 + S2
#   2MS6 = 2*M2 + S2  MSF  = S2 - M2
#
# Each identity is confirmed by the constituent's angular speed and equilibrium
# argument both being the same combination of its parents'.
SHALLOW = {
    "M4": {"M2": 2},
    "MN4": {"M2": 1, "N2": 1},
    "MS4": {"M2": 1, "S2": 1},
    "2MS6": {"M2": 2, "S2": 1},
    "MSF": {"S2": 1, "M2": -1},
}
