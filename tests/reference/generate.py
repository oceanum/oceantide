"""Generate the pyTMD reference fixture for oceantide's prediction path.

Needs pyTMD installed, which oceantide does not depend on:

    python tests/reference/generate.py

Rerun only when the reference genuinely needs to change, and say why in
the commit -- the point of the fixture is that it does not move.
"""
import json, datetime, numpy as np, pandas as pd, xarray as xr
import pyTMD
from pyTMD.predict import time_series as tmd

# Constituents where oceantide and pyTMD carry identical frequencies and
# equilibrium arguments, so any disagreement is a real defect rather than a
# difference of table.
CONS = ["M2", "S2", "N2", "K2", "K1", "O1", "P1", "Q1", "MS4"]
AMP = [1.5, 0.45, 0.30, 0.12, 0.35, 0.22, 0.11, 0.05, 0.03]
PHA = [10.0, 55.0, 100.0, 145.0, 190.0, 235.0, 280.0, 325.0, 12.0]

# 48 hourly steps at the epoch, then quarterly out to 2040 so the nodal cycle
# is sampled more than twice over.
times = pd.DatetimeIndex(
    list(pd.date_range("1992-01-01", periods=48, freq="h"))
    + list(pd.date_range("1992-01-01", "2040-01-01", freq="91D"))
)

z = (np.array(AMP) * np.exp(-1j * np.deg2rad(PHA))).reshape(-1, 1, 1)
ds = xr.Dataset({c.lower(): (("y", "x"), z[i]) for i, c in enumerate(CONS)},
                coords={"y": [-35.0], "x": [175.0]})
for c in CONS:
    ds[c.lower()].attrs["units"] = "m"

tdays = (times.values - np.datetime64("1992-01-01")) / np.timedelta64(1, "D")
eta = np.asarray(tmd(tdays, ds, corrections="OTIS").values).ravel()

out = {
    "description": (
        "Reference tide elevations for oceantide's prediction path. Generated "
        "with pyTMD, an independent implementation of the same OTIS constituent "
        "tables and nodal corrections, so that the whole chain -- epoch, "
        "angular speeds, equilibrium arguments, time-varying nodal f and u, and "
        "the complex amplitude convention -- is checked against something other "
        "than itself."
    ),
    "source": f"pyTMD {pyTMD.version.full_version}, "
              "pyTMD.predict.time_series(corrections='OTIS')",
    "generated": datetime.date.today().isoformat(),
    "epoch": "1992-01-01T00:00:00Z",
    "constituents": CONS,
    "amplitude_m": AMP,
    "phase_deg_gmt": PHA,
    "times": [str(t) for t in times],
    "elevation_m": [float(v) for v in eta],
}
with open("tests/reference/pytmd_elevation.json", "w") as stream:
    json.dump(out, stream, indent=1)
print("wrote", len(eta), "values, range", eta.min().round(4), eta.max().round(4))
