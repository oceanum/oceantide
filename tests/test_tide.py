"""Test tide accessor."""
from pathlib import Path
import pytest
import datetime
import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr

from oceantide import read_oceantide
from oceantide.constituents import OMEGA
from oceantide.core.utils import nodal, nodal_corrections


FILES_DIR = Path(__file__).parent / "test_files"


@pytest.fixture(scope="module")
def dset():
    filename = FILES_DIR / "oceantide.zarr"
    _dset = read_oceantide(filename)
    yield _dset


def test_accessor_created(dset):
    assert hasattr(dset, "tide")
    assert hasattr(dset.tide, "predict")


def test_predict_elevation_only(dset):
    times = [datetime.datetime(2001, 1, 1, H) for H in range(6)]
    eta = dset.tide.predict(times, components=["h"])
    assert "u" not in eta and "v" not in eta and "h" in eta


def test_predict_current_only(dset):
    times = [datetime.datetime(2001, 1, 1, H) for H in range(6)]
    eta = dset.tide.predict(times, components=["u", "v"])
    assert "h" not in eta and "u" in eta and "v" in eta


def test_always_predict_something(dset):
    times = [datetime.datetime(2001, 1, 1, H) for H in range(6)]
    with pytest.raises(ValueError):
        dset.tide.predict(times, components=[])


def test_nodal_corrections_vary_along_time():
    """f and u track the 18.6 year cycle; v0 is a constant of the constituent."""
    times = pd.date_range("2020-01-01", "2035-01-01", freq="30D")
    tsec = xr.DataArray(
        (times.values - np.datetime64("1992-01-01")) / np.timedelta64(1, "s"),
        dims="time",
        coords={"time": times},
    )
    pu, pf, v0u = nodal_corrections(tsec, ["M2", "K1", "MF"])

    assert pf.dims == ("time", "con")
    # M2 spans about +-3.7% over the cycle and Mf about +-40%.
    assert float(pf.sel(con="M2").max() - pf.sel(con="M2").min()) == pytest.approx(
        0.075, abs=0.005
    )
    assert float(pf.sel(con="MF").max() - pf.sel(con="MF").min()) > 0.6
    assert float(pu.sel(con="M2").std("time")) > 0.0
    # v0 is referenced to the 1992 epoch; the advance with time is omega*t.
    assert (v0u == v0u.isel(time=0)).all()


def test_predict_uses_the_nodal_correction_of_each_time(dset):
    """Regression: the correction used to be frozen at the first time step.

    Reproduces the sum explicitly from ``nodal`` evaluated per time, which is
    the definition the accessor is meant to implement.
    """
    times = pd.date_range("2020-06-01", periods=8, freq="90D")
    got = dset.tide.predict(times, components=["h"]).h.compute()

    tsec = (times.values - np.datetime64("1992-01-01")) / np.timedelta64(1, "s")
    cons = list(dset.con.values)
    expected = np.zeros((len(times),) + dset.h.shape[1:])
    for it, seconds in enumerate(tsec):
        pu, pf, v0u = nodal(seconds / 86400.0 + 48622.0, cons)
        for ic, con in enumerate(cons):
            theta = OMEGA[con] * seconds + v0u[ic] + pu[ic]
            z = dset.h.isel(con=ic).values
            expected[it] += pf[ic] * (np.cos(theta) * z.real - np.sin(theta) * z.imag)

    assert np.allclose(got.transpose("time", "lat", "lon").values, expected, equal_nan=True)


def test_predict_would_differ_from_a_frozen_correction(dset):
    """The fix has to actually move the answer, or it is not being applied."""
    times = pd.date_range("2020-01-01", periods=4, freq="1000D")
    got = dset.tide.predict(times, components=["h"]).h.compute()

    tsec = (times.values - np.datetime64("1992-01-01")) / np.timedelta64(1, "s")
    cons = list(dset.con.values)
    pu, pf, v0u = nodal(tsec[0] / 86400.0 + 48622.0, cons)  # the old behaviour
    frozen = np.zeros((len(times),) + dset.h.shape[1:])
    for it, seconds in enumerate(tsec):
        for ic, con in enumerate(cons):
            theta = OMEGA[con] * seconds + v0u[ic] + pu[ic]
            z = dset.h.isel(con=ic).values
            frozen[it] += pf[ic] * (np.cos(theta) * z.real - np.sin(theta) * z.imag)

    diff = np.abs(got.transpose("time", "lat", "lon").values - frozen)
    assert np.nanmax(diff) > 0.01, "prediction is still using a frozen correction"


def test_predict_stays_lazy(dset):
    """A chunked dataset must give back a dask graph, not a materialised array."""
    times = pd.date_range("2020-01-01", periods=240, freq="h")
    out = dset.chunk({"lat": 10, "lon": 10}).tide.predict(times, components=["h"])
    assert isinstance(out.h.data, da.Array)
    assert out.h.chunks is not None
