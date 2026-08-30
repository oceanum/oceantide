"""End to end check of the prediction path against an independent library.

The other tests check the constituent tables and the nodal corrections against
published coefficients. This one checks the whole chain -- epoch, angular
speeds, equilibrium arguments, time-varying f and u, the complex amplitude
convention and the sum over constituents -- against elevations produced by
pyTMD, a separate implementation of the same OTIS lineage.

The fixture is committed so the test needs no optional dependency. It was
written by ``tests/reference/generate.py``; regenerate it only when there is a
reason to, and say so in the commit.

Only constituents whose frequency and equilibrium argument are identical in
both libraries are included, so the tolerance can be tight enough to catch a
real defect. Where the two legitimately differ -- pyTMD's Mm drops a term this
carries, and its M4 and MN4 frequencies are truncated where these are derived
from M2 and N2 -- the constituent-level tests in ``test_constituents.py`` are
the reference instead.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import oceantide  # noqa: F401  (registers the accessor)


REFERENCE = Path(__file__).parent / "reference" / "pytmd_elevation.json"

# pyTMD accumulates in float64 from the same coefficients, so agreement should
# be at rounding level. A micrometre is far below anything tidally meaningful
# and far above float64 noise on a metre-scale sum.
ATOL = 1e-6


@pytest.fixture(scope="module")
def reference():
    with open(REFERENCE) as stream:
        return json.load(stream)


@pytest.fixture(scope="module")
def constituents(reference):
    """The reference constants as an oceantide dataset at a single point."""
    amp = np.array(reference["amplitude_m"])
    pha = np.deg2rad(reference["phase_deg_gmt"])
    z = (amp * np.exp(-1j * pha)).reshape(-1, 1, 1)
    return xr.Dataset(
        {
            "h": (("con", "lat", "lon"), z),
            "dep": (("lat", "lon"), np.array([[50.0]])),
        },
        coords={
            "con": np.array(reference["constituents"], dtype="U4"),
            "lat": [-35.0],
            "lon": [175.0],
        },
    )


def test_predicted_elevation_matches_pytmd(reference, constituents):
    times = pd.DatetimeIndex(reference["times"])
    expected = np.array(reference["elevation_m"])

    got = constituents.tide.predict(times, components=["h"]).h.compute()
    got = got.transpose("time", "lat", "lon").values.ravel()

    assert np.abs(got - expected).max() < ATOL


def test_predicted_elevation_matches_pytmd_over_a_chunked_time_axis(
    reference, constituents
):
    """The lazy path must give the same answer as the eager one."""
    times = pd.DatetimeIndex(reference["times"])
    expected = np.array(reference["elevation_m"])

    got = constituents.chunk().tide.predict(
        times, components=["h"], time_chunk=7
    )
    assert got.h.chunks is not None
    got = got.h.compute().transpose("time", "lat", "lon").values.ravel()

    assert np.abs(got - expected).max() < ATOL


def test_reference_spans_more_than_one_nodal_cycle(reference):
    """A reference confined to one epoch would not test the corrections."""
    times = pd.DatetimeIndex(reference["times"])
    assert (times[-1] - times[0]).days > 2 * 6798
