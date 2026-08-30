=========
oceantide
=========

Ocean tide prediction on top of xarray.

oceantide adds a ``tide`` accessor to ``xarray.Dataset`` objects holding tidal
harmonic constituents, and predicts elevations and currents from them. The
prediction is vectorised over constituents, time and space at once, and stays
lazy when the constituents are dask backed, so a long series over a large grid
can be streamed to disk without being held in memory.

.. code-block:: python

    import pandas as pd
    import oceantide

    dset = oceantide.read_otis_netcdf("DATA/Model_tpxo9")

    times = pd.date_range("2026-01-01", "2026-02-01", freq="30min")
    eta = dset.sel(lon=slice(172, 179), lat=slice(-42, -34)).tide.predict(times)

    eta.h.isel(lat=0, lon=0).plot()

``predict`` returns a dataset of ``h`` (elevation, m), ``u`` and ``v``
(depth averaged velocity, m/s) over a new ``time`` dimension, keeping whatever
spatial dimensions the constituents had. Pass ``components`` to compute only
some of them.

For a time series at a location, interpolate first and predict after. The
constituents are stored as complex amplitudes, so interpolating them is a
linear operation and does the right thing, which is not true of amplitude and
phase held separately:

.. code-block:: python

    site = dset.interp(lon=174.8, lat=-36.8)
    eta = site.tide.predict(times, components=["h"])

Amplitude and Greenwich phase are available directly:

.. code-block:: python

    dset.tide.amplitude("h")   # metres
    dset.tide.phase("h")       # degrees, relative to GMT


Installation
------------

.. code-block:: console

    pip install oceantide


Formats
-------

Reading:

============================ ====================================================
``read_otis_binary``         OTIS binary, as distributed with TPXO
``read_otis_netcdf``         OTIS netcdf
``read_otis_atlas_netcdf``   OTIS atlas netcdf, one file per constituent
``read_oceantide``           the oceantide format, netcdf or zarr
============================ ====================================================

Writing, as accessor methods: ``to_oceantide``, ``to_otis_binary``,
``to_otis_netcdf``.

The oceantide format keeps the complex constituents split into real and
imaginary variables packed as int16 against a fixed scale, which makes the
files roughly a quarter the size of the equivalent complex64 and lets netcdf
hold them. The packing covers amplitudes within ±20 and depths from 0 to
12000 m; writing anything outside that raises rather than wrapping silently.

There is also a small command line interface:

.. code-block:: console

    oceantide convert -r otis_netcdf DATA/Model_tpxo9 tpxo9.zarr
    oceantide info tpxo9.zarr


Constituents and conventions
----------------------------

23 constituents are supported: Z0, Sa, Ssa, Mm, MSf, Mf, Q1, O1, P1, S1, K1,
2N2, mu2, N2, nu2, M2, T2, S2, K2, MN4, M4, MS4 and 2MS6. Constituents present
in a dataset but not in that list are dropped with a warning rather than
predicted at an arbitrary phase.

Angular speeds and equilibrium arguments follow the OTIS ``constit.h`` tables
that the TPXO atlases are fit against. Nodal corrections use the OTIS
coefficients and are evaluated at every prediction time rather than held fixed
at the start of the series, which matters for anything longer than a few days:
on a 3 m tide a fixed correction costs about 1 cm RMS over a year and 10 cm
over a full 18.6 year nodal cycle. Phases are relative to GMT. The prediction
path is tested against `pyTMD <https://github.com/pyTMD/pyTMD>`_, an
independent implementation of the same tables, and agrees to below a
micrometre.


Scope
-----

oceantide is deliberately narrow. It reads OTIS-family constituents, predicts
from them, and writes a compact archive format, and it is built to do that
lazily and quickly over large grids.

It does not do harmonic analysis, and it does not read the FES, GOT, EOT,
HAMTIDE or DTU families. If you need broader model coverage, minor constituent
inference, load tides or solid earth tides, use
`pyTMD <https://github.com/pyTMD/pyTMD>`_, which is more comprehensive and
actively maintained; `eo-tides <https://github.com/GeoscienceAustralia/eo-tides>`_
builds satellite Earth observation workflows on top of it. For harmonic
analysis of a gauge record, use `UTide <https://github.com/wesleybowman/UTide>`_.


License
-------

MIT.
