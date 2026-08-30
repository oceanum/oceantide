=====
Usage
=====

Reading constituents
--------------------

Each reader returns an ``xarray.Dataset`` carrying the complex constituents
``h``, ``u`` and ``v`` over a ``con`` dimension, the depth ``dep``, and the
``tide`` accessor:

.. code-block:: python

    import oceantide

    dset = oceantide.read_otis_netcdf("DATA/Model_tpxo9")
    dset = oceantide.read_otis_binary("DATA/Model_tpxo9")
    dset = oceantide.read_otis_atlas_netcdf("DATA/atlas", x0=170, x1=180, y0=-45, y1=-33)
    dset = oceantide.read_oceantide("tpxo9.zarr")

The atlas reader covers the whole globe at high resolution, so bound it with
``x0``, ``x1``, ``y0`` and ``y1`` unless you have the memory for all of it.

Slice before predicting. Everything downstream is lazy, but the smaller the
grid the smaller the graph:

.. code-block:: python

    dset = dset.sel(lon=slice(172, 179), lat=slice(-42, -34))


Predicting
----------

.. code-block:: python

    import pandas as pd

    times = pd.date_range("2026-01-01", "2026-02-01", freq="30min")
    pred = dset.tide.predict(times)

``pred`` has ``h`` in metres and ``u``, ``v`` as depth averaged velocity in m/s,
over a new ``time`` dimension. Restrict the work with ``components``:

.. code-block:: python

    eta = dset.tide.predict(times, components=["h"])

``times`` may be a list of datetimes, a ``pandas.DatetimeIndex``, a numpy array
of ``datetime64``, or a ``DataArray`` of times. A ``DataArray`` lets the time
axis carry other dimensions, for instance a different set of times per site.

Nodal corrections are evaluated at every time in the series, so predictions
over years are as valid as predictions over days.


Working at sites
----------------

Interpolate the constituents, then predict. They are stored as complex
amplitudes, so interpolation is linear and correct; interpolating amplitude and
phase separately is not:

.. code-block:: python

    import xarray as xr

    sites = dset.interp(
        lon=xr.DataArray([174.8, 177.9], dims="site"),
        lat=xr.DataArray([-36.8, -39.5], dims="site"),
    )
    pred = sites.tide.predict(times)

Use ``sel(..., method="nearest")`` instead to take the containing cell without
interpolating, which avoids pulling in land values near the coast.


Staying lazy
------------

``predict`` builds a dask graph and does not compute. Chunk the constituents
and write the result straight out:

.. code-block:: python

    pred = dset.chunk({"lat": 100, "lon": 100}).tide.predict(times)
    pred.to_zarr("tide_2026.zarr")

``time_chunk`` sets the chunk size along the time axis and defaults to 50.
Raise it when the spatial grid is small and the series long.


Amplitude and phase
-------------------

.. code-block:: python

    dset.tide.amplitude("h")   # metres
    dset.tide.phase("h")       # degrees relative to GMT

Both take ``"h"``, ``"u"`` or ``"v"``.


Writing
-------

.. code-block:: python

    dset.tide.to_oceantide("tpxo9.zarr")
    dset.tide.to_oceantide("tpxo9.nc")
    dset.tide.to_otis_binary("outdir")
    dset.tide.to_otis_netcdf("outdir")

The oceantide format packs into int16 against a fixed scale: amplitudes within
±20 and depths from 0 to 12000 m. Values outside that cannot be represented and
raise a ``ValueError`` naming the variable. Pass ``check_range=False`` to skip
the check, and the pass over the data it costs, where the range is known to be
safe.


Command line
------------

.. code-block:: console

    oceantide convert -r otis_netcdf DATA/Model_tpxo9 tpxo9.zarr
    oceantide info tpxo9.zarr

``convert`` reads any supported input format and writes the oceantide format,
choosing netcdf or zarr from the output extension. ``info`` prints the
constituents, variables and extent.
