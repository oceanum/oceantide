=======
History
=======

********
Releases
********


0.8.0 (2026-08-29)
__________________

New Features
------------

* Add support for the SA, SSA, MSF, MU2, NU2, T2 and 2MS6 constituents, taking the
  supported set from 16 to 23
  (`PR2 <https://github.com/oceanum/oceantide/pull/2>`_).
* Declare compound constituents as their parents and exponents in the new SHALLOW
  table and derive their nodal factors as the product of the parents'
  (`PR2 <https://github.com/oceanum/oceantide/pull/2>`_).

Bug Fixes
---------

* Fix the angular speed of MS4 in OMEGA, previously 1.7% low and inconsistent with
  the value implied by this constituent's own entry in PERIODS
  (`PR2 <https://github.com/oceanum/oceantide/pull/2>`_).
* Define the equilibrium arguments and nodal corrections of 2N2, M4, MN4, MS4, MF,
  MM and S1, previously defined in OMEGA but missing from V0U and the nodal table
  and therefore predicted with a zero equilibrium argument and no nodal correction
  (`PR2 <https://github.com/oceanum/oceantide/pull/2>`_).
* Warn and drop constituents with no equilibrium argument in `nodal` instead of
  predicting them at an arbitrary phase
  (`PR2 <https://github.com/oceanum/oceantide/pull/2>`_).
* Support both zarr v2 and v3 when writing the oceantide zarr format. The
  FixedScaleOffset filters are now defined from numcodecs or from zarr codecs
  according to the zarr format being written, and codec encodings inherited from
  the file being read are dropped from the coordinates so they cannot conflict with
  the format being written.

Internal Changes
----------------

* Writing the oceantide format no longer modifies the encoding of the dataset
  passed to the accessor.
* Test the values of datasets read back after writing the oceantide formats, and
  test writing the zarr format in both zarr formats 2 and 3.


0.9.0 (unreleased)
__________________

Bug Fixes
---------

* Nodal corrections are evaluated at every prediction time instead of being
  held at the value they take at the first time step. On a 3 m tide the frozen
  correction cost 1.4 cm RMS over a year, 5.9 cm over five years and 9.5 cm
  over a full 18.6 year nodal cycle. Predictions over more than a few days will
  change, by design.
* Fix the K1 nodal amplitude factor, which used 0.01554 where its own phase
  correction and every published table use 0.1554. Worth up to 1.2% of K1
  amplitude.
* Fix the MSF nodal amplitude factor, which was inverted: compound constituents
  raise the parent factors to the magnitude of their exponents, and only the
  phase correction carries the sign. Worth up to 7.8% of MSF amplitude.
* Derive compound constituent frequencies and equilibrium arguments from their
  parents. M4 and MN4 were tabulated a part in 1e10 away from twice M2 and from
  M2 + N2, which is 6 and 12 degrees of phase 34 years from the 1992 epoch, and
  2MS6 carried M4's equilibrium argument verbatim.
* ``to_oceantide`` raises instead of silently wrapping when a value falls
  outside the range the int16 packing can represent. A 25 m amplitude used to
  round-trip as -15 m and a 15000 m depth as 2999 m. Pass ``check_range=False``
  to skip the check and the pass over the data it costs.
* The OTIS readers no longer produce infinite currents at nodes bordering land,
  where transport was divided by a zero node depth. h, u and v now share one
  mask.
* ``to_oceantide`` and ``read_oceantide`` handle datasets carrying only some of
  h, u and v, which the accessor has supported since 0.7.0. Writing one without
  ``dep`` raised ``KeyError``.
* Predictions no longer inherit the codec encoding of the file the constituents
  were read from, which made writing one back to zarr fail across format
  versions.
* The release workflow installed the published package over the working tree,
  so the tests gating a release ran against the previous release.

Internal Changes
----------------

* ``nodal()`` accepts a time of any shape; a scalar time behaves as before.
* New ``nodal_corrections()`` returns the corrections as DataArrays, lazily
  when the time axis is dask backed.
* ``PERIODS`` is derived from ``OMEGA`` rather than tabulated alongside it.
* Prediction is checked end to end against a committed pyTMD reference, and the
  constituent tables against the Doodson expansion.
* ``oceantide`` on the command line is a working ``convert`` and ``info``
  command rather than the cookiecutter placeholder.
* Declare the ``pyyaml`` dependency, which was only ever present transitively.
* ``requires-python`` is now >=3.11, matching what the dependencies resolve for.
* Tests run on every push and pull request across Python 3.11 to 3.13, plus a
  job pinning ``zarr<3`` so both branches of the zarr writer are exercised.
* Remove the vendored ``oceantide.ellipse`` module, which nothing imported and
  which had not been runnable since matplotlib 3.0.
* Write a real README and usage documentation; the docs build again.


0.7.0 (2024-01-31)
__________________

Internal Changes
----------------

* Support datasets in tide accessors with only a subset of h, u, v variables
  (`PR1 <https://github.com/oceanum/oceantide/pull/1>`_).


0.6.0 (2023-12-18)
__________________

Internal Changes
----------------

* Make gcsfs an extra dependency.
* Add extra dependencies for docs.
* Change sytle of docstrings to numpydoc.

0.5.0 (2023-07-26)
__________________

Internal Changes
----------------

* Redefine packaging via pyproject.toml.

Bug Fixes
---------

* Fix bug when in oceantide.core.otis when indexing with a dask boolean array.


0.1.0 (2020-10-02)
__________________

* First release on PyPI.
