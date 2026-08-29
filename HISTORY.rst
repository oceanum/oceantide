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
