"""Test tide accessor."""
from pathlib import Path
import pytest
import datetime
import numpy as np
import xarray as xr

from oceantide import read_otis_netcdf, read_otis_binary, read_oceantide
from oceantide.output.oceantide import ZARR_VERSION


FILES_DIR = Path(__file__).parent / "test_files"

# Tolerances for the int16 packing applied when writing the oceantide format
DEP_ATOL = 0.2
AMP_ATOL = 1e-3


def _assert_roundtrip(dset, dset2):
    assert set(dset2.data_vars) == set(dset.data_vars)
    for varname, dvar in dset.data_vars.items():
        atol = DEP_ATOL if varname == "dep" else AMP_ATOL
        other = dset2[varname].values
        assert np.allclose(dvar.values, other, atol=atol, equal_nan=True)


def test_read_otis_binary_with_filename():
    read_otis_binary(FILES_DIR / "otis_binary/Model_rag")


def test_read_otis_binary_without_filename():
    read_otis_binary(
        gfile=FILES_DIR / "otis_binary/grid_rag",
        hfile=FILES_DIR / "otis_binary/h_rag",
        ufile=FILES_DIR / "otis_binary/u_rag",
    )


def test_otis_binary_correct_args():
    with pytest.raises(ValueError):
        read_otis_netcdf(
            hfile=FILES_DIR / "otis_binary/h_rag",
            ufile=FILES_DIR / "otis_binary/u_rag",
        )


def test_read_otis_netcdf_with_filename():
    read_otis_netcdf(FILES_DIR / "otis_netcdf/Model_test")


def test_read_otis_netcdf_without_filename():
    read_otis_netcdf(
        gfile=FILES_DIR / "otis_netcdf/grid.test.nc",
        hfile=FILES_DIR / "otis_netcdf/hf.test.nc",
        ufile=FILES_DIR / "otis_netcdf/uv.test.nc",
    )


def test_otis_netcdf_correct_args():
    with pytest.raises(ValueError):
        read_otis_netcdf(
            hfile=FILES_DIR / "otis_netcdf/hf.test.nc",
            ufile=FILES_DIR / "otis_netcdf/uv.test.nc",
        )


def test_read_write_oceantide_netcdf(tmpdir):
    dset = read_oceantide(FILES_DIR / "oceantide.nc")
    filename = str(tmpdir / "newoceantide.nc")
    dset.tide.to_oceantide(filename)
    _assert_roundtrip(dset, read_oceantide(filename))


def test_read_write_oceantide_zarr(tmpdir):
    dset = read_oceantide(FILES_DIR / "oceantide.zarr")
    filename = str(tmpdir / "newoceantide.zarr")
    dset.tide.to_oceantide(filename)
    _assert_roundtrip(dset, read_oceantide(filename))


@pytest.mark.parametrize("zarr_format", [2, 3])
def test_read_write_oceantide_zarr_format(tmpdir, zarr_format):
    """Both zarr formats are supported when writing with zarr>=3."""
    if ZARR_VERSION < 3:
        pytest.skip("The zarr_format option requires zarr>=3")
    dset = read_oceantide(FILES_DIR / "oceantide.zarr")
    filename = str(tmpdir / f"newoceantide{zarr_format}.zarr")
    dset.tide.to_oceantide(filename, zarr_format=zarr_format)
    _assert_roundtrip(dset, read_oceantide(filename))


def test_supported_oceantide_formats(tmpdir):
    dset = read_oceantide(FILES_DIR / "oceantide.zarr")
    with pytest.raises(ValueError):
        dset.tide.to_oceantide(tmpdir / "newoceantide.zarr", file_format="txt")
    with pytest.raises(ValueError):
        dset.tide.to_oceantide(tmpdir / "newoceantide.txt")
    dset.tide.to_oceantide(tmpdir / "newoceantide.txt", file_format="netcdf")


def test_write_otis_netcdf(tmpdir):
    dset = read_otis_netcdf(FILES_DIR / "otis_netcdf/Model_test")
    dset.tide.to_otis_netcdf(dirname=tmpdir, suffix="")
    dset2 = read_otis_netcdf(tmpdir / "model")


def test_write_otis_binary(tmpdir):
    dset = read_otis_binary(FILES_DIR / "otis_binary/Model_rag")
    dset.tide.to_otis_binary(dirname=tmpdir, suffix="")
    dset2 = read_otis_binary(tmpdir / "model")


def _one_point(h=1.0 + 0.5j, dep=50.0):
    """Smallest dataset the oceantide writer accepts."""
    return xr.Dataset(
        {
            "h": (("con", "lat", "lon"), np.array([[[h]]])),
            "u": (("con", "lat", "lon"), np.array([[[0.1 + 0.0j]]])),
            "v": (("con", "lat", "lon"), np.array([[[0.1 + 0.0j]]])),
            "dep": (("lat", "lon"), np.array([[dep]])),
        },
        coords={"con": np.array(["M2"], dtype="U4"), "lat": [0.0], "lon": [0.0]},
    )


@pytest.mark.parametrize("ext", [".nc", ".zarr"])
@pytest.mark.parametrize(
    "kwargs,offender",
    [
        ({"h": 25.0 + 1.0j}, "h_real"),
        ({"h": 1.0 - 40.0j}, "h_imag"),
        ({"dep": 15000.0}, "dep"),
        ({"dep": -1.0}, "dep"),
    ],
)
def test_out_of_range_values_are_refused(tmpdir, ext, kwargs, offender):
    """Regression: these used to wrap silently, 25 m coming back as -15 m."""
    dset = _one_point(**kwargs)
    with pytest.raises(ValueError, match="outside the range"):
        dset.tide.to_oceantide(str(tmpdir / f"out{ext}"))
    with pytest.raises(ValueError, match=offender):
        dset.tide.to_oceantide(str(tmpdir / f"out{ext}"))


@pytest.mark.parametrize("ext", [".nc", ".zarr"])
def test_in_range_values_round_trip(tmpdir, ext):
    """The check must not reject anything the packing can actually hold."""
    dset = _one_point(h=-19.9 + 19.9j, dep=11999.0)
    filename = str(tmpdir / f"edge{ext}")
    dset.tide.to_oceantide(filename)
    back = read_oceantide(filename)
    assert back.h.values[0, 0, 0].real == pytest.approx(-19.9, abs=1e-3)
    assert back.h.values[0, 0, 0].imag == pytest.approx(19.9, abs=1e-3)
    assert float(back.dep.values[0, 0]) == pytest.approx(11999.0, abs=DEP_ATOL)


@pytest.mark.parametrize("ext", [".nc", ".zarr"])
def test_range_check_can_be_disabled(tmpdir, ext):
    """The opt-out exists for known-safe data; it must not itself raise."""
    _one_point().tide.to_oceantide(str(tmpdir / f"skip{ext}"), check_range=False)


def test_all_missing_variable_is_not_flagged(tmpdir):
    """An entirely masked variable has no values to pack."""
    dset = _one_point()
    dset["h"] = dset.h.where(False)
    dset.tide.to_oceantide(str(tmpdir / "masked.nc"))


@pytest.mark.parametrize("ext", [".nc", ".zarr"])
@pytest.mark.parametrize("keep", [["h"], ["u", "v"], ["h", "dep"], ["h", "u", "v"]])
def test_partial_variables_round_trip(tmpdir, ext, keep):
    """The accessor supports subsets of h, u, v, so the writer must too.

    Regression: to_oceantide indexed self._obj[["dep"]] unconditionally and
    raised KeyError on a dataset without it.
    """
    dset = _one_point()[keep]
    filename = str(tmpdir / f"partial{ext}")
    dset.tide.to_oceantide(filename)
    back = read_oceantide(filename)

    assert set(back.data_vars) == set(keep)
    for varname in keep:
        atol = DEP_ATOL if varname == "dep" else AMP_ATOL
        assert np.allclose(back[varname].values, dset[varname].values, atol=atol)
