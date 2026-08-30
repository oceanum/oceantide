#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the oceantide command line interface."""

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import oceantide
from oceantide import cli


FILES_DIR = Path(__file__).parent / "test_files"


@pytest.fixture
def runner():
    return CliRunner()


def test_help_lists_the_commands(runner):
    result = runner.invoke(cli.main, ["--help"])
    assert result.exit_code == 0
    assert "convert" in result.output
    assert "info" in result.output


def test_version(runner):
    result = runner.invoke(cli.main, ["--version"])
    assert result.exit_code == 0
    assert oceantide.__version__ in result.output


def test_no_arguments_shows_help(runner):
    """Regression: this used to print 'Replace this message by putting your
    code into oceantide.cli.main'."""
    result = runner.invoke(cli.main, [])
    assert "Replace this message" not in result.output
    assert "Usage:" in result.output


def test_info(runner):
    result = runner.invoke(cli.main, ["info", str(FILES_DIR / "oceantide.zarr")])
    assert result.exit_code == 0, result.output
    assert "constituents" in result.output
    assert "M2" in result.output
    assert "lon" in result.output and "lat" in result.output


@pytest.mark.parametrize("ext", [".nc", ".zarr"])
def test_convert_oceantide(runner, tmpdir, ext):
    outfile = str(tmpdir / f"out{ext}")
    result = runner.invoke(
        cli.main, ["convert", str(FILES_DIR / "oceantide.zarr"), outfile]
    )
    assert result.exit_code == 0, result.output

    original = oceantide.read_oceantide(FILES_DIR / "oceantide.zarr")
    converted = oceantide.read_oceantide(outfile)
    assert set(converted.data_vars) == set(original.data_vars)
    assert np.allclose(
        converted.h.values, original.h.values, atol=1e-3, equal_nan=True
    )


def test_convert_otis_netcdf(runner, tmpdir):
    outfile = str(tmpdir / "otis.zarr")
    result = runner.invoke(
        cli.main,
        [
            "convert",
            "-r",
            "otis_netcdf",
            str(FILES_DIR / "otis_netcdf/Model_test"),
            outfile,
        ],
    )
    assert result.exit_code == 0, result.output
    assert oceantide.read_oceantide(outfile).con.size > 0


def test_convert_rejects_an_unknown_extension(runner, tmpdir):
    result = runner.invoke(
        cli.main, ["convert", str(FILES_DIR / "oceantide.zarr"), str(tmpdir / "x.txt")]
    )
    assert result.exit_code != 0


def test_missing_input_is_reported(runner, tmpdir):
    result = runner.invoke(cli.main, ["info", str(tmpdir / "nope.zarr")])
    assert result.exit_code != 0
    assert "does not exist" in result.output
