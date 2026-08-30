"""Console script for oceantide."""
import click

import oceantide


READERS = {
    "oceantide": "read_oceantide",
    "otis_binary": "read_otis_binary",
    "otis_netcdf": "read_otis_netcdf",
    "otis_atlas_netcdf": "read_otis_atlas_netcdf",
}

reader_option = click.option(
    "--reader",
    "-r",
    type=click.Choice(sorted(READERS)),
    default="oceantide",
    show_default=True,
    help="Format of the input dataset.",
)


def _read(infile: str, reader: str):
    """Open a constituents dataset with the named reader."""
    try:
        read = getattr(oceantide, READERS[reader])
    except AttributeError as err:  # a reader whose optional deps are missing
        raise click.ClickException(
            f"Reader '{READERS[reader]}' is not available: {err}"
        ) from err
    return read(infile)


@click.group()
@click.version_option(oceantide.__version__)
def main():
    """Ocean tide prediction from harmonic constituents."""


@main.command()
@click.argument("infile", type=click.Path(exists=True))
@click.argument("outfile", type=click.Path())
@reader_option
@click.option(
    "--no-check-range",
    is_flag=True,
    help="Skip the check that all values fit the int16 packing.",
)
def convert(infile, outfile, reader, no_check_range):
    """Convert a constituents dataset into the oceantide format.

    INFILE is a file or, for otis_atlas_netcdf, a directory. The output format
    is taken from the OUTFILE extension, either .nc or .zarr.

    \b
    oceantide convert -r otis_netcdf DATA/Model_tpxo9 tpxo9.zarr
    """
    dset = _read(infile, reader)
    dset.tide.to_oceantide(outfile, check_range=not no_check_range)
    click.echo(f"Wrote {outfile}")


@main.command()
@click.argument("infile", type=click.Path(exists=True))
@reader_option
def info(infile, reader):
    """Summarise a constituents dataset."""
    dset = _read(infile, reader)

    cons = [str(c) for c in dset.con.values]
    click.echo(f"constituents  {len(cons)}: {' '.join(cons)}")
    click.echo(f"variables     {' '.join(sorted(dset.data_vars))}")
    for axis in ("lon", "lat"):
        if axis in dset.coords:
            values = dset[axis].values
            click.echo(
                f"{axis:<13} {values.min():.4f} to {values.max():.4f} "
                f"({values.size} points)"
            )
