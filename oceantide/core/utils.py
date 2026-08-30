"""Core tools for tidal analysis and prediction."""
import warnings
from pathlib import Path
import yaml
import numpy as np
import xarray as xr

from oceantide.constituents import SHALLOW, V0U


HERE = Path(__file__).parent


def nodal(time: np.ndarray, con: np.ndarray):
    """Nodal correction.

    Parameters
    ----------
    time (1darray)
        Time given as the number of days since 01 Jan 1992 + 48622, or equivalently
        the number of days since 17 Nov 1858.
    con (1darray)
        Constituents to compute.

    Returns
    -------
    pu (1darray)
        Nodal correction pu.
    pf (1darray)
        Nodal correction pf.
    v0u (1darray)
        Nodal correction v0u.

    """
    rad = np.pi / 180.0
    s, h, p, omega = astrol(time)
    sinn = np.sin(omega * rad)
    cosn = np.cos(omega * rad)
    sin2n = np.sin(2 * omega * rad)
    cos2n = np.cos(2 * omega * rad)
    sin3n = np.sin(3 * omega * rad)

    def polar(real, imag):
        """Nodal factor and angle as the modulus and argument of one term.

        The nodal correction of these constituents is a single complex
        modulation term: f is its modulus and u its argument. Writing the two
        out separately means spelling the coefficients twice, which is how they
        drift apart -- K1's f carried 0.01554 where its own u used 0.1554, a
        dropped digit worth 1.2% of K1 amplitude at the nodal extremes.

        """
        z = real + 1j * imag
        return {"f": np.abs(z), "u": np.angle(z) / rad}

    # Semidiurnal lunar group, Schureman via the OTIS `nodal.f` coefficients.
    m2 = polar(
        1.0 - 0.03731 * cosn + 0.00052 * cos2n,
        -0.03731 * sinn + 0.00052 * sin2n,
    )

    ndict = {
        "M2": m2,
        "S2": {"f": 1.0, "u": 0.0},
        "K1": polar(
            1.0 + 0.1158 * cosn - 0.0029 * cos2n,
            -0.1554 * sinn + 0.0029 * sin2n,
        ),
        # O1 and Q1 keep OTIS's hand-written angle series rather than the
        # argument of their own f term; the two are not interchangeable, and
        # Q1 genuinely uses 0.188 in f against 0.189 in u.
        "O1": {
            "f": np.hypot(
                1.0 + 0.189 * cosn - 0.0058 * cos2n,
                0.189 * sinn - 0.0058 * sin2n,
            ),
            "u": 10.8 * sinn - 1.3 * sin2n + 0.2 * sin3n,
        },
        "N2": m2,
        "P1": {"f": 1.0, "u": 0.0},
        "K2": polar(
            1.0 + 0.2852 * cosn + 0.0324 * cos2n,
            -0.3108 * sinn - 0.0324 * sin2n,
        ),
        "Q1": {
            "f": np.hypot(1.0 + 0.188 * cosn, 0.188 * sinn),
            "u": np.arctan(0.189 * sinn / (1.0 + 0.189 * cosn)) / rad,
        },
        # Long period. Mm and Mf carry the standard Schureman factors; Mf's
        # nodal modulation is large -- f swings by about 40% over the 18.6 year
        # cycle -- so leaving it at unity is a material error.
        "MM": {"f": 1.0 - 0.1300 * cosn + 0.0013 * cos2n, "u": 0.0},
        "MF": {
            "f": 1.043 + 0.414 * cosn,
            "u": -23.7 * sinn + 2.7 * sin2n - 0.4 * sin3n,
        },
        # Solar constituents have no lunar nodal modulation. Listed explicitly
        # rather than left to the default so that unity is visibly deliberate.
        "SA": {"f": 1.0, "u": 0.0},
        "SSA": {"f": 1.0, "u": 0.0},
        "S1": {"f": 1.0, "u": 0.0},
        "T2": {"f": 1.0, "u": 0.0},
    }

    # 2N2, mu2 and nu2 sit in M2's nodal group and share its factors.
    for name in ("2N2", "MU2", "NU2"):
        ndict[name] = ndict["M2"]

    # Compound constituents are products of their parents. A compound tide
    # arises from a product of the parent potentials in the nonlinear terms, so
    # its amplitude scales with the product of the parents' amplitudes however
    # their arguments combine: f takes the magnitude of each exponent, and only
    # the phase correction u carries the sign. Foreman (1977) and t_tide do the
    # same -- t_vuf.m raises f to abs(coef) and sums u with coef.
    #
    # This matters only for MSF (S2 - M2), the one entry with a negative
    # exponent, where the signed form gave 1/f(M2) instead of f(M2) -- a 7.8%
    # error on MSF amplitude across the nodal cycle.
    for name, parents in SHALLOW.items():
        f = 1.0
        u = 0.0
        for parent, power in parents.items():
            f = f * ndict[parent]["f"] ** abs(power)
            u = u + ndict[parent]["u"] * power
        ndict[name] = {"f": f, "u": u}

    ncon = len(con)
    pu = np.zeros(ncon)
    pf = np.ones(ncon)
    v0u = np.zeros(ncon)

    for ii, vv in enumerate(con):
        if vv in ndict:
            pu[ii] = ndict[vv]["u"] * rad
            pf[ii] = ndict[vv]["f"]
        if vv in V0U:
            v0u[ii] = V0U[vv]
        else:
            # Falling through with v0u = 0 references the constituent to an
            # arbitrary phase, which produces a plausible looking but wrong
            # tide. Say so, and contribute nothing, rather than let it pass.
            warnings.warn(
                f"No equilibrium argument for constituent {vv!r}; it would be "
                f"predicted with an arbitrary phase reference and is dropped. "
                f"Add it to oceantide.constituents.V0U to support it.",
                stacklevel=2,
            )
            pf[ii] = 0.0

    return pu, pf, v0u


def astrol(time: np.ndarray):
    """Mean astronomical longitudes  s, h, p, N.

    Parameters
    ----------
    time (1darray)
        Time given as the number of days since 01 Jan 1992 + 48622, or equivalently
        the number of days since 17 Nov 1858.

    Returns
    -------
    s (1darray)
        The mean longitude of moon.
    h (1darray)
        The mean longitude of sun.
    p (1darray)
        The mean longitude of lunar perigee.
    N (1darray)
        The mean longitude of ascending lunar node.

    """
    T = time - 51544.4993
    s = (218.3164 + 13.17639648 * T) % 360
    h = (280.4661 + 0.98564736 * T) % 360
    p = (83.3535 + 0.11140353 * T) % 360
    N = (125.0445 - 0.05295377 * T) % 360
    return s, h, p, N


def arakawa_grid(
    nx: int, ny: int, x0: float, x1: float, y0: float, y1: float, variable: str
):
    """Arakawa grid coordinates for variable.

    Parameters
    ----------
    nx (int)
        Number of grid point along the x direction.
    ny (int)
        Number of grid point along the y direction.
    x0 (float)
        Left limit (degrees).
    x1 (float)
        Right limit (degrees).
    y0 (float)
        Bottom limit (degrees).
    y1 (float)
        Left limit (degrees).
    variable (str)
        Model variable, one of 'h', 'u', 'v'.

    Returns
    -------
    lon (1darray)
        Longitude coordinates in Arakawa grid for variable.
    lat (1darray)
        Latitude coordinates in Arakawa for variable.

    """
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    xcorner = np.clip(np.arange(x0, x1, dx), x0, x1)[0:nx]
    ycorner = np.clip(np.arange(y0, y1, dy), y0, y1)[0:ny]
    if variable == "h":
        lon = xcorner + dx / 2
        lat = ycorner + dy / 2
    elif variable == "u":
        lon = xcorner
        lat = ycorner + dy / 2
    elif variable == "v":
        lon = xcorner + dx / 2
        lat = ycorner
    else:
        raise ValueError(f"'variable' must be one of 'h', 'u', 'v', got {variable}")
    return lon, lat


def set_attributes(dset: xr.Dataset, dataset_type: str):
    """Set variable attributes for variables from given dataset type.

    Parameters
    ----------
    dset (Dataset)
        Dataset to set attributes from.docs
    dataset_type (str)
        Dataset type, e.g., 'otis'.

    """
    all_attrs = yaml.load(open(HERE / "attributes.yml"), Loader=yaml.Loader)
    if isinstance(dset, xr.Dataset):
        for varname, darr in dset.variables.items():
            darr.attrs.update(all_attrs.get(dataset_type).get(varname, {}))
    elif isinstance(dset, xr.DataArray):
        dset.attrs.update(all_attrs.get(dataset_type).get(dset.name, {}))


def compute_scale_and_offset(vmin: float, vmax: float, nbit: int = 16):
    """Returns the scale_factor and add_offset for packing data.

    Parameters
    ----------
    vmin (float)
        The minumum value in the array to pack.
    vmax (float)
        The maximum value in the array to pack.
    nbit (int)
        The number of bits into which you wish to pack.

    Returns
    -------
    scale_factor (float)
        The scale factor value to use for packing array.
    add_offset (float)
        The add_offset value to use for packing array.

    """
    if vmax == vmin:
        scale_factor = 1
    else:
        scale_factor = (vmax - vmin) / (2**nbit - 1)
    add_offset = vmin + 2 ** (nbit - 1) * scale_factor
    return scale_factor, add_offset
