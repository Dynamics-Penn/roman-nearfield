"""
healpix.py
----------
Utilities for working with HEALPix maps, including coordinate transformations
and loading of survey footprint and Gaia stellar density maps.
"""

import numpy as np
import healpy as hp
from scipy import stats

from . import DATA_DIR


def change_coord(m, coord):
    """Change the coordinate system of a HEALPix map.

    Rotates a HEALPix map (or array of maps) from one coordinate system to
    another by reassigning pixel values according to the transformed angular
    coordinates.

    Parameters
    ----------
    m : numpy.ndarray
        HEALPix map or array of maps to be rotated. The last axis must have
        length ``npix``.
    coord : sequence of str
        Two-character sequence specifying the input and output coordinate
        systems. The first character is the coordinate system of ``m`` and
        the second is the desired output system. Allowed values are
        ``'G'`` (Galactic), ``'E'`` (Ecliptic), and ``'C'`` (Equatorial/ICRS),
        following the HEALPy convention.

    Returns
    -------
    m_rot : numpy.ndarray
        Rotated HEALPix map with the same shape as ``m``.

    Examples
    --------
    Rotate a map from Galactic to Equatorial coordinates:

    >>> m_eq = change_coord(m_gal, ['G', 'C'])

    Notes
    -----
    Implementation adapted from
    https://stackoverflow.com/questions/44443498/how-to-convert-and-save-healpy-map-to-different-coordinate-system
    """
    # Basic HEALPix parameters
    npix = m.shape[-1]
    nside = hp.npix2nside(npix)
    ang = hp.pix2ang(nside, np.arange(npix))

    # Select the coordinate transformation
    rot = hp.Rotator(coord=reversed(coord))

    # Convert the coordinates
    new_ang = rot(*ang)
    new_pix = hp.ang2pix(nside, *new_ang)

    return m[..., new_pix]


def get_footprint(coords='G'):
    """Load and combine the Roman HLWAS wide and medium survey footprint maps.

    Reads the HLWAS WIDE and HLWAS MEDIUM HEALPix maps built from APT outputs,
    combines them into a single map with distinct integer values for each tier,
    and optionally converts to Galactic coordinates.

    Parameters
    ----------
    coords : str, optional
        Output coordinate system. If ``'G'`` (default), the combined map is
        rotated from Equatorial to Galactic coordinates using `change_coord`.
        Any other value leaves the map in its native Equatorial system.

    Returns
    -------
    combined_map : numpy.ndarray
        HEALPix map where pixels covered by HLWAS WIDE have value 1.0,
        pixels covered by HLWAS MEDIUM have value 2.0, and uncovered pixels
        are set to ``hp.UNSEEN``.

    Examples
    --------
    Load the footprint in Galactic coordinates (default):

    >>> fp = get_footprint()

    Load the footprint in Equatorial coordinates:

    >>> fp = get_footprint(coords='C')
    """
    # Read maps with the tiles covered -- these are built from APT's outputs
    HLWAS_wide = hp.read_map(DATA_DIR / 'HLWAS_WIDE_0328_APT994_1024.fits.gz')
    HLWAS_medium = hp.read_map(DATA_DIR / 'HLWAS_MEDIUM_0328_APT994_1024.fits.gz')
    maps = [HLWAS_wide, HLWAS_medium]

    # Combine maps: assign integer tier values
    combined_map = np.zeros_like(HLWAS_wide)
    for i, m in enumerate(maps):
        m[m == hp.UNSEEN] = 0.
        combined_map[m != 0] = (i + 1) * 1.0

    combined_map[combined_map <= 0] = hp.UNSEEN  # Fill back with UNSEEN

    if coords == 'G':
        combined_map = change_coord(combined_map, ['C', 'G'])

    return combined_map


def get_gaia_map(verbose=False):
    """Load and rescale the Gaia stellar density HEALPix map.

    Reads the Gaia G<21 stellar density map, subtracts the minimum value,
    and applies an asinh stretch scaled to 20% of the standard deviation
    to compress the dynamic range for visualization.

    Parameters
    ----------
    verbose : bool, optional
        If True, prints the scaling factor and descriptive statistics of
        the stretched map. Default is False.

    Returns
    -------
    scaled_gmap : numpy.ndarray
        Asinh-stretched HEALPix stellar density map suitable for use as a
        background in sky plots.

    Examples
    --------
    Load the Gaia map and use it as a grayscale background:

    >>> gmap = get_gaia_map()
    >>> sp.draw_hpxmap(gmap, cmap='gray_r', vmin=0.5, vmax=5.0)
    """
    filename = DATA_DIR / 'gaia_stellar_density_map_G_21_gal_nside_128_v0.fits'
    gaiamap = hp.read_map(filename)

    # Rescale the map so it looks good
    gaiamap -= np.min(gaiamap)
    scale = np.std(gaiamap) * 0.2
    if verbose:
        print(scale)
    scaled_gmap = np.asinh(gaiamap / scale)
    if verbose:
        stats.describe(scaled_gmap)

    return scaled_gmap


def get_custom_colorbar():
    """Build a discrete colormap and labels for the Roman survey tier footprint.

    Creates a three-color ``ListedColormap`` mapping the integer tier values
    in the combined footprint map (as produced by `get_footprint`) to distinct
    colors, along with corresponding tier labels.

    Returns
    -------
    cm : matplotlib.colors.ListedColormap
        Colormap with colors for HLWAS WIDE (palegreen), HLWAS MEDIUM
        (limegreen), and GPS (firebrick), in tier order.
    labels : numpy.ndarray of str
        Array of tier label strings: ``['HLWAS WIDE', 'HLWAS MEDIUM', 'GPS']``.
    col_dict : dict
        Mapping from integer tier value (1, 2, 3) to color name string.

    Examples
    --------
    >>> cm, labels, col_dict = get_custom_colorbar()
    >>> sp.draw_hpxmap(fp, cmap=cm, vmin=1, vmax=3)
    """
    from matplotlib.colors import ListedColormap

    col_dict = {1: "palegreen",
                2: "limegreen",
                3: "firebrick"}

    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])
    labels = np.array(["HLWAS WIDE", "HLWAS MEDIUM", "GPS"])

    return cm, labels, col_dict
