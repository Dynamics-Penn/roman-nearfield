"""
footprint_select.py
-------------------
Functions for selecting objects that fall within the Roman survey footprint
or specific sub-regions such as the Galactic Plane Survey (GPS).
"""

import numpy as np
import healpy as hp
from astropy import units as u

from .readers import read_gcs, read_satellites, get_sky_coords


def select_gcs_in_footprint(fp, gcs=None):
    """Select globular clusters that fall within the Roman HLWAS footprint.

    Tests each cluster's position against the HEALPix footprint map and
    splits the result into northern and southern Galactic hemisphere subsets.

    Parameters
    ----------
    fp : numpy.ndarray
        HEALPix footprint map as returned by `get_footprint`, with
        ``hp.UNSEEN`` for uncovered pixels.
    gcs : astropy.table.Table, optional
        Globular cluster catalog as returned by `read_gcs`. If None, loads
        via `read_gcs()`.

    Returns
    -------
    gcs_north : numpy.ndarray of int
        Indices of clusters in the northern Galactic hemisphere (b > 0)
        that fall within the footprint.
    gcs_south : numpy.ndarray of int
        Indices of clusters in the southern Galactic hemisphere (b < 0)
        that fall within the footprint.

    Examples
    --------
    >>> fp = get_footprint()
    >>> gcs_north, gcs_south = select_gcs_in_footprint(fp)
    >>> print(len(gcs_north), 'clusters in northern footprint')
    """
    if gcs is None:
        gcs = read_gcs()
    gc_coords = get_sky_coords(gcs)

    gn = []
    gs = []
    for i, gc in enumerate(gc_coords):
        ipix = hp.ang2pix(hp.get_nside(fp),
                          gc.galactic.l.value, gc.galactic.b.value,
                          lonlat=True)
        if not fp[ipix] == hp.UNSEEN:
            if gc.galactic.b.value > 0:
                gn.append(i)
            else:
                gs.append(i)

    gcs_north = np.array(gn)
    gcs_south = np.array(gs)

    print(len(gcs_north), 'gcs in northern footprint')
    print(len(gcs_south), 'gcs in southern footprint')

    return gcs_north, gcs_south


def select_gcs_in_gps_main(gcs=None):
    """Select globular clusters within the Roman Galactic Plane Survey (GPS) footprint.

    Tests each cluster's Galactic coordinates against the approximate
    polygon defining the main GPS region.

    Parameters
    ----------
    gcs : astropy.table.Table, optional
        Globular cluster catalog as returned by `read_gcs`. If None, loads
        via `read_gcs()`.

    Returns
    -------
    gcs_gps : numpy.ndarray of int
        Indices of clusters that fall within the GPS polygon.

    Examples
    --------
    >>> gcs_gps = select_gcs_in_gps_main()
    >>> print(len(gcs_gps), 'clusters in GPS')
    """
    import matplotlib.patches as mpatches

    if gcs is None:
        gcs = read_gcs()
    gc_coords = get_sky_coords(gcs)

    ggps = []
    l_RGPS = [50.1, 30, 30, 26.5, 26.5, 10, 10, -10, -10, -67, -67,
              -79, -79, -67, -67, -10, -10, 10, 10, 50.1, 50.1] * u.deg
    b_RGPS = [2, 2, 4.5, 4.5, 2, 2, 6, 6, 2, 2, 2.0, 2.0, -2.5,
              -2.5, -2, -2, -6, -6, -2, -2, 2] * u.deg
    GPS_coords = np.array([l_RGPS, b_RGPS]).T
    GPS_patch = mpatches.Polygon(GPS_coords, closed=True)

    for i, gc in enumerate(gc_coords):
        if GPS_patch.contains_point(
                (gc.galactic.l.wrap_at(180 * u.deg).value, gc.galactic.b.value)):
            ggps.append(i)

    gcs_gps = np.array(ggps)
    print(len(gcs_gps), 'gcs in GPS')

    return gcs_gps


def select_galaxies_in_footprint(fp, nbg_coords, nbgs):
    """Select nearby galaxies that fall within the Roman HLWAS footprint.

    Tests each galaxy's Galactic coordinates against the footprint map,
    splits results by Galactic hemisphere, and excludes MW satellites
    (galaxies closer than 0.3 Mpc).

    Parameters
    ----------
    fp : numpy.ndarray
        HEALPix footprint map as returned by `get_footprint`.
    nbg_coords : astropy.coordinates.SkyCoord
        Sky coordinates of the galaxies to test.
    nbgs : pandas.DataFrame or astropy.table.Table
        Table of galaxy data with a ``'D (Mpc)'`` column used to identify
        and exclude MW satellites.

    Returns
    -------
    nbgs_in_hlwas_north : list of int
        Sorted indices of non-satellite galaxies in the northern footprint.
    nbgs_in_hlwas_south : list of int
        Sorted indices of non-satellite galaxies in the southern footprint.

    Examples
    --------
    >>> fp = get_footprint()
    >>> nbg, nbg_coords = read_galaxies()
    >>> north, south = select_galaxies_in_footprint(fp, nbg_coords, nbg)
    """
    gn = []
    gs = []
    for i, sat in enumerate(nbg_coords):
        ipix = hp.ang2pix(hp.get_nside(fp),
                          sat.galactic.l.value, sat.galactic.b.value,
                          lonlat=True)
        if not fp[ipix] == hp.UNSEEN:
            if sat.galactic.b.value > 0:
                gn.append(i)
            else:
                gs.append(i)

    in_mw = set(np.where(nbgs['D (Mpc)'] < 0.3)[0].astype(int))

    nbgs_in_hlwas_north = sorted(set(gn).difference(in_mw))
    nbgs_in_hlwas_south = sorted(set(gs).difference(in_mw))

    print(len(nbgs_in_hlwas_north), 'galaxies in northern footprint')
    print(len(nbgs_in_hlwas_south), 'galaxies in southern footprint')

    return nbgs_in_hlwas_north, nbgs_in_hlwas_south


def select_satellites_in_footprint(fp, nbglm=None):
    """Select MW satellite galaxies that fall within the Roman HLWAS footprint.

    Tests each satellite's sky position against the footprint map and splits
    the result into northern and southern Galactic hemisphere subsets.

    Parameters
    ----------
    fp : numpy.ndarray
        HEALPix footprint map as returned by `get_footprint`.
    nbglm : pandas.DataFrame, optional
        Table of satellite data with ``'ra'`` and ``'dec'`` columns in
        decimal degrees. If None, loads via `read_satellites()`.

    Returns
    -------
    nbglm_in_hlwas_north : list of int
        Indices of satellites in the northern Galactic hemisphere footprint.
    nbglm_in_hlwas_south : list of int
        Indices of satellites in the southern Galactic hemisphere footprint.

    Examples
    --------
    >>> fp = get_footprint()
    >>> north, south = select_satellites_in_footprint(fp)
    >>> print(len(north), 'satellites in northern footprint')
    """
    import astropy.coordinates

    if nbglm is None:
        nbglm = read_satellites()

    nbglm_in_hlwas_north = []
    nbglm_in_hlwas_south = []

    for i in range(len(nbglm)):
        coord = astropy.coordinates.SkyCoord(
            nbglm['ra'][i] * u.deg, nbglm['dec'][i] * u.deg)
        ipix = hp.ang2pix(hp.get_nside(fp),
                          coord.galactic.l.value, coord.galactic.b.value,
                          lonlat=True)
        if not fp[ipix] == hp.UNSEEN:
            if coord.galactic.b.value > 0:
                nbglm_in_hlwas_north.append(i)
            else:
                nbglm_in_hlwas_south.append(i)

    print(len(nbglm_in_hlwas_north), 'satellites in northern footprint')
    print(len(nbglm_in_hlwas_south), 'satellites in southern footprint')

    return nbglm_in_hlwas_north, nbglm_in_hlwas_south


def get_gptds_fields():
    """Return the bounding boxes and vertices for Roman GPTDS target fields.

    Defines six Galactic Plane Time Domain Survey (GPTDS) fields by name,
    Galactic longitude/latitude ranges, and corner vertices, as provided
    by Bob Benjamin.

    Returns
    -------
    gptds : dict
        Dictionary keyed by field name, where each value is a dict with
        keys ``'lmin'``, ``'lmax'``, ``'bmin'``, ``'bmax'``, and
        ``'vertices'`` (a ``(2, 4)`` astropy Quantity array in degrees
        giving the field corners).

    Examples
    --------
    >>> fields = get_gptds_fields()
    >>> for name, field in fields.items():
    ...     print(name, field['lmin'], 'to', field['lmax'])
    """
    tds_name = [
        'TDS_Carina',
        'TDS_NGC_6334_6357',
        'TDS_Galactic_Center_Q4',
        'TDS_Galactic_Center_Q1',
        'TDS_Serpens_South_W40',
        'TDS_W43']
    tdsL = [
        [-73.8, -71.2],
        [-9.2, -6.6],
        [-2.8, -0.2],
        [0.1, 2.6],
        [26.5, 30.0],
        [29.3, 31.9]]
    tdsB = [
        [-1.0, -0.2],
        [0.3, 1.1],
        [-0.5, 0.3],
        [-0.5, 0.3],
        [2.0, 4.5],
        [-0.5, 0.3]]

    gptds = {}
    for i, k in enumerate(tds_name):
        lmm = tdsL[i]
        bmm = tdsB[i]
        gptds[k] = {
            'lmin': lmm[0],
            'lmax': lmm[1],
            'bmin': bmm[0],
            'bmax': bmm[1],
            'vertices': np.array([[lmm[0], lmm[0], lmm[1], lmm[1]],
                                   [bmm[0], bmm[1], bmm[1], bmm[0]]]) * u.deg
        }

    return gptds
