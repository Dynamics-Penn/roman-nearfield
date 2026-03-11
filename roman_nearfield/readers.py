"""
readers.py
----------
Functions for reading galaxy, globular cluster, and satellite catalogs
used in the Roman Near-Field Cosmology notebooks.
"""

import numpy as np
import astropy.coordinates
import astropy.table as table
from astropy import units as u

from . import DATA_DIR


def read_galaxies():
    """Read the Karachentsev et al. nearby galaxy catalog within 10 Mpc.

    Parses the pipe-delimited ``nbg.cat`` file, converts the sexagesimal
    RA/Dec strings to an ``astropy.coordinates.SkyCoord`` array, and
    renames columns to a consistent set of names.

    Returns
    -------
    nbg : pandas.DataFrame
        Table of galaxy data with standardized column names including
        ``'Name'``, ``'RAJ2000'``, ``'DEJ2000'``, ``'D (Mpc)'``, and others.
    nbgs : astropy.coordinates.SkyCoord
        Sky coordinates for each galaxy in ``nbg``, in ICRS.

    Examples
    --------
    >>> nbg, nbgs = read_galaxies()
    >>> print(nbg['Name'][:5])
    """
    import pandas as pd

    nbg = pd.read_csv(DATA_DIR / 'nbg.cat', sep="|", usecols=range(1, 13))
    ras = []
    decs = []
    hms_str = ['h', 'm', 's']
    dms_str = ['d', 'm', 's']
    for i in range(len(nbg)):
        ra_str = nbg['RA J2000  '][i].split(' ')
        dec_str = nbg['DEC J2000'][i].split(' ')
        ras.append("".join(x + y for x, y in zip(ra_str, hms_str)))
        decs.append("".join(x + y for x, y in zip(dec_str, dms_str)))
    nbgs = astropy.coordinates.SkyCoord(ras, decs)

    # Rename columns to be consistent with the other galaxy table
    nbg = nbg.rename(columns={
        'Gal. Name    ': 'Name',
        'RA J2000  ':    'RAJ2000',
        'DEC J2000':     'DEJ2000',
        'Bmag ':         'Bmag',
        ' Kmag ':        'Kmag',
        ' D (Mpc)':      'D (Mpc)',
        '  HRV (km/s)   ': 'HRV (km/s)',
        'AbsBmag ':      'BMag',
        'Maj_arcmin':    'Maj. Diam (arcmin)'
    })

    return nbg, nbgs


def read_galaxies_all(dmax=10 * u.Mpc):
    """Read all galaxies within a specified distance from the Karachentsev et al. 2013 catalog.

    Reads the ``nbgs.csv`` file and returns galaxies with distances below
    ``dmax``, together with their sky coordinates.

    Parameters
    ----------
    dmax : astropy.units.Quantity, optional
        Maximum distance for the returned sample. Default is 10 Mpc.

    Returns
    -------
    nbg_all : pandas.DataFrame
        Subset of the full catalog with distances below ``dmax``.
    nbgs_all : astropy.coordinates.SkyCoord
        Sky coordinates for each returned galaxy, in ICRS.

    Examples
    --------
    Load all galaxies within 5 Mpc:

    >>> nbg, nbgs = read_galaxies_all(dmax=5*u.Mpc)
    """
    import pandas as pd

    nbg_all = pd.read_csv(DATA_DIR / 'nbgs.csv')
    ras = nbg_all['RAJ2000'].values * u.deg
    decs = nbg_all['DEJ2000'].values * u.deg
    dsel = nbg_all['D (Mpc)'] < dmax.to(u.Mpc)
    nbgs_all = astropy.coordinates.SkyCoord(ras[dsel], decs[dsel])

    return nbg_all[dsel], nbgs_all


def read_gcs(age_min=0.0, dist_max=500., vb=False):
    """Read the Milky Way globular cluster catalog.

    By default reads from the Local Volume Database (LVDB; Pace 2025).
    Optionally reads from Vasiliev & Baumgardt 2021 (V&B) via VizieR for
    backward compatibility. In both cases the returned table uses a
    consistent set of column names.

    Parameters
    ----------
    age_min : float, optional
        Minimum cluster age in Gyr. Default is 0.0 (all clusters included).
        Currently only applied when ``vb=False``.
    dist_max : float, optional
        Maximum heliocentric distance in kpc. Default is 500 kpc, which
        encompasses the full MW globular cluster system.
    vb : bool, optional
        If True, reads the Vasiliev & Baumgardt 2021 catalog from VizieR
        (requires an internet connection). If False (default), reads the
        current LVDB release.

    Returns
    -------
    gcs : astropy.table.Table
        Table of globular cluster data with columns including ``'name'``,
        ``'ra'``, ``'dec'``, ``'distance'``, ``'pmra'``, ``'pmdec'``,
        and others.

    Examples
    --------
    Load the default LVDB catalog:

    >>> gcs = read_gcs()

    Load the Vasiliev & Baumgardt 2021 catalog for backward compatibility:

    >>> gcs = read_gcs(vb=True)
    """
    if vb:
        from astroquery.vizier import Vizier
        vizier = Vizier()
        CATALOGUE = "J/MNRAS/505/5978/tablea1"
        vizier.ROW_LIMIT = -1
        gcs = vizier.get_catalogs(CATALOGUE)[0]
        vizier_names = ['Name', 'OName', 'RAJ2000', 'DEJ2000', 'pmRA', 'e_pmRA',
                        'pmDE', 'e_pmDE', 'corr', 'plx', 'e_plx', 'Rscale',
                        'Nstar', 'SimbadName']
        lvdb_names = ['name', 'oname', 'ra', 'dec', 'pmra', 'pmra_em', 'pmdec',
                      'pmdec_em', 'corr', 'parallax', 'parallax_em', 'rhalf',
                      'nstar', 'simbad_name']
        gcs.rename_columns(names=vizier_names, new_names=lvdb_names)
        dists = 1.0 / gcs['parallax']  # distances in kpc; some may be negative
        gcs.add_column(dists, name='distance')
        # Mask parallaxes smaller than their uncertainty (including negative)
        gcs['distance'].mask = gcs['parallax'] < gcs['parallax_em']
    else:
        # Load the current release from the LVDB GitHub
        version_number_string = table.Table.read(
            'https://raw.githubusercontent.com/apace7/local_volume_database/main/code/release_version.txt',
            format='ascii.fast_no_header')['col1'][0]
        comb_all = table.Table.read(
            'https://github.com/apace7/local_volume_database/releases/download/'
            + version_number_string + '/comb_all.csv')
        igc = (comb_all['confirmed_star_cluster'] == 1)
        igc &= (comb_all['distance'] < dist_max)
        gcs = comb_all[igc]

    return gcs


def read_satellites():
    """Read the Roman MW Dwarf Satellite Targets catalog.

    Reads the ``Roman_MW_Dwarf_Targets.csv`` file from the data directory.

    Returns
    -------
    nbglm : pandas.DataFrame
        Table of MW satellite/dwarf galaxy targets with columns including
        ``'ra'``, ``'dec'``, and ``'distance'``.

    Examples
    --------
    >>> sats = read_satellites()
    >>> print(len(sats), 'satellites loaded')
    """
    import pandas as pd
    nbglm = pd.read_csv(DATA_DIR / 'Roman_MW_Dwarf_Targets.csv')
    return nbglm


def get_sky_coords(gcs):
    """Extract ICRS sky coordinates from a globular cluster table.

    Handles tables where the RA/Dec columns may or may not carry
    ``astropy.units`` attached, attaching degrees units as needed.

    Parameters
    ----------
    gcs : astropy.table.Table or pandas.DataFrame
        Table of globular cluster data with ``'ra'`` and ``'dec'`` columns
        in decimal degrees.

    Returns
    -------
    gcradec : astropy.coordinates.SkyCoord
        Sky coordinates for each cluster in ICRS.

    Examples
    --------
    >>> gcs = read_gcs()
    >>> coords = get_sky_coords(gcs)
    >>> print(coords.galactic.l[:5])
    """
    from astropy.coordinates import SkyCoord

    if gcs['ra'].unit == u.deg:
        gcradec = SkyCoord(ra=gcs['ra'], dec=gcs['dec'])
    else:
        gcradec = SkyCoord(ra=gcs['ra'] * u.deg, dec=gcs['dec'] * u.deg)

    return gcradec
