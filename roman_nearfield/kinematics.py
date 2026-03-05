"""
kinematics.py
-------------
Functions for converting between proper motions, distances, and transverse
velocities, and for computing proper motion detection limits as a function
of observing baseline and distance.
"""

from astropy import units as u

k = 4.74 *((u.km/u.s)/(u.pc)/(u.arcsec/u.yr))  ## conversion const for PMs


def mu(dist, vt):
    """Convert transverse velocity and distance to proper motion.

    Computes the proper motion corresponding to a transverse velocity ``vt``
    at heliocentric distance ``dist``. Accepts plain numbers (assumed to be
    in kpc and km/s respectively) or astropy Quantities with any compatible
    units.

    Parameters
    ----------
    dist : float, array-like, or astropy.units.Quantity
        Heliocentric distance. Plain numbers are assumed to be in kpc.
    vt : float, array-like, or astropy.units.Quantity
        Transverse velocity. Plain numbers are assumed to be in km/s.

    Returns
    -------
    pm : astropy.units.Quantity
        Proper motion in microarcseconds per year (uas/yr).

    Examples
    --------
    Using plain floats (kpc and km/s assumed):

    >>> mu(50., 100.)   # distance 50 kpc, velocity 100 km/s
    <Quantity ... uas / yr>

    Using astropy Quantities:

    >>> mu(50.*u.kpc, 100.*u.km/u.s)
    <Quantity ... uas / yr>
    """

    if not isinstance(dist, u.Quantity):
        dist = dist * u.kpc
    if not isinstance(vt, u.Quantity):
        vt = vt * (u.km / u.s)
    return (vt / (k * dist)).to(u.uas/u.yr)


def vt(dist, pm):
    """Convert proper motion and distance to transverse velocity.

    Computes the transverse velocity corresponding to a proper motion ``pm``
    at heliocentric distance ``dist``. Accepts plain numbers (assumed to be
    in kpc and uas/yr respectively) or astropy Quantities with any compatible
    units.

    Parameters
    ----------
    dist : float, array-like, or astropy.units.Quantity
        Heliocentric distance. Plain numbers are assumed to be in kpc.
    pm : float, array-like, or astropy.units.Quantity
        Proper motion. Plain numbers are assumed to be in uas/yr.

    Returns
    -------
    vel : astropy.units.Quantity
        Transverse velocity in km/s.

    Examples
    --------
    Using plain floats (kpc and uas/yr assumed):

    >>> vt(50., 10.)   # distance 50 kpc, PM 10 uas/yr
    <Quantity ... km / s>

    Using astropy Quantities:

    >>> vt(50.*u.kpc, 10.*u.uas/u.yr)
    <Quantity ... km / s>
    """
    if not isinstance(dist, u.Quantity):
        dist = dist * u.kpc
    if not isinstance(pm, u.Quantity):
        pm = pm * (u.uas / u.yr)
    return (k * dist * pm).to(u.km/u.s)


def velocity_limit(baseline, dist, loc_hls=100. * u.uas, Nobs_hls=5):
    """Compute the transverse velocity detection limit for a given observing baseline.

    Calculates the minimum detectable transverse velocity at distance ``dist``
    given a proper motion uncertainty of ``loc_hls / sqrt(Nobs_hls)`` per
    measurement over the provided time baseline.

    Parameters
    ----------
    baseline : astropy.units.Quantity
        Total time baseline of the observations (e.g. in years).
    dist : float or astropy.units.Quantity
        Heliocentric distance. Plain numbers are assumed to be in kpc.
    loc_hls : astropy.units.Quantity, optional
        Single-epoch astrometric uncertainty. Default is 100 uas.
    Nobs_hls : int, optional
        Number of observations. Default is 5.

    Returns
    -------
    vel : astropy.units.Quantity
        Minimum detectable transverse velocity in km/s.

    Examples
    --------
    >>> velocity_limit(5.*u.yr, 100.*u.kpc)
    <Quantity ... km / s>
    """
    pm_limit = loc_hls / baseline / Nobs_hls ** 0.5
    return vt(dist, pm_limit).to(u.km / u.s)


def velocity_limit_year(year, dist, hls_year=2027 * u.yr,
                        loc_hls=100. * u.uas, Nobs_hls=5):
    """Compute the transverse velocity detection limit as a function of catalog epoch.

    Computes the minimum detectable transverse velocity using the baseline
    from a reference catalog epoch (``year``) to the Roman HLS start year
    (``hls_year``).

    Parameters
    ----------
    year : astropy.units.Quantity
        Epoch of the reference catalog in years (e.g. a Gaia data release
        epoch).
    dist : float or astropy.units.Quantity
        Heliocentric distance. Plain numbers are assumed to be in kpc.
    hls_year : astropy.units.Quantity, optional
        Start year of the Roman HLS. Default is 2027 yr.
    loc_hls : astropy.units.Quantity, optional
        Single-epoch astrometric uncertainty. Default is 100 uas.
    Nobs_hls : int, optional
        Number of Roman HLS observations. Default is 5.

    Returns
    -------
    vel : astropy.units.Quantity
        Minimum detectable transverse velocity in km/s.

    Examples
    --------
    >>> velocity_limit_year(2016.*u.yr, 50.*u.kpc)
    <Quantity ... km / s>
    """
    pm_limit = loc_hls / (hls_year - year) / Nobs_hls ** 0.5
    return vt(dist, pm_limit).to(u.km / u.s)


def get_year(s):
    """Extract the year as a float from a string beginning with a 4-digit year.

    Parameters
    ----------
    s : str
        String whose first four characters are a calendar year (e.g. a
        date string like ``'2016.5'`` or ``'2025-03-01'``).

    Returns
    -------
    year : float
        The year parsed from the first four characters of ``s``.

    Examples
    --------
    >>> get_year('2016.5')
    2016.0
    >>> get_year('2025-03-01')
    2025.0
    """
    return float(s[:4])
