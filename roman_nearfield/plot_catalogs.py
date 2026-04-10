"""
plot_catalogs.py
----------------
Functions for plotting astronomical catalogs (galaxies, globular clusters,
stellar streams) onto skyproj map objects.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy import units as u

from .readers import read_galaxies, read_galaxies_all, read_gcs, get_sky_coords


def spread_clusters_kernel(points, radius, k=0.1, max_iter=100, tol=1e-6,
                           min_repulsion=0.01):
    """Spread apart annotation label positions that are too close together.

    Uses a kernel-based repulsion approach: points within ``radius`` of each
    other exert a repulsive force scaled by an inverse-square kernel and by
    local point density. Iterates until convergence or ``max_iter`` steps.

    Parameters
    ----------
    points : numpy.ndarray
        Array of shape ``(N, 2)`` containing the x, y positions to adjust.
    radius : float
        Points within this distance are considered too close and will be
        repelled.
    k : float, optional
        Strength of the repulsion force. Default is 0.1.
    max_iter : int, optional
        Maximum number of iterations. Default is 100.
    tol : float, optional
        Convergence tolerance on maximum per-step displacement. Default is 1e-6.
    min_repulsion : float, optional
        Minimum repulsion weight applied per neighbor, to avoid numerical
        issues at very small separations. Default is 0.01.

    Returns
    -------
    points : numpy.ndarray
        Adjusted point positions with the same shape as the input.

    Notes
    -----
    This function was originally generated with DeepSeek and has not been
    exhaustively tested. Use with care for unusual point configurations.

    Examples
    --------
    >>> adjusted = spread_clusters_kernel(label_positions, radius=10, k=5)
    """
    from scipy.spatial import cKDTree

    points = points.copy().astype(float)
    n_points = len(points)
    tree = cKDTree(points)

    for iteration in range(max_iter):
        max_displacement = 0
        new_points = points.copy()

        for i in range(n_points):
            indices = tree.query_ball_point(points[i], radius)

            if len(indices) <= 1:
                continue

            vectors = points[indices] - points[i]
            distances = np.linalg.norm(vectors, axis=1)

            mask = distances > 0
            vectors = vectors[mask]
            distances = distances[mask]

            if len(distances) == 0:
                continue

            # Inverse-square kernel: stronger repulsion for closer points
            repulsion_weights = k * (1.0 - distances / radius) ** 2

            # Density factor: more neighbors = stronger overall repulsion
            density_factor = len(distances) / (np.pi * radius ** 2)
            repulsion_weights *= (1 + density_factor)

            repulsion_weights = np.maximum(repulsion_weights, min_repulsion)

            norm_vectors = vectors / distances[:, np.newaxis]
            displacement = np.sum(norm_vectors * repulsion_weights[:, np.newaxis], axis=0)

            new_points[i] += displacement
            max_displacement = max(max_displacement, np.linalg.norm(displacement))

        points = new_points
        tree = cKDTree(points)

        if max_displacement < tol:
            break

    return points


def plot_galaxies(sp, gals=None, glist=None, annotate=False, glabels=None,
                  radec=False, spread_scale=4, angles=None, radii=None,
                  annotation_style='auto'):
    """Plot galaxies from the Karachentsev et al. nearby galaxy catalog.

    Plots galaxies on a skyproj map, colored by distance on a logarithmic
    scale. Optionally annotates a subset of galaxies with labels and
    leader lines, with several strategies available for placing the labels.

    Parameters
    ----------
    sp : skyproj.Skyproj
        Map object on which to plot the galaxies.
    gals : pandas.DataFrame or astropy.table.Table, optional
        Table of galaxy data. If None, loads the default catalog via
        `read_galaxies`. Accepts tables with either ``'RAJ2000'``/
        ``'DEJ2000'`` columns (HMS/DMS string format) or ``'ra'``/``'dec'``
        columns (decimal degrees).
    glist : array-like of int, optional
        Indices into ``gals`` selecting the subset of galaxies to plot and
        optionally annotate. If None, all galaxies are plotted.
    annotate : bool or dict, optional
        If True, annotates the galaxies selected by ``glist`` with leader
        lines. If a dict, its contents are passed as keyword arguments to
        ``ax.annotate``. Default is False.
    glabels : array-like of str, optional
        Labels to use for annotation. If None, galaxies are numbered 1 to N.
    radec : bool, optional
        If True, plots in RA/Dec coordinates. If False (default), plots in
        Galactic coordinates.
    spread_scale : float, optional
        Scaling factor used by the ``'spread'`` annotation strategy to
        separate crowded labels. Default is 4.
    angles : array-like of float, optional
        Angles in radians for manual label placement. Only used when
        ``annotation_style='manual'``.
    radii : array-like of float, optional
        Radii in points for manual label placement. Only used when
        ``annotation_style='manual'``.
    annotation_style : {'auto', 'spread', 'random', 'manual'}, optional
        Strategy for placing annotation labels. ``'spread'`` uses a
        pole-based layout refined by `spread_clusters_kernel` to separate
        crowded labels; ``'random'`` uses random angular offsets;
        ``'manual'`` uses the ``angles`` and ``radii`` arrays. ``'auto'``
        (default) infers the style from the format of ``gals``, replicating
        the per-notebook behavior prior to this package.

    Returns
    -------
    cbh : matplotlib.cm.ScalarMappable
        Scalar mappable containing the colormap and normalization, suitable
        for passing to ``plt.colorbar()``.

    Examples
    --------
    Plot all galaxies on a predefined skyproj in Galactic coordinates:

    >>> cbh = plot_galaxies(sp)

    Plot a labeled subset using the spread annotation strategy:

    >>> cbh = plot_galaxies(sp, gals=my_table, glist=[0, 3, 7],
    ...                     annotate=True, glabels=['M31', 'M33', 'LMC'],
    ...                     annotation_style='spread')
    """
    from matplotlib.colors import LogNorm
    from astropy.coordinates import SkyCoord

    if gals is None:
        nbg, nbgs = read_galaxies(kara=True)
        dists = nbg['D (Mpc)'] * u.Mpc
        norm = LogNorm(vmin=0.2, vmax=10)
        r = 80
    else:
        if 'RAJ2000' in gals.keys():  #nearby galaxies
            ras = []
            decs = []
            hms_str = ['h', 'm', 's']
            dms_str = ['d', 'm', 's']
            for i in range(len(gals)):
                ra_str = gals['RAJ2000'][i].split(' ')
                dec_str = gals['DEJ2000'][i].split(' ')
                ras.append("".join(x + y for x, y in zip(ra_str, hms_str)))
                decs.append("".join(x + y for x, y in zip(dec_str, dms_str)))
            nbgs = SkyCoord(ras, decs)
            dists = gals['D (Mpc)']
            norm = LogNorm(vmin=0.2, vmax=10)
            r = 80
        else:  #milky way satellites
            nbgs = SkyCoord(gals['ra'].values * u.deg, gals['dec'].values * u.deg)
            dists = ((gals['distance'].values * u.kpc)).value #.to(u.Mpc)).value
            norm = LogNorm(vmin=15, vmax=350)
            r = 60

    cmap = plt.get_cmap('cividis')

    if glist is None:
        if radec:
            cbh = sp.scatter(nbgs.ra.value, nbgs.dec.value, c=dists,
                             edgecolors="black", cmap=cmap, norm=norm, s=40, zorder=2.4)
        else:
            cbh = sp.scatter(nbgs.galactic.l.value, nbgs.galactic.b.value, c=dists,
                             edgecolors="black", cmap=cmap, norm=norm, s=40, zorder=2.4)
    else:
        if radec:
            cbh = sp.scatter(nbgs.ra.value[glist], nbgs.dec.value[glist],
                             c=dists[glist], edgecolors="black", cmap=cmap,
                             norm=norm, s=40, zorder=2.4)
        else:
            cbh = sp.scatter(nbgs.galactic.l.value[glist],
                             nbgs.galactic.b.value[glist],
                             c=dists[glist], edgecolors="black", cmap=cmap,
                             norm=norm, s=40, zorder=2.4)

        if annotate:
            if glabels is None:
                glabels = np.arange(1, len(glist) + 1).astype('str')

            xy_positions = []
            for i, l in enumerate(glabels):
                if radec:
                    x, y = sp.proj(nbgs[glist[i]].ra.deg, nbgs[glist[i]].dec.deg)
                else:
                    x, y = sp.proj(nbgs[glist[i]].galactic.l.deg,
                                   nbgs[glist[i]].galactic.b.deg)
                xy_positions.append([x[0], y[0]])
            xy_positions = np.array(xy_positions)

            # Determine annotation strategy
            if annotation_style == 'manual' or (
                    annotation_style == 'auto' and radii is not None and angles is not None):
                xyt_positions = np.array(
                    xy_positions.T + [radii * np.cos(angles), radii * np.sin(angles)]).T

            elif annotation_style == 'random' or (
                    annotation_style == 'auto' and gals is not None and 'ra' in gals.keys()):
                theta = np.random.random_sample(size=len(xy_positions)) * 2.0 * np.pi
                xyt_positions = np.array(
                    xy_positions.T + [r * np.cos(theta), r * np.sin(theta)]).T

            else:
                # 'spread' or 'auto' with RAJ2000-style table
                xy_center = np.mean(xy_positions, axis=0)
                theta_pole = np.arctan2(*tuple((xy_positions - xy_center).T))
                theta = theta_pole + np.pi / 2
                xyt_positions = (xy_positions.T + [r * np.cos(theta), r * np.sin(theta)]).T
                xyt_positions = [spread_scale, 1] * spread_clusters_kernel(
                    xyt_positions * [1 / spread_scale, 1], radius=10, k=10)

            xyt_positions -= xy_positions
            for i, l in enumerate(glabels):
                x, y = xy_positions[i]
                xt, yt = xyt_positions[i]
                sp.ax.annotate(l.strip(' '), (x, y), (xt, yt),
                               textcoords='offset points',
                               arrowprops=dict(arrowstyle='-', color='k'),
                               va='top', ha='left',
                               **(annotate if isinstance(annotate, dict) else {}))

    return cbh


def add_extra_galaxies(sp, cbh=None, nbg_coords=None, dists=None, glist=None,
                       radec=False, dmax=10 * u.Mpc):
    """Add unlabeled background galaxies to an existing skyproj map.

    Plots additional galaxies (without annotation) out to distance ``dmax``.
    If coordinates and distances are supplied they are plotted with a shared
    color scale; otherwise they are plotted as small gray dots.

    Parameters
    ----------
    sp : skyproj.Skyproj
        Map object on which to plot the galaxies.
    cbh : matplotlib.cm.ScalarMappable, optional
        Scalar mappable from a prior ``plot_galaxies`` call. If provided
        along with ``dists``, the same colormap and normalization are reused.
        If None, galaxies are plotted as gray dots.
    nbg_coords : astropy.coordinates.SkyCoord, optional
        Sky coordinates of the galaxies to plot. If None, loads via
        `read_galaxies_all` with the given ``dmax``.
    dists : array-like, optional
        Distances in Mpc for each galaxy in ``nbg_coords``. Required for
        color-coded plotting.
    glist : array-like of int, optional
        Indices selecting a subset of ``nbg_coords`` to plot. If None, all
        coordinates are plotted.
    radec : bool, optional
        If True, plots in RA/Dec coordinates. If False (default), plots in
        Galactic coordinates.
    dmax : astropy.units.Quantity, optional
        Maximum distance for loading galaxies when ``nbg_coords`` is None.
        Default is 10 Mpc.

    Returns
    -------
    sp : skyproj.Skyproj
        The updated map object (returned when plotting gray dots).
    sp, cbh : skyproj.Skyproj, matplotlib.cm.ScalarMappable
        The updated map and scalar mappable (returned when color-coding).

    Examples
    --------
    Add gray background galaxies to an existing map:

    >>> sp = add_extra_galaxies(sp)

    Add color-coded background galaxies sharing the scale of a prior plot:

    >>> sp, cbh = add_extra_galaxies(sp, cbh=cbh, nbg_coords=coords, dists=d)
    """
    if nbg_coords is None:
        nbg, nbg_coords = read_galaxies_all(dmax=dmax,kara=True)
        dists = nbg['D (Mpc)']

    if glist is None:
        nbg_to_plot = nbg_coords
        dists_to_plot = dists
    else:
        nbg_to_plot = nbg_coords[glist]
        dists_to_plot = dists[glist] if dists is not None else None

    if dists is None or cbh is None:
        if radec:
            sp.plot(nbg_to_plot.ra.value, nbg_to_plot.dec.value,
                    'o', mfc='gray', mec='none', size=2, zorder=2.3)
        else:
            sp.plot(nbg_to_plot.galactic.l.value, nbg_to_plot.galactic.b.value,
                    'o', mfc='gray', mec='none', size=2, zorder=2.3)
        return sp
    else:
        norm = cbh.norm
        cmap = cbh.get_cmap()
        if radec:
            cbh = sp.scatter(nbg_to_plot.ra.value, nbg_to_plot.dec.value,
                             c=dists_to_plot, edgecolors="none", cmap=cmap,
                             norm=norm, s=30, zorder=2.3)
        else:
            cbh = sp.scatter(nbg_to_plot.galactic.l.value,
                             nbg_to_plot.galactic.b.value,
                             c=dists_to_plot, edgecolors="none", cmap=cmap,
                             norm=norm, s=30, zorder=2.3)
        return sp, cbh


def plot_gcs(sp, gcs=None, gclist=None, annotate=False, cmap_name='plasma',
             norm=None, angles=None, radii=None, vb=False):
    """Plot Milky Way globular clusters on a skyproj map.

    Plots globular clusters colored by heliocentric distance on a logarithmic
    scale. Optionally annotates clusters with names and leader lines.

    Parameters
    ----------
    sp : skyproj.Skyproj
        Map object on which to plot the clusters.
    gcs : astropy.table.Table, optional
        Table of globular cluster data as returned by `read_gcs`. If None,
        loads via `read_gcs(vb=vb)`.
    gclist : array-like of int, optional
        Indices into ``gcs`` selecting a subset of clusters to plot. If None,
        all clusters are plotted.
    annotate : bool, optional
        If True, annotates each cluster with its name and a leader line.
        Default is False.
    cmap_name : str, optional
        Name of the matplotlib colormap to use for distance coloring.
        Default is ``'plasma'``.
    norm : matplotlib.colors.Normalize, optional
        Normalization for the colormap. If None, uses
        ``LogNorm(vmin=1.0, vmax=100.)``.
    angles : array-like of float, optional
        Angles in radians for annotation leader line directions, one per
        cluster. If None (and ``annotate=True``), random offsets are used.
    radii : array-like of float, optional
        Radii in points for annotation offsets. If None (and
        ``annotate=True``), a fixed radius of 30 points is used.
    vb : bool, optional
        Passed to `read_gcs` when ``gcs`` is None. If True, loads the
        Vasiliev & Baumgardt 2021 catalog instead of LVDB. Default is False.

    Returns
    -------
    cbh : matplotlib.cm.ScalarMappable
        Scalar mappable containing the colormap and normalization, suitable
        for passing to ``plt.colorbar()``.

    Examples
    --------
    Plot all globular clusters on a predefined skyproj:

    >>> cbh = plot_gcs(sp)

    Plot a subset with annotations using the V&B catalog:

    >>> cbh = plot_gcs(sp, gclist=[0, 5, 10], annotate=True, vb=True)
    """
    from matplotlib.colors import LogNorm
    from astropy.coordinates import SkyCoord

    if gcs is None:
        gcs = read_gcs(vb=vb)

    if gclist is not None:
        gcs = gcs[gclist]

    gcradec = get_sky_coords(gcs)

    cmap = plt.get_cmap(cmap_name)
    if norm is None:
        norm = LogNorm(vmin=1.0, vmax=100.)

    dists = gcs['distance'] * u.kpc

    cbh = sp.scatter(gcradec.galactic.l.value, gcradec.galactic.b.value,
                     c=dists, edgecolors="black", cmap=cmap, norm=norm, s=20)

    if annotate:
        for i, l in enumerate(gcs['name']):
            x, y = sp.proj(gcradec[i].galactic.l.deg, gcradec[i].galactic.b.deg)
            if angles is None or radii is None:
                r = 30
                costh = 2 * np.random.rand() - 1
            else:
                r = radii[i]
                costh = np.cos(angles[i])
            xt, yt = r * costh, r * np.sqrt(1 - costh * costh)
            sp.ax.annotate(l, (x, y), (xt, yt),
                           textcoords='offset points',
                           arrowprops=dict(arrowstyle='-', color='k'),
                           va='top', ha='left')

    return cbh


def plot_streams(sp, mws=None, dist=False, streamlist=None, annotate=False):
    """Plot stellar streams from the galstreams catalog on a skyproj map.

    Plots stream tracks in Galactic coordinates, either in alternating colors
    or color-coded by heliocentric distance. Streams without distance
    measurements are shown in gray when ``dist=True``.

    Parameters
    ----------
    sp : skyproj.Skyproj
        Map object on which to plot the streams.
    mws : galstreams.MWStreams, optional
        Pre-loaded stream catalog instance. If None, the catalog is loaded
        internally and deleted afterward to save memory.
    dist : bool, optional
        If True, colors stream tracks by heliocentric distance using a
        logarithmic ``'magma'`` colormap (5--60 kpc). If False (default),
        streams are plotted in alternating colors.
    streamlist : list of str, optional
        Names of streams to plot. If None, all streams in the catalog are
        plotted.
    annotate : bool, optional
        If True, labels each stream with its catalog ID number at one of its
        endpoints. Default is False.

    Returns
    -------
    cbh : matplotlib.cm.ScalarMappable or None
        Scalar mappable for the distance colormap when ``dist=True``, suitable
        for passing to ``plt.colorbar()``. Returns None if ``dist=False``.

    Examples
    --------
    Plot all streams with default alternating colors:

    >>> cbh = plot_streams(sp)

    Plot a subset of streams colored by distance from a pre-loaded catalog:

    >>> mws = galstreams.MWStreams(verbose=False, implement_Off=False,
    ...                            print_topcat_friendly_files=False)
    >>> s_pole = mws.get_track_names_in_sky_window(
    ...     [0, 359]*u.deg, [-90, -30]*u.deg, frame=astropy.coordinates.Galactic)
    >>> cbh = plot_streams(sp, mws=mws, dist=True, streamlist=s_pole, annotate=True)
    """
    import galstreams

    if mws is None:
        mws = galstreams.MWStreams(verbose=False, implement_Off=False,
                                   print_topcat_friendly_files=False)
        cleanup = True
    else:
        cleanup = False

    if streamlist is None:
        streamlist = mws.keys()

    cbh = None
    for st in streamlist:
        if dist:
            cmap = plt.get_cmap('magma')
            norm = matplotlib.colors.LogNorm(vmin=5.0, vmax=60.)
            if mws.summary.has_D[st] == 1:
                cbh = sp.scatter(
                    mws[st].track.galactic.l.wrap_at(180 * u.deg).value,
                    mws[st].track.galactic.b.value,
                    c=mws[st].track.distance, alpha=0.75, s=10, marker='.',
                    cmap=cmap, norm=norm,
                    label="{ID:.0f}={Name}".format(ID=mws[st].ID,
                                                   Name=mws[st].stream_name))
            else:
                sp.scatter(
                    mws[st].track.galactic.l.wrap_at(180 * u.deg).value,
                    mws[st].track.galactic.b.value,
                    c='gray', alpha=0.75, s=10, marker='.',
                    label="{ID:.0f}={Name}".format(ID=mws[st].ID,
                                                   Name=mws[st].stream_name))
        else:
            sp.scatter(
                mws[st].track.galactic.l.wrap_at(180 * u.deg).value,
                mws[st].track.galactic.b.value,
                alpha=0.75, s=10, marker='.',
                label="{ID:.0f}={Name}".format(ID=mws[st].ID,
                                               Name=mws[st].stream_name))

        if annotate:
            xo = mws[st].end_points.galactic.l.wrap_at(180 * u.deg)[0].value
            yo = mws[st].end_points.galactic.b[0].value
            ext = sp.get_extent()
            lmin, lmax = min(ext[:2]), max(ext[:2])
            bmin, bmax = min(ext[2:]), max(ext[2:])
            if xo < lmin or xo > lmax or yo < bmin or yo > bmax:
                xo = mws[st].end_points.galactic.l.wrap_at(180 * u.deg)[1].value
                yo = mws[st].end_points.galactic.b[1].value
            sp.ax.text(xo, yo, mws[st].ID)

    if cleanup:
        del mws

    return cbh
