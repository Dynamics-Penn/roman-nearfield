"""
plot_layouts.py
---------------
High-level functions for producing standard polar-projection sky figures
overlaying the Roman survey footprint with different object catalogs.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from .plot_catalogs import plot_gcs, plot_galaxies


def plot_polar_projection_gcs(fp, whichpole, gclist=None, annotate=False,
                               alph=1.0, angles=None, radii=None):
    """Plot a polar projection of the Roman footprint overlaid with globular clusters.

    Creates a Lambert equal-area polar projection centered on the north or
    south Galactic pole, draws the Roman HLWAS footprint, and overlays
    globular clusters colored by distance.

    Parameters
    ----------
    fp : numpy.ndarray
        HEALPix footprint map as returned by `get_footprint`, with integer
        tier values and ``hp.UNSEEN`` for uncovered pixels.
    whichpole : {'N', 'S'}
        Which Galactic pole to center the projection on. ``'N'`` shows
        latitudes 10 to 90 deg; ``'S'`` shows -90 to -10 deg.
    gclist : array-like of int, optional
        Indices into the globular cluster catalog selecting clusters to plot
        and annotate. If None, all clusters are plotted.
    annotate : bool, optional
        If True, annotates each cluster with its name. Default is False.
    alph : float, optional
        Alpha (transparency) for the footprint map. Default is 1.0.
    angles : array-like of float, optional
        Angles in radians for annotation leader line directions. Passed
        through to `plot_gcs`.
    radii : array-like of float, optional
        Radii in points for annotation offsets. Passed through to `plot_gcs`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The completed figure.

    Examples
    --------
    >>> fp = get_footprint()
    >>> fig = plot_polar_projection_gcs(fp, 'N')
    >>> fig.savefig('north_pole_gcs.png')
    """
    import skyproj
    from matplotlib.colors import ListedColormap

    col_dict = {1: "palegreen", 2: "limegreen"}
    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

    patch_HLWAS_w = mpatches.Patch(color=col_dict[1], label='HLWAS WIDE')
    patch_HLWAS_m = mpatches.Patch(color=col_dict[2], label='HLWAS MEDIUM')

    fig, ax = plt.subplots(figsize=(8, 8))
    if whichpole == 'S':
        sp = skyproj.LaeaSkyproj(ax=ax, lat_0=-90.0, galactic=True,
                                  extent=[0, 360, -90, -10])
        im, _, _, _ = sp.draw_hpxmap(fp, lon_range=(0, 360),
                                      lat_range=(-90, -10), alpha=alph,
                                      vmin=1, vmax=3, cmap=cm, zoom=False)
    else:
        sp = skyproj.LaeaSkyproj(ax=ax, lat_0=90.0, galactic=True,
                                  extent=[0, 360, 10, 90])
        im, _, _, _ = sp.draw_hpxmap(fp, lon_range=(0, 360),
                                      lat_range=(10, 90), alpha=alph,
                                      vmin=1, vmax=3, cmap=cm, zoom=False)

    cbh = plot_gcs(sp, gclist=gclist, annotate=annotate,
                   angles=angles, radii=radii)

    shandles = [patch_HLWAS_m, patch_HLWAS_w]
    handles, labels = sp.ax.get_legend_handles_labels()
    shandles.extend(handles)

    sp.legend(handles=shandles, loc='upper left', ncol=2, fontsize=10,
              bbox_to_anchor=(0.0, 0.0, 1.0, 0.05))

    sp.ax.tick_params(axis="x", labelsize=14)
    sp.ax.tick_params(axis="y", labelsize=14)

    cb = plt.colorbar(cbh, shrink=0.8)
    cb.set_label(label='Distance (kpc)', size=16)
    cb.ax.tick_params(labelsize=16)

    plt.tight_layout(h_pad=0.1, w_pad=0.1)

    return fig


def plot_polar_projection_nbgs(fp, gals, glist, whichpole, glabels=None,
                                alph=1.0):
    """Plot a polar projection of the Roman footprint overlaid with nearby galaxies.

    Creates a Lambert equal-area polar projection centered on the north or
    south Galactic pole, draws the Roman HLWAS footprint, and overlays nearby
    galaxies colored by distance in Mpc.

    Parameters
    ----------
    fp : numpy.ndarray
        HEALPix footprint map as returned by `get_footprint`.
    gals : pandas.DataFrame or astropy.table.Table
        Table of galaxy data passed through to `plot_galaxies`.
    glist : array-like of int
        Indices selecting the subset of galaxies to plot and annotate.
    whichpole : {'N', 'S'}
        Which Galactic pole to center the projection on.
    glabels : array-like of str, optional
        Labels for the annotated galaxies. If None, galaxies are numbered.
    alph : float, optional
        Alpha (transparency) for the footprint map. Default is 1.0.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The completed figure.

    Examples
    --------
    >>> fp = get_footprint()
    >>> fig = plot_polar_projection_nbgs(fp, nbg, glist_north, 'N',
    ...                                  glabels=names_north)
    """
    import skyproj
    from matplotlib.colors import ListedColormap

    col_dict = {1: "palegreen", 2: "limegreen"}
    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

    patch_HLWAS_w = mpatches.Patch(color=col_dict[1], label='HLWAS WIDE')
    patch_HLWAS_m = mpatches.Patch(color=col_dict[2], label='HLWAS MEDIUM')

    fig, ax = plt.subplots(figsize=(8, 8))
    if whichpole == 'S':
        sp = skyproj.LaeaSkyproj(ax=ax, lat_0=-90.0, galactic=True,
                                  extent=[0, 360, -90, -10])
        im, _, _, _ = sp.draw_hpxmap(fp, lon_range=(0, 360),
                                      lat_range=(-90, -10), alpha=alph,
                                      vmin=1, vmax=3, cmap=cm, zoom=False)
    else:
        sp = skyproj.LaeaSkyproj(ax=ax, lat_0=90.0, galactic=True,
                                  extent=[0, 360, 10, 90])
        im, _, _, _ = sp.draw_hpxmap(fp, lon_range=(0, 360),
                                      lat_range=(10, 90), alpha=alph,
                                      vmin=1, vmax=3, cmap=cm, zoom=False)

    cbh = plot_galaxies(sp, gals, glist=glist, annotate=True, glabels=glabels)

    shandles = [patch_HLWAS_m, patch_HLWAS_w]
    handles, labels = sp.ax.get_legend_handles_labels()
    shandles.extend(handles)

    sp.legend(handles=shandles, loc='upper left', ncol=2, fontsize=10,
              bbox_to_anchor=(0.0, 0.0, 1.0, 0.05))

    sp.ax.tick_params(axis="x", labelsize=14)
    sp.ax.tick_params(axis="y", labelsize=14)

    cb = plt.colorbar(cbh, shrink=0.8)
    cb.set_label(label='Distance (Mpc)', size=16)
    cb.ax.tick_params(labelsize=16)

    plt.tight_layout(h_pad=0.1, w_pad=0.1)

    return fig


def plot_polar_projection_nbgs_new(fp, gals, glist, whichpole, glabels=None,
                                    alph=1.0, annotate=True, scaled_gmap=None,
                                    ax=None, no_legend=False, no_cbar=False):
    """Plot a polar projection of the Roman footprint overlaid with nearby galaxies (extended version).

    Extended version of `plot_polar_projection_nbgs` that additionally
    supports an optional Gaia stellar density background, an externally
    supplied axes object, and flags to suppress the legend and colorbar.

    Parameters
    ----------
    fp : numpy.ndarray
        HEALPix footprint map as returned by `get_footprint`.
    gals : pandas.DataFrame or astropy.table.Table
        Table of galaxy data passed through to `plot_galaxies`.
    glist : array-like of int
        Indices selecting the subset of galaxies to plot and annotate.
    whichpole : {'N', 'S'}
        Which Galactic pole to center the projection on.
    glabels : array-like of str, optional
        Labels for the annotated galaxies. If None, galaxies are numbered.
    alph : float, optional
        Alpha (transparency) for the footprint map. Default is 1.0.
    annotate : bool, optional
        Whether to annotate galaxies with labels. Default is True.
    scaled_gmap : numpy.ndarray, optional
        Asinh-stretched Gaia stellar density map (as returned by
        `get_gaia_map`) to draw as a grayscale background. If None, no
        background is drawn.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw into. If None, a new 8x8 inch figure is created.
    no_legend : bool, optional
        If True, suppresses the survey tier legend. Default is False.
    no_cbar : bool, optional
        If True, suppresses the distance colorbar. Default is False.

    Returns
    -------
    sp : skyproj.Skyproj
        The skyproj map object, for further customization.
    cbh : matplotlib.cm.ScalarMappable
        Scalar mappable for the galaxy distance colormap.

    Examples
    --------
    >>> fp = get_footprint()
    >>> gmap = get_gaia_map()
    >>> sp, cbh = plot_polar_projection_nbgs_new(fp, nbg, glist_south, 'S',
    ...                                           scaled_gmap=gmap,
    ...                                           glabels=names_south)
    """
    import skyproj
    from matplotlib.colors import ListedColormap

    col_dict = {1: "palegreen", 2: "limegreen"}
    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

    patch_HLWAS_w = mpatches.Patch(color=col_dict[1], label='HLWAS WIDE')
    patch_HLWAS_m = mpatches.Patch(color=col_dict[2], label='HLWAS MEDIUM')

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    if whichpole == 'S':
        sp = skyproj.LaeaSkyproj(ax=ax, lat_0=-90.0, galactic=True,
                                  extent=[0, 360, -90, -10])
        lon_range = (0, 360)
        lat_range = (-90, -10)
    else:
        sp = skyproj.LaeaSkyproj(ax=ax, lat_0=90.0, galactic=True,
                                  extent=[0, 360, 10, 90])
        lon_range = (0, 360)
        lat_range = (10, 90)

    if scaled_gmap is not None:
        sp.draw_hpxmap(scaled_gmap, alpha=0.9, cmap='gray_r', vmin=0.5,
                       vmax=5.0, lon_range=lon_range, lat_range=lat_range)

    im, _, _, _ = sp.draw_hpxmap(fp, lon_range=lon_range, lat_range=lat_range,
                                   alpha=alph, vmin=1, vmax=3, cmap=cm,
                                   zoom=False)

    cbh = plot_galaxies(sp, gals, glist=glist, annotate=annotate,
                        glabels=glabels)

    shandles = [patch_HLWAS_m, patch_HLWAS_w]
    handles, labels = sp.ax.get_legend_handles_labels()
    shandles.extend(handles)

    if not no_legend:
        sp.legend(handles=shandles, loc='upper left', ncol=2, fontsize=10,
                  bbox_to_anchor=(0.0, 0.0, 1.0, 0.05))

    sp.ax.tick_params(axis="x", labelsize=14)
    sp.ax.tick_params(axis="y", labelsize=14)

    if not no_cbar:
        cb = plt.colorbar(cbh, shrink=0.8)
        cb.set_label(label='Distance (Mpc)', size=16)
        cb.ax.tick_params(labelsize=16)

    plt.tight_layout(h_pad=0.1, w_pad=0.1)

    return sp, cbh


def plot_polar_projection_sats(fp, gals, glist, whichpole, glabels=None,
                                alph=1.0, cbar=True, angles=None, radii=None):
    """Plot a polar projection of the Roman footprint overlaid with MW satellites.

    Similar to `plot_polar_projection_nbgs` but tailored to MW satellite
    galaxies: uses a distance colorbar in kpc rather than Mpc, and supports
    explicit annotation placement via ``angles`` and ``radii``.

    Parameters
    ----------
    fp : numpy.ndarray
        HEALPix footprint map as returned by `get_footprint`.
    gals : pandas.DataFrame or astropy.table.Table
        Table of satellite data passed through to `plot_galaxies`.
    glist : array-like of int
        Indices selecting the subset of satellites to plot and annotate.
    whichpole : {'N', 'S'}
        Which Galactic pole to center the projection on.
    glabels : array-like of str, optional
        Labels for the annotated satellites. If None, satellites are numbered.
    alph : float, optional
        Alpha (transparency) for the footprint map. Default is 1.0.
    cbar : bool, optional
        If True (default), adds a distance colorbar in kpc.
    angles : array-like of float, optional
        Angles in radians for annotation leader line directions, passed
        through to `plot_galaxies`.
    radii : array-like of float, optional
        Radii in points for annotation offsets, passed through to
        `plot_galaxies`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The completed figure.

    Examples
    --------
    >>> fp = get_footprint()
    >>> fig = plot_polar_projection_sats(fp, sats, sat_list_north, 'N',
    ...                                   glabels=sat_names_north)
    """
    import skyproj
    import matplotlib.ticker as mticker
    from matplotlib.colors import ListedColormap

    col_dict = {1: "palegreen", 2: "limegreen"}
    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])

    patch_HLWAS_w = mpatches.Patch(color=col_dict[1], label='HLWAS WIDE')
    patch_HLWAS_m = mpatches.Patch(color=col_dict[2], label='HLWAS MEDIUM')

    fig, ax = plt.subplots(figsize=(8, 9))
    if whichpole == 'S':
        sp = skyproj.LaeaSkyproj(ax=ax, lat_0=-90.0, galactic=True,
                                  extent=[0, 360, -90, -10])
        im, _, _, _ = sp.draw_hpxmap(fp, lon_range=(0, 360),
                                      lat_range=(-90, -10), alpha=alph,
                                      vmin=1, vmax=3, cmap=cm, zoom=False)
    else:
        sp = skyproj.LaeaSkyproj(ax=ax, lat_0=90.0, galactic=True,
                                  extent=[0, 360, 10, 90])
        im, _, _, _ = sp.draw_hpxmap(fp, lon_range=(0, 360),
                                      lat_range=(10, 90), alpha=alph,
                                      vmin=1, vmax=3, cmap=cm, zoom=False)

    cbh = plot_galaxies(sp, gals, glist=glist, annotate=True,
                        glabels=glabels, annotation_style='manual',
                        angles=angles, radii=radii)

    shandles = [patch_HLWAS_m, patch_HLWAS_w]
    handles, labels = sp.ax.get_legend_handles_labels()
    shandles.extend(handles)

    sp.legend(handles=shandles, loc='upper left', ncol=2, fontsize=10,
              bbox_to_anchor=(0.05, 0.05, 1.0, 0.05))

    sp.ax.tick_params(axis="x", labelsize=14, direction='in', pad=-15)
    sp.ax.tick_params(axis="y", labelsize=14, direction='in', pad=-22)

    if cbar:
        cb = plt.colorbar(cbh, fraction=0.046, pad=0.04,
                          ticks=[30, 50, 100, 200, 300],
                          format=mticker.FixedFormatter(
                              ['30', '50', '100', '200', '300']))
        cb.set_label(label='Distance [kpc]', size=16)
        cb.ax.tick_params(labelsize=16)
        fig.subplots_adjust(left=0.22)

    return fig
