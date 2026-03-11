"""
roman_nearfield
===============
A utility package for Roman Near-Field Cosmology notebooks.

Provides shared functions for reading astronomical catalogs, plotting
objects on HEALPix sky maps, selecting objects within survey footprints,
and computing proper motion / velocity conversions.

The data directory is resolved at import time from the environment variable
``ROMAN_NEARFIELD_DATA_DIR``. If the variable is not set, it defaults to
``../data/`` relative to the current working directory.

To override, set the environment variable before importing::

    export ROMAN_NEARFIELD_DATA_DIR=/path/to/your/data

Or at runtime before the first import::

    import os
    os.environ["ROMAN_NEARFIELD_DATA_DIR"] = "/path/to/your/data"
    import roman_nearfield
"""

import os
from pathlib import Path

__version__ = "0.1.0"

DATA_DIR = Path(os.environ.get("ROMAN_NEARFIELD_DATA_DIR", "../data/"))

# -- healpix ------------------------------------------------------------------
from .healpix import (
    change_coord,
    get_footprint,
    get_gaia_map,
    get_custom_colorbar,
)

# -- readers ------------------------------------------------------------------
from .readers import (
    read_galaxies,
    read_galaxies_all,
    read_gcs,
    read_satellites,
    get_sky_coords,
)

# -- plot_catalogs ------------------------------------------------------------
from .plot_catalogs import (
    spread_clusters_kernel,
    plot_galaxies,
    add_extra_galaxies,
    plot_gcs,
    plot_streams,
)

# -- plot_layouts -------------------------------------------------------------
from .plot_layouts import (
    plot_polar_projection_gcs,
    plot_polar_projection_nbgs,
    plot_polar_projection_nbgs_new,
    plot_polar_projection_sats,
)

# -- footprint_select ---------------------------------------------------------
from .footprint_select import (
    select_gcs_in_footprint,
    select_gcs_in_gps_main,
    select_galaxies_in_footprint,
    select_satellites_in_footprint,
    get_gptds_fields,
)

# -- kinematics ---------------------------------------------------------------
from .kinematics import (
    mu,
    vt,
    velocity_limit,
    velocity_limit_year,
    get_year,
)

__all__ = [
    # config
    "DATA_DIR",
    "__version__",
    # healpix
    "change_coord",
    "get_footprint",
    "get_gaia_map",
    "get_custom_colorbar",
    # readers
    "read_galaxies",
    "read_galaxies_all",
    "read_gcs",
    "read_satellites",
    "get_sky_coords",
    # plot_catalogs
    "spread_clusters_kernel",
    "plot_galaxies",
    "add_extra_galaxies",
    "plot_gcs",
    "plot_streams",
    # plot_layouts
    "plot_polar_projection_gcs",
    "plot_polar_projection_nbgs",
    "plot_polar_projection_nbgs_new",
    "plot_polar_projection_sats",
    # footprint_select
    "select_gcs_in_footprint",
    "select_gcs_in_gps_main",
    "select_galaxies_in_footprint",
    "select_satellites_in_footprint",
    "get_gptds_fields",
    # kinematics
    "mu",
    "vt",
    "velocity_limit",
    "velocity_limit_year",
    "get_year",
]
