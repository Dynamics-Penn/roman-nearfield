# roman-nearfield

This repository reproduces plots from "A field guide to the near field" (Sanderson et al in prep).
Footprints for the core community surveys were provided by Saurabh Jha and Javier Sanchez (HLWAS and GBTDS) and Bob Benjamin (GPS).

# UPDATED 11 March 2026!

This update refactors the functions used in the notebooks into a pip-installable package.
If you want the updated version in your local clone, you need a hard reset (stash your local changes first before you do it!):
```
git stash
git fetch origin
git reset --hard origin/main
git stash pop
```

Or, you can simply `pip install git+https://github.com/Dynamics-Penn/roman-nearfield.git` to get the new version and update your copies of the notebooks by hand.

## Reference

If you use this repo please cite [Sanderson et al in prep](https://www.overleaf.com/read/xgxwjmkmctnc#b4a1df) (will be on arxiv by March 9 2026). Thank you!

The package consolidates functions that were previously duplicated across
multiple Jupyter notebooks into a single importable library, covering:

- HEALPix coordinate transformations and survey footprint loading
- Reading nearby galaxy, globular cluster, and MW satellite catalogs
- Plotting catalogs on skyproj sky maps
- Selecting objects within survey footprints
- Proper motion / transverse velocity conversions

## Installation

You can now either install with `pip` or clone the repo, as described below.

### Cloning the repository

When you clone this repository you'll need to have Git's [Large File Storage](https://git-lfs.com/) extension installed.
Cloning the repo will download only pointers to many of the data files. To actually download them you will need to run

```bash
git lfs checkout
```

from within the cloned repo. Once you have run this command one time, you should not have to do anything differently than your normal git workflow unless you add a new large file with an extension not tracked by `git lfs`. See the [LFS wiki](https://github.com/git-lfs/git-lfs/wiki/Tutorial) for more about the extension.

### Installing the `roman_nearfield` package

The shared utility functions used across the notebooks are collected in the
`roman_nearfield` package. Install it in editable mode from the repo root:

```bash
pip install -e .
```

To also install the test dependencies:

```bash
pip install -e .[test]
```

You can also install directly from GitHub without cloning:

```bash
pip install git+https://github.com/Dynamics-Penn/roman-nearfield.git
```

## Usage

Import the package at the top of a notebook in place of the inline function
definitions:

```python
from roman_nearfield import (
    change_coord, get_footprint, read_gcs, plot_gcs,
    plot_galaxies, mu, vt
)
```

All public functions are also available directly from the top-level namespace:

```python
import roman_nearfield as rn

fp = rn.get_footprint()
gcs = rn.read_gcs()
```

## Data directory

The package reads data files from a configurable directory. By default this
is `../data/` relative to the current working directory, which matches the
layout used by the project notebooks. Override it by setting the environment
variable before importing:

```bash
export ROMAN_NEARFIELD_DATA_DIR=/path/to/your/data
```

Or at runtime before the first import:

```python
import os
os.environ["ROMAN_NEARFIELD_DATA_DIR"] = "/path/to/your/data"
import roman_nearfield
```

## Package structure

| Module | Contents |
|---|---|
| `healpix.py` | `change_coord`, `get_footprint`, `get_gaia_map`, `get_custom_colorbar` |
| `readers.py` | `read_galaxies`, `read_galaxies_all`, `read_gcs`, `read_satellites`, `get_sky_coords` |
| `plot_catalogs.py` | `plot_galaxies`, `plot_gcs`, `plot_streams`, `add_extra_galaxies`, `spread_clusters_kernel` |
| `plot_layouts.py` | `plot_polar_projection_gcs`, `plot_polar_projection_nbgs`, `plot_polar_projection_nbgs_new`, `plot_polar_projection_sats` |
| `footprint_select.py` | `select_gcs_in_footprint`, `select_gcs_in_gps_main`, `select_galaxies_in_footprint`, `select_satellites_in_footprint`, `get_gptds_fields` |
| `kinematics.py` | `mu`, `vt`, `velocity_limit`, `velocity_limit_year`, `get_year` |

## Running the tests

```bash
pytest test/ -v
```

## Dependencies

- numpy, scipy, matplotlib, pandas
- astropy, astroquery
- healpy, skyproj, galstreams

## License

This repo is released under the GNU Affero General Public License v3.0. See [LICENSE](LICENSE) for details.

## AI acknowledgment

Code in this repo was refactored with help from Claude Sonnet 4.6.
