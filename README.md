# NLLGrid

Python class for reading and writing [NonLinLoc] grid files.

(c) 2015-2025 Claudio Satriano, Natalia Poiata, Robert Pickle

[![changelog-badge]][changelog-link]
[![cf-badge]][cf-link]
[![PyPI-badge]][PyPI-link]
[![license-badge]][license-link]
[![docs-badge]][docs-link]

## Installation

### Using Anaconda

If you use [Anaconda], the latest release of nllgrid is available via
[conda-forge][cf-link].

To install, simply run:

    conda install -c conda-forge nllgrid

### Using pip and PyPI

The latest release of nllgrid is available on the
[Python Package Index][PyPI-link].

You can install it easily through `pip`:

    pip install nllgrid

### From nllgrid GitHub releases

Download the latest release from the
[releases page][releases-link],
in `zip` or `tar.gz` format, then:

    pip install nllgrid-X.Y.zip

or

    pip install nllgrid-X.Y.tar.gz

Where, `X.Y` is the version number (e.g., `1.3`).
You don't need to uncompress the release files yourself.

### Installing a development snapshot

If you need a recent feature that is not in the latest release (see the
`unreleased` section in [CHANGELOG][changelog-link]), you want to use the
more recent development snapshot from the
[nllgrid GitHub repository][github-repo].

#### Using pip

The easiest way to install the most recent development snapshot is to download
and install it through `pip`, using its builtin `git` client:

    pip install git+https://github.com/claudiodsf/nllgrid.git

Run this command again, from times to times, to keep NLLGrid updated with
the development version.

#### Cloning the NLLGrid GitHub repository

If you want to take a look at the source code (and possibly modify it 😉),
clone the project using `git`:

    git clone https://github.com/claudiodsf/nllgrid.git

or, using SSH:

    git clone git@github.com:claudiodsf/nllgrid.git

(avoid using the "Download ZIP" option from the green "Code" button, since
version number is lost).

Then, go into the `nllgrid` main directory and install the code in "editable
mode" by running:

    pip install -e .

You can keep your local NLLGrid repository updated by running `git pull`
from times to times. Thanks to `pip`'s "editable mode", you don't need to
reinstall NLLGrid after each update.

## Getting Started

### Reading a NLL grid

A NLL grid is composed of two files (`.hdr` and `.buf`).

To read a NLL grid, do:

```python
>>> from nllgrid import NLLGrid
>>> grd = NLLGrid('somegrid.hdr')
```

or, using the `.buf` filename:

```python
>>> grd = NLLGrid('somegrid.buf')
```

or even without any extension:

```python
>>> grd = NLLGrid('somegrid')
```

A grid description can be obtained by:

```python
>>> print(grd)
```

The grid data array is accessed by `grd.array`.
The grid can be plotted doing:

```python
>>> grd.plot()
```

Use Python introspection (e.g. `dir(grd)`) to see all the available
methods and attributes.


### Creating a NLL grid

Suppose that you have a 3D data array stored into a NumPy array
called `mydata`.

First, create an empty NLL grid object:

```python
>>> from nllgrid import NLLGrid
>>> grd = NLLGrid()
```

then, add the data array and information on grid sampling and grid
origin, e.g.:

```python
>>> grd.array = mydata
>>> grd.dx = 0.5  #km
>>> grd.dy = 0.5  #km
>>> grd.dz = 0.5  #km
>>> grd.x_orig = -10  #km
>>> grd.y_orig = -20. #km
>>> grd.z_orig = -1.  #km
```

Optionally, add a grid type and/or a geographic transformation:

```python
>>> grd.type = 'VELOCITY'
>>> grd.orig_lat = 40.63
>>> grd.orig_lon = 15.80
>>> grd.proj_name = 'LAMBERT'
>>> grd.first_std_paral = 38.
>>> grd.second_std_paral = 42.
>>> grd.proj_ellipsoid = 'WGS-84'
```

Finally, give a basename and write to disk:

```python
>>> grd.basename = 'mygrid'
>>> grd.write_hdr_file()
>>> grd.write_buf_file()
```

This will create the two files `mygrid.hdr` and `mygrid.buf`.

If you want to save your grid in double precision (required for
instance by NLDiffLoc), change `grd.float_type` to `'DOUBLE'` before
saving the grid (default is `'FLOAT'`):

```python
>>> grd.float_type = 'DOUBLE'
```

Note that if you want to use your grid as input for NonLinLoc
`Grid2Time` code, the grid type has to be `SLOW_LEN` and your grid
array has to be transformed into slowness (in s/km) and multiplied
by the grid step (in km).

[changelog-badge]: https://img.shields.io/badge/Changelog-136CB6.svg
[changelog-link]: CHANGELOG.md
[cf-badge]: http://img.shields.io/conda/vn/conda-forge/nllgrid.svg
[cf-link]: https://anaconda.org/conda-forge/nllgrid
[PyPI-badge]: http://img.shields.io/pypi/v/nllgrid.svg
[PyPI-link]: https://pypi.python.org/pypi/nllgrid
[license-badge]: https://img.shields.io/badge/license-GPLv3-green
[license-link]: https://www.gnu.org/licenses/gpl-3.0.html
[docs-badge]: https://readthedocs.org/projects/nllgrid/badge/?version=latest
[docs-link]: https://nllgrid.readthedocs.io/en/latest/?badge=latest

[NonLinLoc]: http://alomax.free.fr/nlloc
[Anaconda]: https://www.anaconda.com/products/individual
[releases-link]: https://github.com/claudiodsf/nllgrid/releases
[github-repo]: https://github.com/claudiodsf/nllgrid
