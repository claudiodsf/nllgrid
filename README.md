# NLLGrid

Python class for reading and writing
[NonLinLoc](http://alomax.free.fr/nlloc) grid files.

(c) 2015-2021 Claudio Satriano, Natalia Poiata


## Installation

### Using pip and PyPI (preferred method)

The latest release of NLLGrid is available on the
[Python Package Index](https://pypi.org/project/nllgrid/).

You can install it easily through `pip`:

    pip install nllgrid

### From nllgrid GitHub releases

Download the latest release from the
[releases page](https://github.com/claudiodsf/nllgrid/releases),
in `zip` or `tar.gz` format, then:

    pip install nllgrid-X.Y.zip

or

    pip install nllgrid-X.Y.tar.gz

Where, `X.Y` is the version number (e.g., `1.3`).
You don't need to uncompress the release files yourself.

### From nllgrid GitHub repository

If you need a recent feature that is not in the latest release (see the
`unreleased` section in [CHANGELOG](CHANGELOG.md)), you want to use the source
code from the
[nllgrid GitHub repository](https://github.com/claudiodsf/nllgrid).

For that, clone the project:

    git clone https://github.com/claudiodsf/nllgrid.git

(avoid using the "Download ZIP" option from the green "Code" button, since
version number is lost), then install the code from within the `nllgrid`
main directory by running:

    pip install .

(use `pip install -e .` to install in developer mode).



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
