# NLLGrid

Python class for reading and writing
[NonLinLoc](http://alomax.free.fr/nlloc) grid files.

(c) 2015-2018 Claudio Satriano, Natalia Poiata

## Installation

Clone or download the project from GitHub, uncompress the archive
(if needed), and install the codes by running (from within the main
directory):

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
