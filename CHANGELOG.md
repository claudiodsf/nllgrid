# NLLGrid

Python class for reading and writing
[NonLinLoc](http://alomax.free.fr/nlloc) grid files.

(c) 2015-2021 Claudio Satriano, Natalia Poiata


## v1.4.2 - 2021-10-19

- Still another fix for version string


## v1.4.1 - 2021-10-12

- Fix version string reported by `nllgrid.__version__`


## v1.4 - 2021-10-12

- `NLLGrid.get_value()`: support for 2D grids
- `NLLGrid.plot()`: support for 2D grids
- Fixed: `NLLGrid.get_value()` now correctly returns angles,
  instead of `None`, for `ANGLE` and `ANGLE2D` grids


## v1.3 - 2021-07-09

- `NLLGrid.iproject()` method for inverse projection


## v1.2.3 - 2021-05-20

- Move to `python-versioneer`


## v1.2.2 - 2021-03-05

- Fixes for PyPI packaging


## v1.2.1 - 2021-03-03

- Installation instructions via pip


## v1.2 - 2021-03-03

- Bug fixes


## v1.1 - 2018-07-13

- Python 3 compatibility
- `NLLGrid.project()` method
- Bug fixes


## v1.0 - 2018-06-15

- Initial release
