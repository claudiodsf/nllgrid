# NLLGrid

Python class for reading and writing
[NonLinLoc](http://alomax.free.fr/nlloc) grid files.

(c) 2015-2024 Claudio Satriano, Natalia Poiata, Robert Pickle

## v1.5.2 - 2024-09-30

- Fix for missing test data in PyPI package

## v1.5.1 - 2024-09-30

- Fix for PyPI publishing

## v1.5 - 2024-09-30

- `NLLGrid.horizontal_recenter()` method to move the (x, y) origin of grid
  cartesian system to the grid center
- `NLLGrid.horizontal_rotate()` method to rotate the grid horizontally around
  its center
- `NLLGrid.nudge()` method to expand or contract in any direction a grid
- Added `TRANS_MERC` and `AZIMUTHAL_EQUIDIST` projections
- Fix projections when `map_rot != 0`

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
