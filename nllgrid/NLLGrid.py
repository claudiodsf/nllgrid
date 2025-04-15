# -*- coding: utf8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Reading and writing of NonLinLoc grid files.

:copyright:
    2013-2025 Claudio Satriano <satriano@ipgp.fr>,
              Natalia Poiata <poiata@ipgp.fr>,
              Robert Pickle <rpickle@gmail.com>
:license:
    GNU General Public License v3.0 or later
    (https://www.gnu.org/licenses/gpl-3.0-standalone.html)
"""
import contextlib
import math
from ctypes import Union, c_float, c_ushort
from copy import deepcopy
from collections.abc import Iterable
import numpy as np
from scipy.ndimage import zoom, rotate
from pyproj import Proj


valid_grid_types = (
    'VELOCITY',
    'VELOCITY_METERS',
    'SLOWNESS',
    'VEL2',
    'SLOW2',
    'SLOW2_METERS',
    'SLOW_LEN',
    'STACK',
    'TIME',
    'TIME2D',
    'PROB_DENSITY',
    'MISFIT',
    'ANGLE',
    'ANGLE2D'
)

valid_float_types = {
    # NLL_type: numpy_type
    'FLOAT': 'float32',
    'DOUBLE': 'float64'
}

valid_projections = (
    'NONE',
    'SIMPLE',
    'LAMBERT',
    'TRANS_MERC',
    'AZIMUTHAL_EQUIDIST'
)

valid_ellipsoids = (
    'WGS-84',
    'GRS-80',
    'WGS-72',
    'Australian',
    'Krasovsky',
    'International',
    'Hayford-1909',
    'Clarke-1880',
    'Clarke-1866',
    'Airy',
    'Bessel',
    'Hayford-1830',
    'Sphere'
)

# Dictionary mapping ellipsoid names between NLL (keys) and proj (values)
ellipsoid_name_mapping = {
    'WGS-84': 'WGS84',
    'GRS-80': 'GRS80',
    'WGS-72': 'WGS72',
    'Australian': 'aust_SA',
    'Krasovsky': 'krass',
    'International': 'new_intl',
    'Hayford-1909': 'intl',
    'Clarke-1880': 'clrk80',
    'Clarke-1866': 'clrk66',
    'Airy': 'airy',
    'Bessel': 'bessel',
    # 'Hayford-1830':
    'Sphere': 'sphere'
    }


class TakeOffAngles(Union):
    """
    Union-style class for decoding take off angles.

    The take off angles are stored in the grid file as a 16-bit unsigned
    integer. The first 8 bits are the take off angle in degrees, the
    second 8 bits are the take off angle in minutes.

    The class has two attributes, fval and ival, which are the float
    and integer representations of the take off angle.
    """

    _fields_ = [('fval', c_float),
                ('ival', c_ushort*2)]


class NLLGrid(object):
    """
    Class for manipulating NonLinLoc grid files.

    It provides methods to read and write grid files,
    compute statistics and plot.

    Parameters
    ----------
    basename : str, optional
        The name of the grid file.
    nx : int, optional
        Number of grid points in x direction.
    ny : int, optional
        Number of grid points in y direction.
    nz : int, optional
        Number of grid points in z direction.
    x_orig : float, optional
        X coordinate of the lower left point in map view.
    y_orig : float, optional
        Y coordinate of the lower left point in map view.
    z_orig : float, optional
        Z coordinate of the shallowest point.
    dx : float, optional
        Spacing of grid points in x direction.
    dy : float, optional
        Spacing of grid points in y direction.
    dz : float, optional
        Spacing of grid points in z direction.
    """
    #: The name of the grid file
    basename = None
    #: Number of grid points in x direction
    nx = 1.
    #: Number of grid points in y direction
    ny = 1.
    #: Number of grid points in z direction
    nz = 1.
    #: X coordinate of the lower left point in map view
    x_orig = 0.
    #: Y coordinate of the lower left point in map view
    y_orig = 0.
    #: Z coordinate of the shallowest point
    z_orig = 0.
    #: Spacing of grid points in x direction
    dx = 1.
    #: Spacing of grid points in y direction
    dy = 1.
    #: Spacing of grid points in z direction
    dz = 1.
    #: The grid data as a 3D numpy array
    __array = None
    #: The type of the grid as a string
    __type = None
    #: Datatype for floating point numbers (FLOAT or DOUBLE)
    __float_type = 'FLOAT'
    __np_float_type = valid_float_types[__float_type]
    #: Longitude of the grid origin
    __orig_lon = None
    #: Latitude of the grid origin
    __orig_lat = None
    #: The name of the projection
    __proj_name = None
    #: Latitude of the first standard parallel for LAMBERT projection
    first_std_paral = None
    #: Latitude of the second standard parallel for LAMBERT projection
    second_std_paral = None
    #: Rotation angle of the grid in map view, counter-clockwise
    map_rot = 0.
    #: The name of the projection ellipsoid
    __proj_ellipsoid = None
    #: The projection function to perform direct and inverse projections
    __proj_function = None
    #: Name of the station
    station = None
    #: X coordinate of the station
    sta_x = None
    #: Y coordinate of the station
    sta_y = None
    #: Z coordinate of the station
    sta_z = None
    #: Coordinates of the grid mean point
    xyz_mean = None
    #: Covariance matrix of the grid
    xyz_cov = None
    #: The 68% confidence ellipsoid for the grid
    ellipsoid = None

    def __init__(self,
                 basename=None,
                 nx=1, ny=1, nz=1,
                 x_orig=0., y_orig=0., z_orig=0.,
                 dx=1., dy=1., dz=1.):
        """
        Initialize a NLLGrid object.

        If `basename` is not `None`, `read_hdr_file` and `read_buf_file`
        methods are called.
        """
        if basename is not None:
            self.basename = self.remove_extension(basename)
            self.read_hdr_file()
            self.read_buf_file()
        else:
            self.basename = None
            self.nx = nx
            self.ny = ny
            self.nz = nz
            self.x_orig = x_orig
            self.y_orig = y_orig
            self.z_orig = z_orig
            self.dx = dx
            self.dy = dy
            self.dz = dz

    def __str__(self):
        """Return a string representation of the object."""
        s = (
            f'basename: {self.basename}\n'
            f'nx: {self.nx} ny: {self.ny} nz: {self.nz}\n'
            f'x_orig: {self.x_orig} y_orig: {self.y_orig} '
            f'z_orig: {self.z_orig}\n'
            f'dx: {self.dx} dy: {self.dy} dz: {self.dz}\n'
            f'grid_type: {self.type}\n'
            f'float_type: {self.float_type}\n'
        )
        if self.station is not None:
            s += (
                f'station: {self.station} sta_x: {self.sta_x} '
                f'sta_y: {self.sta_y} sta_z: {self.sta_z}\n'
            )
        s += f'transform: {self.get_transform_line()}'
        return s

    def _repr_pretty_(self, p, _cycle):
        """Pretty print."""
        p.text(str(self))

    def __getitem__(self, key):
        """
        Return the value at the given index.

        Make the grid object behave like a numpy array.

        Parameters
        ----------
        key : tuple
            The index tuple.

        Returns
        -------
        value : float
            The value at the given index.
        """
        if self.type in ['ANGLE', 'ANGLE2D']:
            return self.dip[key]
        if self.array is not None:
            return self.array[key]

    @property
    def array(self):
        """The grid data as a 3D numpy array."""
        return self.__array

    @array.setter
    def array(self, array_data):
        array_data = np.asarray(array_data)
        if array_data.ndim != 3:
            raise ValueError('Only 3D arrays are supported')
        self.nx, self.ny, self.nz = array_data.shape
        self.__array = array_data

    @property
    def type(self):
        """The type of the grid as a string"""
        return self.__type

    @type.setter
    def type(self, grid_type):
        try:
            grid_type = grid_type.upper()
        except AttributeError as e:
            raise ValueError('Grid type must be a string') from e
        if grid_type not in valid_grid_types:
            msg = f'Invalid grid type: {grid_type}\n'
            msg += f'Valid grid types are: {valid_grid_types}'
            raise ValueError(msg)
        self.__type = grid_type

    @property
    def float_type(self):
        """The datatype for floating point numbers (FLOAT or DOUBLE)"""
        return self.__float_type

    @float_type.setter
    def float_type(self, float_type):
        try:
            float_type = float_type.upper()
        except AttributeError as e:
            raise ValueError('Float type must be a string') from e
        if float_type not in valid_float_types:
            msg = f'Invalid float type: {float_type}\n'
            msg += f'Valid float types are: {tuple(valid_float_types.keys())}'
            raise ValueError(msg)
        self.__np_float_type = valid_float_types[float_type]
        self.__float_type = float_type

    @property
    def proj_name(self):
        """The name of the projection"""
        return self.__proj_name

    @proj_name.setter
    def proj_name(self, pname):
        try:
            pname = pname.upper()
        except AttributeError as e:
            raise ValueError('Projection name must be a string') from e
        if pname not in valid_projections:
            msg = f'Invalid projection name: {pname}\n'
            msg += f'Valid projection names: {valid_projections}'
            raise ValueError(msg)
        self.__proj_name = pname
        # Reset proj_function
        self.__proj_function = None

    @property
    def proj_ellipsoid(self):
        """The name of the projection ellipsoid"""
        return self.__proj_ellipsoid

    @proj_ellipsoid.setter
    def proj_ellipsoid(self, ellipsoid):
        try:
            # here we just use upper() to check if ellipsoid is a string
            ellipsoid.upper()
        except AttributeError as e:
            raise ValueError('Ellipsoid must be a string') from e
        if ellipsoid not in valid_ellipsoids:
            msg = f'Invalid ellipsoid: {ellipsoid}\n'
            msg += f'Valid ellipsoids: {valid_ellipsoids}'
            raise ValueError(msg)
        self.__proj_ellipsoid = ellipsoid
        # Reset proj_function
        self.__proj_function = None

    @property
    def orig_lon(self):
        """The longitude of the grid origin"""
        return self.__orig_lon

    @orig_lon.setter
    def orig_lon(self, orig_lon):
        self.__orig_lon = float(orig_lon)
        # Reset proj_function
        self.__proj_function = None

    @property
    def orig_lat(self):
        """The latitude of the grid origin"""
        return self.__orig_lat

    @orig_lat.setter
    def orig_lat(self, orig_lat):
        self.__orig_lat = float(orig_lat)
        # Reset proj_function
        self.__proj_function = None

    def remove_extension(self, basename):
        """
        Remove '.hdr' or '.buf' suffixes, if present.

        Parameters
        ----------
        basename : str
            The basename of the grid file.

        Returns
        -------
        str
            The basename without the '.hdr' or '.buf' suffixes.

        Example
        -------
        >>> grd = NLLGrid()
        >>> grd.remove_extension('test.hdr')
        'test'
        >>> grd.remove_extension('test.buf')
        'test'
        >>> grd.remove_extension('test')
        'test'
        """
        bntmp = basename.rsplit('.hdr', 1)[0]
        return bntmp.rsplit('.buf', 1)[0]

    def init_array(self):
        """
        Initialize the grid array to zeros.

        Returns
        -------
        None

        Note
        ----
        This method sets the `array` attribute of the `NLLGrid` instance
        to a 3D numpy array of shape `(self.nx, self.ny, self.nz)` and data
        type `float`, filled with zeros.

        Example
        -------
        >>> grd = NLLGrid(nx=20, ny=20, nz=30, dx=2., dy=2., dz=2.)
        >>> grd.init_array()
        >>> grd.array[2, 4, 10] = 3.
        """
        self.array = np.zeros((self.nx, self.ny, self.nz), float)

    def read_hdr_file(self, basename=None):
        """
        Read header file of NLL grid format.

        The header file provides information about the 3D grid such as the
        number of gridpoints in each dimension, the origin, the cell spacing,
        the type of grid, the data type of the values, and the geographic
        projection.

        Parameters
        ----------
        basename : str, optional
            Basename of the header file or full file name.
            If provided, the `basename` attribute of the class instance
            will be updated.
            If not provided, the `basename` attribute of the class instance
            will be used.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If the header file is not found.

        Example
        -------
        >>> grd = NLLGrid()
        >>> grd.read_hdr_file('test.hdr')
        >>> print(grd)
        basename: test
        nx: 2 ny: 301 nz: 61
        x_orig: 0.0 y_orig: 0.0 z_orig: -3.0
        dx: 5.0 dy: 5.0 dz: 5.0
        grid_type: SLOW_LEN
        float_type: FLOAT
        transform: TRANSFORM  LAMBERT RefEllipsoid Clarke-1880
        LatOrig 15.000000  LongOrig -61.000000  FirstStdParal 10.000000
        SecondStdParal 20.000000  RotCW 0.000000
        """
        if basename is not None:
            self.basename = self.remove_extension(basename)
        filename = f'{self.basename}.hdr'

        # read header file
        with open(filename, 'r', encoding='utf8') as fp:
            lines = fp.readlines()

        # extract information
        vals = lines[0].split()
        self.nx = int(vals[0])
        self.ny = int(vals[1])
        self.nz = int(vals[2])
        self.x_orig = float(vals[3])
        self.y_orig = float(vals[4])
        self.z_orig = float(vals[5])
        self.dx = float(vals[6])
        self.dy = float(vals[7])
        self.dz = float(vals[8])
        self.type = vals[9]
        try:
            self.float_type = vals[10]
        except IndexError:
            self.float_type = 'FLOAT'

        lines.pop(0)

        for line in lines:
            vals = line.split()
            if not vals:
                # skip empty lines
                continue
            if vals[0] in ['TRANS', 'TRANSFORM']:
                if vals[1] == 'NONE':
                    self.proj_name = 'NONE'
                if vals[1] == 'SIMPLE':
                    self.proj_name = 'SIMPLE'
                    self.orig_lat = float(vals[3])
                    self.orig_lon = float(vals[5])
                    self.map_rot = float(vals[7])
                if vals[1] == 'LAMBERT':
                    self.proj_name = 'LAMBERT'
                    self.proj_ellipsoid = vals[3]
                    self.orig_lat = float(vals[5])
                    self.orig_lon = float(vals[7])
                    self.first_std_paral = float(vals[9])
                    self.second_std_paral = float(vals[11])
                    self.map_rot = float(vals[13])
                if vals[1] == 'TRANS_MERC':
                    self.proj_name = 'TRANS_MERC'
                    self.proj_ellipsoid = vals[3]
                    self.orig_lat = float(vals[5])
                    self.orig_lon = float(vals[7])
                    self.map_rot = float(vals[9])
                if vals[1] == 'AZIMUTHAL_EQUIDIST':
                    self.proj_name = 'AZIMUTHAL_EQUIDIST'
                    self.proj_ellipsoid = vals[3]
                    self.orig_lat = float(vals[5])
                    self.orig_lon = float(vals[7])
                    self.map_rot = float(vals[9])
            else:
                self.station = vals[0]
                self.sta_x = float(vals[1])
                self.sta_y = float(vals[2])
                self.sta_z = float(vals[3])

    def read_buf_file(self, basename=None):
        """
        Read buf file as a 3d array.

        The buffer file is a binary representation of the 3D array stored
        in the `array` attribute of the class instance.

        Parameters
        ----------
        basename : str, optional
            Basename of the buffer file or full file name.
            If provided, the `basename` attribute of the class instance
            will be updated.
            If not provided, the `basename` attribute of the class instance
            will be used.

        Raises
        ------
        FileNotFoundError
            If the buffer file is not found.
        ValueError
            If there are not enough data values in buf file.

        Example
        -------
        >>> grd = NLLGrid()
        >>> grd.read_buf_file('test.buf')
        >>> print(grd)
        basename: test
        nx: 1 ny: 1 nz: 1
        x_orig: 0.0 y_orig: 0.0 z_orig: 0.0
        dx: 1.0 dy: 1.0 dz: 1.0
        grid_type: None
        float_type: FLOAT
        transform: None
        """
        if basename is not None:
            self.basename = self.remove_extension(basename)
        filename = f'{self.basename}.buf'

        with open(filename, 'rb') as fp:
            nitems = self.nx * self.ny * self.nz
            buf = np.fromfile(fp, dtype=self.__np_float_type, count=nitems)
            if len(buf) < nitems:
                raise ValueError(
                    'Not enough data values in buf file! '
                    f'({len(buf)} < {nitems})')
        if self.type in ['ANGLE', 'ANGLE2D']:
            take_off_angles = (TakeOffAngles * len(buf))()
            for _i, _val in enumerate(buf):
                take_off_angles[_i].fval = _val
            self.azimuth = np.array(
                [t.ival[1]/10. for t in take_off_angles]
            ).reshape((self.nx, self.ny, self.nz))
            self.dip = np.array(
                [(t.ival[0]//16)/10. for t in take_off_angles]
            ).reshape((self.nx, self.ny, self.nz))
            self.quality = np.array(
                [t.ival[0] % 16 for t in take_off_angles]
            ).reshape((self.nx, self.ny, self.nz))
            self.azimuth[self.quality == 0] = np.nan
            self.dip[self.quality == 0] = np.nan
        else:
            self.array = np.array(buf).reshape((self.nx, self.ny, self.nz))

    def write_hdr_file(self, basename=None):
        """
        Write header file in NLL grid format.

        The header file provides information about the 3D grid such as the
        number of gridpoints in each dimension, the origin, the cell spacing,
        the type of grid, the data type of the values, and the geographic
        projection.

        Parameters
        ----------
        basename : str, optional
            Base name of the header file. If not provided, the `basename`
            attribute of the class instance will be used.

        Returns
        -------
        None
        """
        if basename is not None:
            self.basename = basename
        filename = f'{self.basename}.hdr'

        lines = [
            f'{self.nx} {self.ny} {self.nz}  '
            f'{self.x_orig:.6f} {self.y_orig:.6f} {self.z_orig:.6f}  '
            f'{self.dx:.6f} {self.dy:.6f} {self.dz:.6f} {self.type} '
            f'{self.float_type}\n'
        ]
        if self.station is not None:
            lines.append(
                f'{self.station} '
                f'{self.sta_x:.6f} {self.sta_y:.6f} {self.sta_z:.6f}\n'
            )
        line = self.get_transform_line()
        if line is not None:
            lines.append(f'{line}\n')

        with open(filename, 'w', encoding='utf8') as fp:
            for line in lines:
                fp.write(line)

    def write_buf_file(self, basename=None):
        """
        Write buffer file as a 3D array.

        The buffer file is a binary representation of the 3D array stored
        in the `array` attribute of the class instance.

        Parameters
        ----------
        basename : str, optional
            Base name of the buffer file. If not provided, the `basename`
            attribute of the class instance will be used.

        Raises
        ------
        NotImplementedError
            If the grid type is 'ANGLE' or 'ANGLE2D'.
            Writing buf file is not supported for these grid types.

        Returns
        -------
        None
        """
        if self.type in ['ANGLE', 'ANGLE2D']:
            raise NotImplementedError(
                f'Writing buf file not implemented for {self.type} grid.')
        if self.array is None:
            return
        if basename is not None:
            self.basename = basename
        filename = f'{self.basename}.buf'
        with open(filename, 'wb') as fp:
            self.array.astype(self.__np_float_type).tofile(fp)

    def get_transform_line(self):
        """
        Get the transform line in NLL hdr format.

        Returns
        -------
        line : str
            A string representing the transform line in NLL hdr format.
        """
        if self.proj_name == 'NONE':
            return 'TRANSFORM  NONE'
        if self.proj_name == 'SIMPLE':
            return (
                'TRANSFORM  SIMPLE  '
                f'LatOrig {self.orig_lat:.6f}  LongOrig {self.orig_lon:.6f}  '
                f'RotCW {self.map_rot:.6f}'
            )
        if self.proj_name == 'LAMBERT':
            return (
                f'TRANSFORM  LAMBERT RefEllipsoid {self.proj_ellipsoid}  '
                f'LatOrig {self.orig_lat:.6f}  LongOrig {self.orig_lon:.6f}  '
                f'FirstStdParal {self.first_std_paral:.6f}  '
                f'SecondStdParal {self.second_std_paral:.6f}  '
                f'RotCW {self.map_rot:.6f}'
            )
        if self.proj_name == 'TRANS_MERC':
            return (
                f'TRANSFORM  TRANS_MERC RefEllipsoid {self.proj_ellipsoid}  '
                f'LatOrig {self.orig_lat:.6f}  LongOrig {self.orig_lon:.6f}  '
                f'RotCW {self.map_rot:.6f}'
            )
        if self.proj_name == 'AZIMUTHAL_EQUIDIST':
            return (
                'TRANSFORM  AZIMUTHAL_EQUIDIST '
                f'RefEllipsoid {self.proj_ellipsoid}  '
                f'LatOrig {self.orig_lat:.6f}  LongOrig {self.orig_lon:.6f}  '
                f'RotCW {self.map_rot:.6f}'
            )

    def get_xyz(self, i, j, k):
        """
        Get cartesian coordinates (x, y, z) for grid indexes (i, j, k).

        Parameters
        ----------
        i : int
            The index along the x-axis.
        j : int
            The index along the y-axis.
        k : int
            The index along the z-axis.

        Returns
        -------
        x : float
            The x-coordinate in cartesian space.
        y : float
            The y-coordinate in cartesian space.
        z : float
            The z-coordinate in cartesian space.
        """
        x = i * self.dx + self.x_orig
        y = j * self.dy + self.y_orig
        z = k * self.dz + self.z_orig
        return x, y, z

    def get_ijk(self, x, y, z):
        """
        Get grid indexes (i, j, k) for cartesian coordinates (x, y, z).

        Parameters
        ----------
        x : float
            The x-coordinate in cartesian space.
        y : float
            The y-coordinate in cartesian space.
        z : float
            The z-coordinate in cartesian space.

        Returns
        -------
        i : int
            The index along the x-axis.
        j : int
            The index along the y-axis.
        k : int
            The index along the z-axis.
        """
        i = np.floor((x - self.x_orig) / self.dx).astype(int)
        j = np.floor((y - self.y_orig) / self.dy).astype(int)
        k = np.floor((z - self.z_orig) / self.dz).astype(int)
        return i, j, k

    def get_ijk_max(self):
        """
        Return the indexes (i, j, k) of the grid max point.

        Returns
        -------
        tuple of ints or None
            The 3D index of the grid max point.
            Returns None if `self.array` is None.
        """
        if self.array is None:
            return None
        return np.unravel_index(self.array.argmax(), self.array.shape)

    def get_ijk_min(self):
        """
        Return the indexes (i,j,k) of the grid min point.

        Returns
        -------
        tuple of ints or None
            The 3D index of the grid min point.
            Returns None if `self.array` is None.
        """
        if self.array is None:
            return None
        return np.unravel_index(self.array.argmin(), self.array.shape)

    def get_xyz_max(self):
        """
        Return the coordinates (x,y,z) of the grid max point.

        Returns
        -------
        tuple of float or None
            The 3D coordinates of the grid max point.
            Returns None if `self.array` is None.
        """
        ijk_max = self.get_ijk_max()
        return None if ijk_max is None else self.get_xyz(*ijk_max)

    def get_xyz_min(self):
        """
        Return the coordinates (x,y,z) of the grid min point.

        Returns
        -------
        tuple of float or None
            The 3D coordinates of the grid min point.
            Returns None if `self.array` is None.
        """
        ijk_min = self.get_ijk_min()
        return None if ijk_min is None else self.get_xyz(*ijk_min)

    def get_ijk_mean(self):
        """Return the indexes (i,j,k) of the grid mean point."""
        xyz_mean = self.get_xyz_mean()
        return None if xyz_mean is None else self.get_ijk(*xyz_mean)

    def get_xyz_mean(self):
        """
        Calculate and return the mean (x,y,z) coordinate of the grid.

        Returns
        -------
        xmean, ymean, zmean : float
            Mean x, y, and z coordinates, respectively.

        Note
        ----
        If the grid array is not set, None is returned.
        """
        if self.array is None:
            return None
        xx = np.arange(0, self.nx) * self.dx + self.x_orig
        yy = np.arange(0, self.ny) * self.dy + self.y_orig
        zz = np.arange(0, self.nz) * self.dz + self.z_orig
        yarray, xarray, zarray = np.meshgrid(yy, xx, zz)
        array_sum = self.array.sum()
        xmean = (xarray * self.array).sum()/array_sum
        ymean = (yarray * self.array).sum()/array_sum
        zmean = (zarray * self.array).sum()/array_sum
        self.xyz_mean = (xmean, ymean, zmean)
        return (xmean, ymean, zmean)

    def get_xyz_cov(self):
        """
        Return the covariance matrix of the grid with respect to the mean point
        in (x,y,z).

        Returns
        -------
        cov : numpy.ndarray, shape (3,3)
            The covariance matrix of the grid with respect to the mean point
            in (x,y,z). If the grid is None, returns None.
        """
        if self.array is None:
            return None
        xyz_mean = self.get_xyz_mean()
        xx = np.arange(0, self.nx) * self.dx + self.x_orig
        yy = np.arange(0, self.ny) * self.dy + self.y_orig
        zz = np.arange(0, self.nz) * self.dz + self.z_orig
        yarray, xarray, zarray = np.meshgrid(yy, xx, zz)
        array_sum = self.array.sum()
        cov = np.zeros((3, 3))
        cov[0, 0] = (np.power(xarray, 2) * self.array).sum()/array_sum \
            - (xyz_mean[0] * xyz_mean[0])
        cov[0, 1] = cov[1, 0] =\
            (xarray * yarray * self.array).sum()/array_sum \
            - (xyz_mean[0] * xyz_mean[1])
        cov[0, 2] = cov[2, 0] = \
            (xarray * zarray * self.array).sum()/array_sum \
            - (xyz_mean[0] * xyz_mean[2])
        cov[1, 1] = (np.power(yarray, 2) * self.array).sum()/array_sum \
            - (xyz_mean[1] * xyz_mean[1])
        cov[1, 2] = cov[2, 1] = \
            (yarray * zarray * self.array).sum()/array_sum \
            - (xyz_mean[1] * xyz_mean[2])
        cov[2, 2] = (np.power(zarray, 2) * self.array).sum()/array_sum \
            - (xyz_mean[2] * xyz_mean[2])
        self.xyz_cov = cov
        return cov

    def get_xyz_ellipsoid(self):
        """
        Return the 68% confidence ellipsoid.

        Calculates the 68% confidence ellipsoid from the covariance matrix
        obtained by :func:`get_xyz_cov`. The calculated ellipsoid object
        is stored as an instance attribute `ellipsoid`.

        Returns
        -------
        Ellipsoid3D
            The 68% confidence ellipsoid.

        Note
        ----
        This method is a python translation of the CalcErrorEllipsoid()
        function from the NonLinLoc package, written by Anthony Lomax.
        """
        # pylint: disable=import-outside-toplevel
        try:
            from .ellipsoid import Ellipsoid3D
        except ImportError:
            from ellipsoid import Ellipsoid3D
        # The following code is a python translation of the
        # CalcErrorEllipsoid() c-function from the NonLinLoc package,
        # written by Anthony Lomax
        cov = self.get_xyz_cov()
        if cov is None:
            return None

        u, s, _v = np.linalg.svd(cov)

        del_chi_2 = 3.53  # 3.53: value for 68% conf
        ell = Ellipsoid3D()
        ell.az1 = math.degrees(math.atan2(u[0, 0], u[1, 0]))
        if ell.az1 < 0.0:
            ell.az1 += 360.0
        ell.dip1 = math.degrees(math.asin(u[2, 0]))
        ell.len1 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[0])
        ell.az2 = math.degrees(math.atan2(u[0, 1], u[1, 1]))
        if ell.az2 < 0.0:
            ell.az2 += 360.0
        ell.dip2 = math.degrees(math.asin(u[2, 1]))
        ell.len2 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[1])
        ell.len3 = math.sqrt(del_chi_2) / math.sqrt(1.0 / s[2])

        self.ellipsoid = ell
        return ell

    def get_value(self, x, y, z, array=None):
        """
        Get the grid value at specified cartesian coordinates (x, y, z).

        Parameters
        ----------
        x : float
            The x coordinate.
        y : float
            The y coordinate.
        z : float
            The z coordinate.
        array : numpy.ndarray, optional
            The 3D array to use, by default None.
            If not provided, the instance's `array` attribute is used.

        Returns
        -------
        value : float or tuple of 3 float values
            The grid value at the specified cartesian coordinates.
            If the grid type is 'ANGLE' or 'ANGLE2D', a tuple of
            (azimuth, dip, quality) is returned.

        Raises
        ------
        NotImplementedError
            If the grid type is 'ANGLE' or 'ANGLE2D' and an array argument
            is provided.
        ValueError
            If the specified coordinates are outside of the grid's extent.
        """
        if array is None:
            array = self.array
        elif self.type in ['ANGLE', 'ANGLE2D']:
            raise NotImplementedError(
                f'"array" argument not implemented for {self.type} grid.')
        # Special case of 2D grids: y is epicentral distance
        # note: this doesn't work for GLOBAL grids
        if self.nx <= 2:
            # note: sta_x and sta_y are set to 0 if not defined in grid header
            # (e.g., for model grids)
            sta_x = self.sta_x if self.sta_x is not None else 0.
            sta_y = self.sta_y if self.sta_y is not None else 0.
            y = np.sqrt((x-sta_x)**2 + (y-sta_y)**2)
            x = self.x_orig
        min_x, max_x, min_y, max_y, min_z, max_z = self.get_extent()
        if not (min_x <= x <= max_x and min_y <= y <= max_y and
                min_z <= z <= max_z):
            raise ValueError(f'point {(x, y, z)} outside the grid.')
        i, j, k = self.get_ijk(x, y, z)
        if self.type in ['ANGLE', 'ANGLE2D']:
            azimuth = self.azimuth[i, j, k]
            dip = self.dip[i, j, k]
            quality = self.quality[i, j, k]
            return azimuth, dip, quality
        return array[i, j, k]

    def get_extent(self):
        """
        Get the grid extent in cartesian units.

        Returns
        -------
        extent : tuple
            Tuple of x_min, x_max, y_min, y_max, z_min, z_max values
            in cartesian units (generally km).
        """
        return (
            self.x_orig - self.dx / 2,
            self.x_orig + self.nx * self.dx + self.dx / 2,
            self.y_orig - self.dy / 2,
            self.y_orig + self.ny * self.dy + self.dy / 2,
            self.z_orig - self.dz / 2,
            self.z_orig + self.nz * self.dz + self.dz / 2,
        )

    def get_xy_extent(self):
        """
        Get the grid xy extent in cartesian units.

        Returns
        -------
        extent : tuple
            Tuple of x_min, x_max, y_min, y_max values in cartesian units
            (generally km).
        """
        return self.get_extent()[:4]

    def get_xz_extent(self):
        """
        Get the grid xz extent in cartesian units.

        Returns
        -------
        extent : tuple
            Tuple of x_min, x_max, z_min, z_max values in cartesian units
            (generally km).
        """
        return self.get_extent()[:2] + self.get_extent()[4:]

    def get_zx_extent(self):
        """
        Get the grid zx extent in cartesian units.

        Returns
        -------
        extent : tuple
            Tuple of z_min, z_max, x_min, x_max values in cartesian units
            (generally km).
        """
        return self.get_extent()[4:] + self.get_extent()[:2]

    def get_yz_extent(self):
        """
        Get the grid yz extent in cartesian units.

        Returns
        -------
        extent : tuple
            Tuple of y_min, y_max, z_min, z_max values in cartesian units
            (generally km).
        """
        return self.get_extent()[2:]

    def get_zy_extent(self):
        """
        Get the grid zy extent in cartesian units.

        Returns
        -------
        extent : tuple
            Tuple of z_min, z_max, y_min, y_max values in cartesian units
            (generally km).
        """
        return self.get_extent()[4:] + self.get_extent()[2:4]

    def max(self):
        """
        Get the maximum value of the grid.

        Returns
        -------
        float
            The maximum value of the grid.

        Note
        ----
        If the grid type is 'ANGLE' or 'ANGLE2D', the maximum value of
        `self.dip` is returned. Otherwise, the maximum value of `self.array`
        """
        if self.type in ['ANGLE', 'ANGLE2D']:
            return np.nanmax(self.dip)
        if self.array is not None:
            return np.nanmax(self.array)

    def resample(self, dx, dy, dz):
        """
        Resample the grid to the specified resolution.

        Parameters
        ----------
        dx : float
            The new x-resolution of the grid.
        dy : float
            The new y-resolution of the grid.
        dz : float
            The new z-resolution of the grid.

        Raises
        ------
        NotImplementedError
            If the grid type is 'ANGLE' or 'ANGLE2D', as resampling is not
            implemented for these grid types.
        """
        if self.type in ['ANGLE', 'ANGLE2D']:
            raise NotImplementedError(
                f'Resample not implemented for {self.type} grid.')
        zoom_x = self.dx / dx
        zoom_y = self.dy / dy
        zoom_z = self.dz / dz
        self.array = zoom(self.array, (zoom_x, zoom_y, zoom_z))
        self.nx, self.ny, self.nz = self.array.shape
        if self.type == 'SLOW_LEN':
            self.array *= dx / self.dx
        self.dx = dx
        self.dy = dy
        self.dz = dz

    def get_plot_axes(self, figure=None, ax_xy=None):
        """
        Get the axes for the three projections and colorbar axis.

        Parameters
        ----------
        figure : object, optional
            Matplotlib figure object. The default is None.
        ax_xy : object, optional
            Matplotlib axis object for x-y projection. The default is None.

        Returns
        -------
        ax_xy : object
            Matplotlib axis object for x-y projection.
        ax_xz : object
            Matplotlib axis object for x-z projection.
        ax_zy : object
            Matplotlib axis object for z-y projection.
        ax_cb : object
            Matplotlib axis object for colorbar.

        Note
        ----
        Requires matplotlib.
        If `ax_xy` is not provided, a new figure and axis will be created.
        If `figure` is not provided, the figure will be obtained from
        `ax_xy`.
        """
        # pylint: disable=import-outside-toplevel
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        xmin, xmax, ymin, ymax, zmin, zmax = self.get_extent()
        if figure is None and ax_xy is None:
            figure = plt.figure()
        if figure is None and ax_xy is not None:
            figure = ax_xy.get_figure()
        if ax_xy is None:
            ax_xy = figure.add_subplot(111)

        # Special case of 2D grids:
        # x-axis is epicentral distance (y values on grid)
        # y-axis is depth (z values on grid)
        if self.nx <= 2:
            divider = make_axes_locatable(ax_xy)
            ax_xy.set_xlim(ymin, ymax)
            ax_xy.set_ylim(zmax, zmin)
            ax_xy.set_aspect('equal', 'datalim')
            # color-bar
            ax_cb = divider.append_axes('bottom', size=0.2, pad=0.4)
            return ax_xy, ax_cb

        ratio = float(xmax - xmin) / (ymax - ymin)
        plot_xz_size = ((zmax - zmin)/(xmax - xmin))*100
        plot_yz_size = plot_xz_size / ratio
        plot_cbar_size = 5  # percent
        xz_size = f'{plot_xz_size} %'
        yz_size = f'{plot_yz_size} %'
        cb_size = f'{plot_cbar_size} %'

        # ax_xy
        divider = make_axes_locatable(ax_xy)
        plt.setp(ax_xy.get_xticklabels(), visible=False)
        ax_xy.set_xlim(xmin, xmax)
        ax_xy.set_ylim(ymin, ymax)
        ax_xy.set_aspect('equal', 'datalim')
        plt.setp(ax_xy.get_yticklabels(), rotation=0)

        # ax_yz
        ax_yz = divider.append_axes(
            'right', size=yz_size, pad=0.05, sharey=ax_xy)
        plt.setp(ax_yz.get_yticklabels(), visible=False)
        ax_yz.set_xlim(zmin, zmax)
        ax_yz.set_ylim(ymin, ymax)
        plt.setp(ax_yz.get_xticklabels(), rotation=90)
        plt.setp(ax_yz.get_yticklabels(), rotation=90)

        # ax_xz
        ax_xz = divider.append_axes(
            'bottom', size=xz_size, pad=0.05, sharex=ax_xy)
        ax_xz.set_xlim(xmin, xmax)
        ax_xz.set_ylim(zmax, zmin)

        # color-bar
        ax_cb = divider.append_axes('bottom', size=cb_size, pad=0.5)

        return ax_xy, ax_xz, ax_yz, ax_cb

    def plot(self, slice_index=None, handle=False, figure=None, ax_xy=None,
             vmin=None, vmax=None, cmap=None, line_color='white', array=None):
        """
        Plot the grid using three orthogonal projections.

        Parameters
        ----------
        slice_index : int or str, optional
            Index of the slice to plot. Use 'max' or 'min' to plot the slice
            at the grid's maximum or minimum value, respectively.
            Leave it to None to use the grid's middle slice.
        handle : bool, optional
            Whether to return the handle of the plot. The default is False.
        figure : object, optional
            Matplotlib figure object. The default is None.
        ax_xy : object, optional
            Matplotlib axis object for x-y projection. The default is None.
        vmin : float, optional
            Lower limit for the color scale.
            Leave it to None to use the minimum value of the grid.
        vmax : float, optional
            Upper limit for the color scale.
            Leave it to None to use the maximum value of the grid.
        cmap : object, optional
            Colormap to use for the plot.
            Leave it to None to use the default Matplotlib colormap.
        line_color : str, optional
            Color of the grid lines. The default is 'white'.
        array : array_like, optional
            Array to plot.
            Leave it to None to use the grid's `array` attribute.

        Returns
        -------
        fig : object
            Matplotlib figure object.

        Note
        ----
        Requires matplotlib.
        If `ax_xy` is not provided, a new figure and axis will be created.
        If `figure` is not provided, the figure will be obtained from
        `ax_xy`.
        """
        # pylint: disable=import-outside-toplevel
        import matplotlib.pyplot as plt
        from matplotlib import ticker

        if array is None:
            if self.array is None:
                return
            else:
                array = self.array

        # Special case of 2D grids:
        # x-axis is epicentral distance (y values on grid)
        # y-axis is depth (z values on grid)
        if self.nx <= 2:
            ax_xy, ax_cb = self.get_plot_axes(figure, ax_xy)
            if figure is None:
                figure = ax_xy.get_figure()
            hnd = ax_xy.imshow(np.transpose(array[0, :, :]),
                               vmin=vmin, vmax=vmax, cmap=cmap,
                               origin='lower', extent=self.get_yz_extent(),
                               zorder=-10)
            fmt = '%.1e' if np.nanmax(array) <= 0.01 else '%.2f'
            cb = figure.colorbar(
                hnd, cax=ax_cb, orientation='horizontal', format=fmt)
            cb.locator = ticker.LinearLocator(numticks=3)
            cb.update_ticks()
            if handle:
                return ax_xy, cb
            plt.show()
            return

        ax_xy, ax_xz, ax_yz, ax_cb = self.get_plot_axes(figure, ax_xy)
        if figure is None:
            figure = ax_xy.get_figure()

        if slice_index is None:
            slice_index = list(map(int, (self.nx/2, self.ny/2, self.nz/2)))
        if slice_index == 'max':
            slice_index = self.get_ijk_max()
        if slice_index == 'min':
            slice_index = self.get_ijk_min()

        if vmin is None:
            vmin = np.nanmin(array)
        if vmax is None:
            vmax = np.nanmax(array)

        hnd = ax_xy.imshow(np.transpose(array[:, :, slice_index[2]]),
                           vmin=vmin, vmax=vmax, cmap=cmap,
                           origin='lower', extent=self.get_xy_extent(),
                           zorder=-10)
        ax_xz.imshow(np.transpose(array[:, slice_index[1], :]),
                     vmin=vmin, vmax=vmax, cmap=cmap,
                     origin='lower', extent=self.get_xz_extent(),
                     aspect='auto', zorder=-10)
        ax_yz.imshow(array[slice_index[0], :, :],
                     vmin=vmin, vmax=vmax, cmap=cmap,
                     origin='lower', extent=self.get_zy_extent(),
                     aspect='auto', zorder=-10)

        x_slice, y_slice, z_slice = self.get_xyz(*slice_index)
        ax_xy.axhline(y_slice, color=line_color, linestyle='dashed', zorder=-1)
        ax_xy.axvline(x_slice, color=line_color, linestyle='dashed', zorder=-1)
        ax_xz.axhline(z_slice, color=line_color, linestyle='dashed', zorder=-1)
        ax_yz.axvline(z_slice, color=line_color, linestyle='dashed', zorder=-1)

        fmt = '%.1e' if np.nanmax(array) <= 0.01 else '%.2f'
        cb = figure.colorbar(
            hnd, cax=ax_cb, orientation='horizontal', format=fmt)
        cb.locator = ticker.LinearLocator(numticks=3)
        cb.update_ticks()

        if handle:
            return (ax_xy, ax_xz, ax_yz), cb
        else:
            plt.show()

    def plot_3D_point(self, axes, point, color='red'):
        """
        Plot a 3D point on the grid in three different projections.

        Parameters
        ----------
        axes : tuple of matplotlib.axes.Axes
            Tuple of 3 axes objects to plot the point.
        point : tuple of float
            Tuple of the 3 point grid coordinates.
        color : str, optional
            Color of the point. Default is 'red'.

        Raises
        ------
        NotImplementedError
            If the grid is not 3D.
        """
        if self.nx <= 2:
            raise NotImplementedError(
                'This method is supported only for 3D grids')
        ax_xy, ax_xz, ax_yz = axes
        ax_xy.scatter(point[0], point[1], color=color)
        ax_xz.scatter(point[0], point[2], color=color)
        ax_yz.scatter(point[2], point[1], color=color)

    def plot_ellipsoid(self, axes, ellipsoid=None, mean_xyz=None):
        """
        Plot an ellipsoid on the grid.

        Parameters
        ----------
        axes : tuple of matplotlib.axes.Axes
            Tuple of 3 axes objects to plot the ellipsoid.
        ellipsoid : object, optional
            Ellipsoid to plot. Default is `None`, in which case
            `self.get_xyz_ellipsoid()` is called.
        mean_xyz : tuple of floats, optional
            Mean of the ellipsoid. Default is `None`, in which case
            `self.get_xyz_mean()` is called.

        Raises
        ------
        NotImplementedError
            If the grid is not 3D.

        Note
        ----
        This method is supported only for 3D grids.
        The method uses the `Vect3D`, `ellipsiod2Axes`, and `toEllipsoid3D`
        functions from the `ellipsoid` module.
        """
        if self.nx <= 2:
            raise NotImplementedError(
                'This method is supported only for 3D grids')
        # pylint: disable=import-outside-toplevel
        try:
            from .ellipsoid import Vect3D, ellipsiod2Axes, toEllipsoid3D
        except ImportError:
            from ellipsoid import Vect3D, ellipsiod2Axes, toEllipsoid3D
        ax_xy, ax_xz, ax_yz = axes

        if ellipsoid is None:
            ellipsoid = self.get_xyz_ellipsoid()
        expect = Vect3D()
        if mean_xyz is None:
            mean_xyz = self.get_xyz_mean()
        expect.x, expect.y, expect.z = mean_xyz

        pax1, pax2, pax3 = ellipsiod2Axes(ellipsoid)

        ellArray12 = toEllipsoid3D(pax1, pax2, expect, 100)
        ellArray13 = toEllipsoid3D(pax1, pax3, expect, 100)
        ellArray23 = toEllipsoid3D(pax2, pax3, expect, 100)
        ell12 = np.array([(vect.x, vect.y) for vect in ellArray12])
        ell13 = np.array([(vect.x, vect.y) for vect in ellArray13])
        ell23 = np.array([(vect.x, vect.y) for vect in ellArray23])

        ax_xy.plot(ell12[:, 0], ell12[:, 1])
        ax_xy.plot(ell13[:, 0], ell13[:, 1])
        ax_xy.plot(ell23[:, 0], ell23[:, 1])

        ell12 = np.array([(vect.x, vect.z) for vect in ellArray12])
        ell13 = np.array([(vect.x, vect.z) for vect in ellArray13])
        ell23 = np.array([(vect.x, vect.z) for vect in ellArray23])
        ax_xz.plot(ell12[:, 0], ell12[:, 1])
        ax_xz.plot(ell13[:, 0], ell13[:, 1])
        ax_xz.plot(ell23[:, 0], ell23[:, 1])

        ell12 = np.array([(vect.y, vect.z) for vect in ellArray12])
        ell13 = np.array([(vect.y, vect.z) for vect in ellArray13])
        ell23 = np.array([(vect.y, vect.z) for vect in ellArray23])
        ax_yz.plot(ell12[:, 1], ell12[:, 0])
        ax_yz.plot(ell13[:, 1], ell13[:, 0])
        ax_yz.plot(ell23[:, 1], ell23[:, 0])

    @property
    def proj_function(self):
        """The projection function to perform direct and inverse projections"""
        if self.__proj_function is not None:
            return self.__proj_function
        if self.proj_name is None:
            raise RuntimeError('No geographical projection defined.')
        ellps = None  # Placeholder to silence linter warnings
        if self.proj_name != 'SIMPLE':
            try:
                ellps = ellipsoid_name_mapping[self.proj_ellipsoid]
            except KeyError as e:
                raise ValueError(
                    f'Ellipsoid not supported: {self.proj_ellipsoid}'
                ) from e
        if self.proj_name == 'LAMBERT':
            self.__proj_function = Proj(
                proj='lcc', lat_0=self.orig_lat, lon_0=self.orig_lon,
                lat_1=self.first_std_paral, lat_2=self.second_std_paral,
                ellps=ellps)
        elif self.proj_name == 'TRANS_MERC':
            self.__proj_function = Proj(
                proj='tmerc', lat_0=self.orig_lat, lon_0=self.orig_lon,
                ellps=ellps)
        elif self.proj_name == 'AZIMUTHAL_EQUIDIST':
            self.__proj_function = Proj(
                proj='aeqd', lat_0=self.orig_lat, lon_0=self.orig_lon,
                ellps=ellps)
        elif self.proj_name == 'SIMPLE':
            self.__proj_function = Proj(
                proj='eqc', lat_0=self.orig_lat, lon_0=self.orig_lon)
        else:
            raise ValueError(f'Projection not supported: {self.proj_name}')
        return self.__proj_function

    def project(self, lon, lat):
        """
        Project longitude and latitude coordinates into grid coordinates.

        Parameters
        ----------
        lon : float or array-like
            The longitude coordinates to be projected.
        lat : float or array-like
            The latitude coordinates to be projected.

        Returns
        -------
        float or array-like
            The projected grid x-coordinates.
        float or array-like
            The projected grid y-coordinates.
        """
        x, y = self.proj_function(lon, lat)
        x = np.array(x)
        y = np.array(y)
        x[np.isnan(lon)] = np.nan
        y[np.isnan(lat)] = np.nan
        x = x / 1000.
        y = y / 1000.
        # Rotate from East-North to map angle
        theta = np.radians(self.map_rot)
        x1 = x*np.cos(theta) - y*np.sin(theta)
        y1 = x*np.sin(theta) + y*np.cos(theta)
        x = x1
        y = y1
        # Try to return the same type of lon, lat
        if not isinstance(lon, np.ndarray):
            x = type(lon)(x) if isinstance(lon, Iterable) else float(x)
        if not isinstance(lat, np.ndarray):
            y = type(lat)(y) if isinstance(lat, Iterable) else float(y)
        return x, y

    def iproject(self, x, y):
        """
        Convert grid coordinates to longitude and latitude.

        Parameters
        ----------
        x : float or array_like
            x-coordinate(s) in the grid's cartesian coordinate system.
        y : float or array_like
            y-coordinate(s) in the grid's cartesian coordinate system.

        Returns
        -------
        float or numpy.ndarray
            Longitude(s) corresponding to `x` and `y`.
        float or numpy.ndarray
            Latitude(s) corresponding to `x` and `y`.
        """
        x = np.array(x)
        y = np.array(y)
        x = x * 1000.
        y = y * 1000.
        # Rotate back coordinates to East-North
        theta = np.radians(-self.map_rot)
        x1 = x*np.cos(theta) - y*np.sin(theta)
        y1 = x*np.sin(theta) + y*np.cos(theta)
        lon, lat = self.proj_function(x1, y1, inverse=True)
        return lon, lat

    def horizontal_recenter(self):
        """
        Move the origin of the grid's cartesian coordinate system to the grid
        center.

        This operation updates the values of `x_orig` and `y_orig`.

        Returns
        -------
        None

        Note
        ----
        If a geographical projection is available, the geographical coordinates
        of the grid center (the new (0, 0) point) are redefined.
        The absolute position of the grid in space does not change, but the
        grid projection is updated (`orig_lon` and `orig_lat`).
        Vertical coordinates are not modified.
        """
        xlen_half = 0.5 * self.nx * self.dx
        ylen_half = 0.5 * self.ny * self.dy
        # The following code block is to redefine the geographical
        # coordinates of the grid center (the new (0, 0) point).
        # Only executed if a geographical projection is available.
        with contextlib.suppress(RuntimeError):
            # Find coordinates of the grid center before recentering the grid
            # (x_orig, y_orig) is the lower left point, so the grid center is
            # (x_orig + xlen_half, y_orig + y_len_half)
            grid_center_x = self.x_orig + xlen_half
            grid_center_y = self.y_orig + ylen_half
            lon_half, lat_half = self.iproject(grid_center_x, grid_center_y)
            self.orig_lon = lon_half
            self.orig_lat = lat_half
        self.x_orig = -xlen_half
        self.y_orig = -ylen_half

    def horizontal_rotate(self, angle, fill_value=0.0):
        """
        Rotate the grid horizontally around its center counterclockwise by
        a given angle.

        Parameters
        ----------
        angle : float
            Angle in degrees to rotate the grid counterclockwise.
        fill_value : float, optional
            Value to fill with beyond the grid edge, by default 0.0.

        Returns
        -------
        None

        Note
        ----
        The grid is recentered using the `horizontal_recenter()` method before
        rotation. The rotation is performed using the `scipy.ndimage.rotate()`
        function. The `map_rot` attribute is updated by adding `angle`.
        """
        self.horizontal_recenter()
        self.array = rotate(
            self.array, angle, axes=(1, 0), reshape=True,
            mode='constant', cval=fill_value)
        self.x_orig = -0.5 * self.nx * self.dx
        self.y_orig = -0.5 * self.ny * self.dy
        self.map_rot += angle

    def nudge(self, direction, num_layers=1):
        """
        'Nudge' a grid's dimensions by expanding or contracting in any
        direction, by any number of 2D layers. The output grid is also
        recentered.

        Parameters
        ----------
        direction : char
            Cardinal direction to adjust, either "east", "west", "north",
            "south", as well as "up" or "down", or simply "e", "w", "n", "s",
            "u", "d".

            West-East is axis 0, or the first "X" element of the array, with
            east at index -1.

            North-South is axis 1, or the second "Y" element, with north at
            index 0.

            Up-Down is axis 2, the third "Z" element, with the surface being
            index 0.
        num_layers : int, default 1
            Number of layers to add (positive number) or subtract (negative).
            If positive, the layer values are duplicated from the outermost
            layer.

        Returns
        -------
        None

        Example
        -------
        If a grid ``p_vel`` originally has shape ``(37, 175, 70)``,

        >>> p_vel.nudge('north', 3)

        will add three duplicated 2D layers of ``p_vel.array[:, 0, :]`` to the
        "north" side, giving the new array a shape of ``(37, 178, 70)``.
        2D slices ``p_vel.array[:, 0, :]`` to ``p_vel.array[:, 3, :]`` will be
        identical.
        """
        direction = direction.lower()
        valid_directions = [
            'east', 'west', 'north', 'south', 'up', 'down',
            'e', 'w', 'n', 's', 'u', 'd']
        if direction not in valid_directions:
            valid_directions = ', '.join(f'"{d}"' for d in valid_directions)
            msg = f'Invalid direction "{direction}". Must be one of: '
            msg += valid_directions
            raise ValueError(msg)
        if not isinstance(num_layers, int):
            raise TypeError('num_layers must be an integer')
        if num_layers == 0:
            return
        polarity = -1 if num_layers < 0 else 1

        for _ in range(abs(num_layers)):
            if direction in ['east', 'e']:
                if polarity < 0:
                    self.array = self.array[:-1, :, :]
                else:
                    layer = self.array[-1, :, :]
                    m, n = layer.shape
                    layer = layer.reshape(1, m, n)
                    self.array = np.concatenate((self.array, layer), axis=0)
            elif direction in ['west', 'w']:
                if polarity < 0:
                    self.array = self.array[1:, :, :]
                else:
                    layer = self.array[0, :, :]
                    m, n = layer.shape
                    layer = layer.reshape(1, m, n)
                    self.array = np.concatenate((layer, self.array), axis=0)
            elif direction in ['north', 'n']:
                if polarity < 0:
                    self.array = self.array[:, :-1, :]
                else:
                    layer = self.array[:, -1, :]
                    m, n = layer.shape
                    layer = layer.reshape(m, 1, n)
                    self.array = np.concatenate((self.array, layer), axis=1)
            elif direction in ['south', 's']:
                if polarity < 0:
                    self.array = self.array[:, 1:, :]
                else:
                    layer = self.array[:, 0, :]
                    m, n = layer.shape
                    layer = layer.reshape(m, 1, n)
                    self.array = np.concatenate((layer, self.array), axis=1)
            elif direction in ['up', 'u']:
                if polarity < 0:
                    self.array = self.array[:, :, 1:]
                    self.z_orig += self.dz
                else:
                    layer = self.array[:, :, 0]
                    m, n = layer.shape
                    layer = layer.reshape(m, n, 1)
                    self.array = np.concatenate((layer, self.array), axis=2)
                    self.z_orig -= self.dz
            elif direction in ['down', 'd']:
                if polarity < 0:
                    self.array = self.array[:, :, :-1]
                else:
                    layer = self.array[:, :, -1]
                    m, n = layer.shape
                    layer = layer.reshape(m, n, 1)
                    self.array = np.concatenate((self.array, layer), axis=2)

        if direction not in ['up', 'down', 'u', 'd']:
            self.horizontal_recenter()

    def copy(self):
        """
        Generate a copy of the grid.

        Returns
        -------
        Grid
            A copy of the grid.

        Note
        ----
        The copy is a deep copy, so that the array is copied and not just
        referenced.
        """
        return deepcopy(self)


def main():
    """
    Test code.

    Test 1: generate a gaussian grid and compute the 3D ellipsoid
    around the grid mean.

    Test 2: grid rotation
    """
    # pylint: disable=import-outside-toplevel
    import matplotlib.pyplot as plt

    # http://stackoverflow.com/q/17190649
    def gauss3D(shape=(3, 3, 3), sigmax=0.5, sigmay=0.5, sigmaz=0.5, theta=0):
        """
        Generate a 3D Gaussian kernel.
        """
        m, n, k = [(ss-1.)/2. for ss in shape]
        y, x, z = np.ogrid[-m:m+1, -n:n+1, -k:k+1]
        # xy rotation
        theta = np.radians(theta)
        x2 = math.cos(theta)*x - math.sin(theta)*y
        y2 = math.sin(theta)*x + math.cos(theta)*y
        h = np.exp(-(x2*x2) / (2.*sigmax*sigmax)) *\
            np.exp(-(y2*y2) / (2.*sigmay*sigmay)) *\
            np.exp(-(z*z) / (2.*sigmaz*sigmaz))
        h[h < np.finfo(h.dtype).eps*h.max()] = 0
        sumh = h.sum()
        if sumh != 0:
            h /= sumh
        return h

    # Generate the grid
    nx = 101
    ny = 201
    nz = 11
    x_orig = -50
    y_orig = -100
    grd = NLLGrid(
        nx=nx, ny=ny, nz=nz,
        dx=1, dy=1, dz=1,
        x_orig=x_orig, y_orig=y_orig)
    grd.basename = 'gaussian'
    print(grd, '\n')
    grd.array = gauss3D((nx, ny, nz), 20, 10, 2, 30)

    # Compute statistics
    mean_xyz = grd.get_xyz_mean()
    max_ijk = grd.get_ijk_max()
    max_xyz = grd.get_xyz_max()

    # Plotting
    axes, _cb = grd.plot(max_ijk, handle=True)
    grd.plot_3D_point(axes, mean_xyz, color='g')
    grd.plot_3D_point(axes, max_xyz, color='r')
    grd.plot_ellipsoid(axes, mean_xyz=mean_xyz)
    plt.show()

    # Test recenter and rotation
    grd = NLLGrid()
    grd.proj_name = 'LAMBERT'
    grd.proj_ellipsoid = 'WGS-84'
    grd.x_orig = 30
    grd.y_orig = 15
    grd.orig_lat = 44
    grd.orig_lon = 12
    grd.first_std_paral = 43
    grd.second_std_paral = 45
    grd.basename = 'unrotated'
    grd.array = np.ones((200, 100, 50))
    grd.array[:, 50, :] = 2.
    grd.array[100, :, :] = 3.
    print(grd, '\n')

    axes, _cb = grd.plot(vmin=0, vmax=3, handle=True)
    line_xy = np.vstack((
        np.linspace(25, 100, 10),
        np.linspace(10, 50, 10),
    ))
    line_lonlat = grd.iproject(line_xy[0], line_xy[1])
    axes[0].plot(line_xy[0], line_xy[1], color='k')

    grd_rec = grd.copy()
    grd_rec.horizontal_recenter()
    grd_rec.basename = 'recentered'
    print(grd_rec, '\n')
    axes, _cb = grd_rec.plot(vmin=0, vmax=3, handle=True)
    line_xy = grd_rec.project(line_lonlat[0], line_lonlat[1])
    axes[0].plot(line_xy[0], line_xy[1], color='k')
    line_lonlat = grd_rec.iproject(line_xy[0], line_xy[1])
    line_xy = grd_rec.project(line_lonlat[0], line_lonlat[1])
    axes[0].plot(line_xy[0], line_xy[1], color='r')

    grd_rot = grd.copy()
    rot_angle = 10
    grd_rot.horizontal_rotate(rot_angle)
    grd_rot.basename = f'rotated_{rot_angle}'
    print(grd_rot, '\n')
    axes, _cb = grd_rot.plot(vmin=0, vmax=3, handle=True)
    line_xy = grd_rot.project(line_lonlat[0], line_lonlat[1])
    axes[0].plot(line_xy[0], line_xy[1], color='k')
    line_lonlat = grd_rot.iproject(line_xy[0], line_xy[1])
    line_xy = grd_rot.project(line_lonlat[0], line_lonlat[1])
    axes[0].plot(line_xy[0], line_xy[1], color='r')

    plt.show()


if __name__ == '__main__':
    main()
