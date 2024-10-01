"""
Tests for the NLLGrid class.
"""
import sys
import os
import unittest
import tempfile
import numpy as np
from numpy.testing import assert_array_equal
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(dir_path, '..'))
from nllgrid import NLLGrid  # noqa


class TestNLLGrid(unittest.TestCase):
    """
    Tests for the NLLGrid class.
    """
    def __init__(self, *args, **kwargs):
        super(TestNLLGrid, self).__init__(*args, **kwargs)
        self.grd_bname = os.path.join(dir_path, 'data', 'layer.P.mod')

    def test_read(self):
        """
        Test reading a NLL grid file (hdr and buf)
        """
        grd = NLLGrid(self.grd_bname)
        self.assertEqual(grd.basename, self.grd_bname)
        self.assertEqual(grd.nx, 58)
        self.assertEqual(grd.ny, 38)
        self.assertEqual(grd.nz, 20)
        self.assertEqual(grd.x_orig, -146.0)
        self.assertEqual(grd.y_orig, -95.0)
        self.assertEqual(grd.z_orig, -2.0)
        self.assertEqual(grd.dx, 5.0)
        self.assertEqual(grd.dy, 5.0)
        self.assertEqual(grd.dz, 5.0)
        self.assertEqual(grd.type, 'SLOW_LEN')
        self.assertEqual(grd.float_type, 'FLOAT')
        self.assertEqual(grd.proj_name, 'LAMBERT')
        self.assertEqual(grd.proj_ellipsoid, 'Clarke-1880')
        self.assertEqual(grd.orig_lon, 132.988)
        self.assertEqual(grd.orig_lat, 33.618)
        self.assertEqual(grd.first_std_paral, 32.0)
        self.assertEqual(grd.second_std_paral, 35.0)
        self.assertEqual(grd.map_rot, 0.0)
        self.assertAlmostEqual(grd.array[1, 2, 3], 0.8333333, places=5)
        self.assertAlmostEqual(grd.array.sum(), 36733.33, places=2)

    def test_write(self):
        """
        Test writing a NLL grid file (hdr and buf).
        """
        grd = NLLGrid(self.grd_bname)
        # write to a temporary file using tempfile module
        with tempfile.NamedTemporaryFile() as tmp:
            grd.write_hdr_file(tmp.name)
            grd.write_buf_file(tmp.name)
            grd2 = NLLGrid(tmp.name)
        self.assertEqual(grd2.nx, grd.nx)
        self.assertEqual(grd2.ny, grd.ny)
        self.assertEqual(grd2.nz, grd.nz)
        self.assertEqual(grd2.x_orig, grd.x_orig)
        self.assertEqual(grd2.y_orig, grd.y_orig)
        self.assertEqual(grd2.z_orig, grd.z_orig)
        self.assertEqual(grd2.dx, grd.dx)
        self.assertEqual(grd2.dy, grd.dy)
        self.assertEqual(grd2.dz, grd.dz)
        self.assertEqual(grd2.type, grd.type)
        self.assertEqual(grd2.float_type, grd.float_type)
        self.assertEqual(grd2.proj_name, grd.proj_name)
        self.assertEqual(grd2.proj_ellipsoid, grd.proj_ellipsoid)
        self.assertEqual(grd2.orig_lon, grd.orig_lon)
        self.assertEqual(grd2.orig_lat, grd.orig_lat)
        self.assertEqual(grd2.first_std_paral, grd.first_std_paral)
        self.assertEqual(grd2.second_std_paral, grd.second_std_paral)
        self.assertEqual(grd2.map_rot, grd.map_rot)
        self.assertEqual(grd2.array.shape, grd.array.shape)
        self.assertEqual(grd2.array.dtype, grd.array.dtype)
        self.assertTrue((grd2.array == grd.array).all())

    def test_create(self):
        """
        Test creating a NLL grid file (hdr and buf) and writing it.
        """
        grd = NLLGrid()
        grd.array = np.random.rand(10, 20, 30).astype(np.float32)
        grd.dx = 0.5
        grd.dy = 0.5
        grd.dz = 0.5
        grd.x_orig = -10.
        grd.y_orig = -20.
        grd.z_orig = -1.
        grd.type = 'VELOCITY'
        grd.orig_lat = 40.63
        grd.orig_lon = 15.80
        grd.proj_name = 'LAMBERT'
        grd.first_std_paral = 38.
        grd.second_std_paral = 42.
        grd.proj_ellipsoid = 'WGS-84'
        with tempfile.NamedTemporaryFile() as tmp:
            grd.write_hdr_file(tmp.name)
            grd.write_buf_file(tmp.name)
            grd2 = NLLGrid(tmp.name)
        self.assertEqual(grd2.nx, grd.nx)
        self.assertEqual(grd2.ny, grd.ny)
        self.assertEqual(grd2.nz, grd.nz)
        self.assertEqual(grd2.x_orig, grd.x_orig)
        self.assertEqual(grd2.y_orig, grd.y_orig)
        self.assertEqual(grd2.z_orig, grd.z_orig)
        self.assertEqual(grd2.dx, grd.dx)
        self.assertEqual(grd2.dy, grd.dy)
        self.assertEqual(grd2.dz, grd.dz)
        self.assertEqual(grd2.type, grd.type)
        self.assertEqual(grd2.float_type, grd.float_type)
        self.assertEqual(grd2.proj_name, grd.proj_name)
        self.assertEqual(grd2.proj_ellipsoid, grd.proj_ellipsoid)
        self.assertEqual(grd2.orig_lon, grd.orig_lon)
        self.assertEqual(grd2.orig_lat, grd.orig_lat)
        self.assertEqual(grd2.first_std_paral, grd.first_std_paral)
        self.assertEqual(grd2.second_std_paral, grd.second_std_paral)
        self.assertEqual(grd2.map_rot, grd.map_rot)
        self.assertEqual(grd2.array.shape, grd.array.shape)
        self.assertEqual(grd2.array.dtype, grd.array.dtype)
        self.assertAlmostEqual(grd2.array.sum(), grd.array.sum(), places=5)

    def test_nudge(self):
        """Test nudging a grid."""
        grd = NLLGrid(self.grd_bname)
        grd.array = np.random.rand(*grd.array.shape)
        n_east = 5
        n_west = 3
        n_north = 8
        n_south = 2
        n_up = 7
        n_down = 12
        orig_array = grd.array.copy()
        shape = grd.array.shape
        grd.nudge('east', n_east)
        for n in range(n_east):
            self.assertEqual(grd.array.shape[0], shape[0]+n_east)
            idx = shape[0] - 1
            assert_array_equal(grd.array[idx, :, :], grd.array[idx+n, :, :])
        shape = grd.array.shape
        grd.nudge('west', n_west)
        for n in range(n_west):
            self.assertEqual(grd.array.shape[0], shape[0]+n_west)
            idx = 0
            assert_array_equal(grd.array[idx, :, :], grd.array[idx+n, :, :])
        shape = grd.array.shape
        grd.nudge('north', n_north)
        for n in range(n_north):
            self.assertEqual(grd.array.shape[1], shape[1]+n_north)
            idx = shape[1] - 1
            assert_array_equal(grd.array[:, idx, :], grd.array[:, idx+n, :])
        shape = grd.array.shape
        grd.nudge('south', n_south)
        for n in range(n_south):
            self.assertEqual(grd.array.shape[1], shape[1]+n_south)
            idx = 0
            assert_array_equal(grd.array[:, idx, :], grd.array[:, idx+n, :])
        shape = grd.array.shape
        z_orig = grd.z_orig
        grd.nudge('up', n_up)
        for n in range(n_up):
            self.assertEqual(grd.array.shape[2], shape[2]+n_up)
            self.assertEqual(grd.z_orig, z_orig - n_up * grd.dz)
            idx = 0
            assert_array_equal(grd.array[:, :, idx], grd.array[:, :, idx+n])
        shape = grd.array.shape
        grd.nudge('down', n_down)
        for n in range(n_down):
            self.assertEqual(grd.array.shape[2], shape[2]+n_down)
            idx = shape[2] - 1
            assert_array_equal(grd.array[:, :, idx], grd.array[:, :, idx+n])
        grd.nudge('east', -n_east)
        grd.nudge('west', -n_west)
        grd.nudge('north', -n_north)
        grd.nudge('south', -n_south)
        grd.nudge('up', -n_up)
        grd.nudge('down', -n_down)
        assert_array_equal(grd.array, orig_array)


if __name__ == '__main__':
    unittest.main()
