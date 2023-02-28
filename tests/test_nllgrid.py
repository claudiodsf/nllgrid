"""
Tests for the NLLGrid class.
"""
import sys
import os
import unittest
import tempfile
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(dir_path, '..'))
from nllgrid import NLLGrid  # noqa


class TestNLLGrid(unittest.TestCase):
    """
    Tests for the NLLGrid class.
    """

    def test_read(self):
        """
        Test reading a NLL grid file (hdr and buf)
        """
        grd_bname = os.path.join(dir_path, 'data', 'layer.P.mod')
        grd = NLLGrid(grd_bname)
        self.assertEqual(grd.basename, grd_bname)
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
        grd_bname = os.path.join(dir_path, 'data', 'layer.P.mod')
        grd = NLLGrid(grd_bname)
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


if __name__ == '__main__':
    unittest.main()
