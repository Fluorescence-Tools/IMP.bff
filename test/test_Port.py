import utils
import os
import unittest
import json
import numpy as np

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import IMP.bff
from constants import *

class Tests(unittest.TestCase):

    def test_port_init_single(self):
        v1 = 23.0
        v2 = 29.0
        p1 = IMP.bff.Port(v1)
        # check setting of value
        p2 = IMP.bff.Port()
        p2.value = v1
        # check fixing
        p3 = IMP.bff.Port(
            value=v1,
            fixed=True
        )
        p4 = IMP.bff.Port(
            value=v1,
            fixed=False
        )
        # check linking
        p5 = IMP.bff.Port(v2)
        p5.link = p4
        np.testing.assert_allclose(p1.value, p2.value)
        self.assertEqual(p3.fixed, True)
        self.assertEqual(p4.fixed, False)
        np.testing.assert_allclose(p5.value, p4.value)

    def test_port_init_singelton(self):
        """Test chinet Port class set_value and get_value"""
        v1 = 23.0
        v2 = 29.0
        p1 = IMP.bff.Port()
        p1.value = v1
        # check setting of value
        p2 = IMP.bff.Port(v1)
        # check fixing
        p3 = IMP.bff.Port(
            value=v1,
            fixed=True
        )
        p4 = IMP.bff.Port(
            value=v1,
            fixed=False
        )
        # check linking
        p5 = IMP.bff.Port(v2)
        p5.link = p4

        np.testing.assert_allclose(p1.value, p2.value)
        self.assertEqual(p3.fixed, True)
        self.assertEqual(p4.fixed, False)
        np.testing.assert_allclose(p5.value, p4.value)

        # check bounds
        fixed = False
        is_output = False
        is_reactive = False
        is_bounded = True
        lower_bound = 2
        upper_bound = 5
        value = 0
        p6 = IMP.bff.Port(
            value=value,
            fixed=fixed,
            is_output=is_output,
            is_reactive=is_reactive,
            is_bounded=is_bounded,
            lb=lower_bound,
            ub=upper_bound
        )
        self.assertEqual(
            np.all(p6.value <= upper_bound),
            True
        )
        self.assertEqual(
            np.all(p6.value >= lower_bound),
            True
        )
        self.assertAlmostEqual(
            p6.value,
            lower_bound  # the lower bound is not part of the
        )

    def test_port_bounds(self):
        """Test IMP.bff.Port class set_value and get_value"""
        v1 = np.array([1, 2, 3, 6, 5.5, -3, -2, -6.1, -10000, 10000], dtype=np.double)
        p1 = IMP.bff.Port()
        p1.value = v1

        np.testing.assert_allclose(p1.value, v1)

        p1.bounds = 0, 1
        p1.bounded = True

        self.assertEqual((p1.value <= 1).all(), True)
        self.assertEqual((p1.value >= 0).all(), True)

    def test_port_get_set_value(self):
        """Test IMP.bff Port class set_value and get_value"""
        v1 = [1, 2, 3]
        p1 = IMP.bff.Port()
        p1.value = v1

        p2 = IMP.bff.Port()
        p2.value = v1
        self.assertEqual(
            (p1.value == p2.value).all(),
            True
        )

    def test_port_init_vector(self):
        """Test IMP.bff Port class set_value and get_value"""
        v1 = [1, 2, 3, 5, 8]
        v2 = [1, 2, 4, 8, 16]
        # check setting of value
        p1 = IMP.bff.Port()
        p1.value = v1
        p2 = IMP.bff.Port(v1)
        self.assertListEqual(
            list(p2.value),
            list(p1.value)
        )
        # check fixing
        p3 = IMP.bff.Port(v1, True)
        p4 = IMP.bff.Port(v1, False)
        self.assertEqual(p3.fixed, True)
        self.assertEqual(p4.fixed, False)

        # check linking
        p5 = IMP.bff.Port(v2, False)
        p5.link = p4
        np.testing.assert_allclose(p5.value, p4.value)

    def test_port_init_array(self):
        """Test IMP.bff Port class set_value and get_value"""
        array = np.array([1, 2, 3, 5, 8, 13], dtype=np.double)
        p1 = IMP.bff.Port()
        p1.value = array
        p2 = IMP.bff.Port(array)
        np.testing.assert_allclose(p1.value, p2.value)

    def test_set_get_value_1(self):
        """Test IMP.bff Port class set_value and get_value"""
        value = 23.0
        port = IMP.bff.Port(value)
        self.assertEqual(port.value, value)

    def test_set_get_value_2(self):
        """Test IMP.bff Port class set_value and get_value"""
        value = (1,)
        port = IMP.bff.Port()
        port.value = value
        self.assertEqual(port.value, value)

    def test_port_link_value(self):
        value1 = np.array([12], dtype=np.double)
        value2 = np.array([6], dtype=np.double)
        p1 = IMP.bff.Port(value1)
        p2 = IMP.bff.Port(value2)
        np.testing.assert_allclose(p1.value, value1)
        np.testing.assert_allclose(p2.value, value2)
        self.assertEqual(np.allclose(p1.value, p2.value), False)

        p2.link = p1
        np.testing.assert_allclose(p1.value, p2.value)
        p2.unlink()
        np.testing.assert_allclose(p2.value, value2)

    def test_port_fixed(self):
        p1 = IMP.bff.Port(12)
        p1.fixed = True
        self.assertEqual(p1.fixed, True)

        p1.fixed = False
        self.assertEqual(p1.fixed, False)

    def test_port_reactive(self):
        p1 = IMP.bff.Port(12)
        p1.reactive = True
        self.assertEqual(p1.reactive, True)

        p1.reactive = False
        self.assertEqual(p1.reactive, False)

    @unittest.skipUnless(CONNECTS, "Could not connect to DB")
    def test_db_write(self):
        value_array = (1, 2, 3, 5, 8, 13)
        port = IMP.bff.Port(
            value=value_array,
            fixed=True
        )
        connect_success = port.connect_to_db(**DB_DICT)
        write_success = port.write_to_db()

        self.assertEqual(connect_success, True)
        self.assertEqual(write_success, True)

    @unittest.skipUnless(CONNECTS, "Could not connect to DB")
    def test_port_db_restore(self):
        value_array = (1, 2, 3, 5, 8, 13)
        value = 17

        port = IMP.bff.Port()
        port.value = value
        port.value = value_array

        port.connect_to_db(**DB_DICT)
        port.write_to_db()

        port_reload = IMP.bff.Port()
        port_reload.connect_to_db(**DB_DICT)
        self.assertEqual(port_reload.read_from_db(port.oid), True)

        dict_port = json.loads(port.get_json())
        dict_port_restore = json.loads(port.get_json())

        self.assertEqual(dict_port, dict_port_restore)


if __name__ == '__main__':
    unittest.main()
