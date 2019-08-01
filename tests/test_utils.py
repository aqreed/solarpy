# coding: utf-8

"""
    Test of the utilities of the solar radiation model"
"""


from solar_radiation import *

import numpy as np
from numpy.testing import (assert_equal, assert_almost_equal,
                           assert_array_almost_equal)
import unittest as ut


class Test_Ranges(ut.TestCase):
    """
    Tests ranges checks
    """
    def test_nth_day_range(self):
        self.assertRaises(ValueError, check_day_range, 0)
        self.assertRaises(ValueError, check_day_range, 366)

    def test_latitude_range(self):
        self.assertRaises(ValueError, check_lat_range, -91)
        self.assertRaises(ValueError, check_lat_range, 91)

    def test_longitude_range(self):
        self.assertRaises(ValueError, check_long_range, -1)
        self.assertRaises(ValueError, check_long_range, 360)


def test_lla2ecef():
    """
    Test function that returns ecef position from lat, long, altitude
    """
    a = 6378137  # [m] Earth equatorial axis
    b = 6356752.3142  # [m] Earth polar axis

    # OX-axis
    lat = 0
    lng = 0
    h = 0
    expected_value = np.array([a, 0, 0])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    lat = 0
    lng = 180
    h = 0
    expected_value = np.array([-a, 0, 0])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    # OY-axis
    lat = 0
    lng = 90
    h = 0
    expected_value = np.array([0, a, 0])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    lat = 0
    lng = 270
    h = 0
    expected_value = np.array([0, -a, 0])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    # OZ-axis
    lat = 90
    lng = 0
    h = 0
    expected_value = np.array([0, 0, b])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    lat = -90
    lng = 0
    h = 0
    expected_value = np.array([0, 0, -b])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)


def test_ned2ecef():
    """
    Test function that transforms ned-basis vectors to ecef-basis
    """
    lat, lng = 0, 0

    v_ned = np.array([1, 0, 0])
    expected_value = np.array([0, 0, 1])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = np.array([0, 1, 0])
    expected_value = np.array([0, 1, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = np.array([0, 0, 1])
    expected_value = np.array([-1, 0, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    lat, lng = 0, 90

    v_ned = np.array([1, 0, 0])
    expected_value = np.array([0, 0, 1])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = np.array([0, 1, 0])
    expected_value = np.array([-1, 0, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = np.array([0, 0, 1])
    expected_value = np.array([0, -1, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    lat, lng = 90, 0

    v_ned = np.array([1, 0, 0])
    expected_value = np.array([-1, 0, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = np.array([0, 1, 0])
    expected_value = np.array([0, 1, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = np.array([0, 0, 1])
    expected_value = np.array([0, 0, -1])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)
