# coding: utf-8

"""
    Test of the utilities of the solar radiation model"
"""


from solarpy.utils import *
from numpy import sin, cos, deg2rad, rad2deg, array
from numpy.testing import assert_array_almost_equal
import unittest as ut


class Test_ranges(ut.TestCase):
    """
    Tests ranges checks
    """
    def test_latitude_range(self):
        self.assertRaises(ValueError, check_lat_range, -181)
        self.assertRaises(ValueError, check_lat_range, 181.0)
        self.assertRaises(ValueError, check_lat_range, array([-215, 2, 55]))

    def test_longitude_range(self):
        self.assertRaises(ValueError, check_long_range, -181)
        self.assertRaises(ValueError, check_long_range, 181.0)
        self.assertRaises(ValueError, check_long_range, array([326, -180]))

    def test_altitude_range(self):
        self.assertRaises(ValueError, check_alt_range, -1)
        self.assertRaises(ValueError, check_alt_range, 24001.0)
        self.assertRaises(ValueError, check_alt_range, array([150, 101800]))


class Test_day_of_the_year(ut.TestCase):
    """
    Tests day of the year values.
    """
    def test_type(self):
        self.assertRaises(TypeError, day_of_the_year, 1)
        self.assertRaises(TypeError, day_of_the_year, [1, 2])
        self.assertRaises(TypeError, day_of_the_year, 'a')

    def test_Jan1(self):
        date = datetime(2019, 1, 1)
        expected_value = 1  # January 1
        self.assertEqual(day_of_the_year(date), expected_value)

    def test_Feb1(self):
        date = datetime(2019, 2, 1)
        expected_value = 32  # February 1
        self.assertEqual(day_of_the_year(date), expected_value)

    def test_summerSolstice(self):
        date = datetime(2019, 6, 20)
        expected_value = 171  # summer solstice
        self.assertEqual(day_of_the_year(date), expected_value)

    def test_Dec31(self):
        date = datetime(2019, 12, 31)
        expected_value = 365  # December 31
        self.assertEqual(day_of_the_year(date), expected_value)

    def test_array(self):
        date = array([datetime(2019, 12, 29),
                      datetime(2019, 12, 30),
                      datetime(2019, 12, 31)])
        expected_value = array([363, 364, 365])  # December 31
        self.assertTrue((day_of_the_year(date) == expected_value).all())


class Test_exception(ut.TestCase):
    """
    Tests customised exception
    """
    def test_msg(self):
        with self.assertRaises(NoSunsetNoSunrise) as error:
            raise NoSunsetNoSunrise

        self.assertEqual(
            "Permanent night (or day) on this latitude on this day",
            error.exception.msg)


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
    expected_value = array([a, 0, 0])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    lat = 0
    lng = 180
    h = 0
    expected_value = array([-a, 0, 0])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    # OY-axis
    lat = 0
    lng = 90
    h = 0
    expected_value = array([0, a, 0])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    lat = 0
    lng = -90
    h = 0
    expected_value = array([0, -a, 0])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    # OZ-axis
    lat = 90
    lng = 0
    h = 0
    expected_value = array([0, 0, b])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)

    lat = -90
    lng = 0
    h = 0
    expected_value = array([0, 0, -b])
    assert_array_almost_equal(lla2ecef(lat, lng, h), expected_value, 4)


def test_ned2ecef():
    """
    Test function that transforms ned-basis vectors to ecef-basis
    """
    lat, lng = 0, 0

    v_ned = array([1, 0, 0])
    expected_value = array([0, 0, 1])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = array([0, 1, 0])
    expected_value = array([0, 1, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = array([0, 0, 1])
    expected_value = array([-1, 0, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    lat, lng = 0, 90

    v_ned = array([1, 0, 0])
    expected_value = array([0, 0, 1])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = array([0, 1, 0])
    expected_value = array([-1, 0, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = array([0, 0, 1])
    expected_value = array([0, -1, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    lat, lng = 90, 0

    v_ned = array([1, 0, 0])
    expected_value = array([-1, 0, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = array([0, 1, 0])
    expected_value = array([0, 1, 0])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)

    v_ned = array([0, 0, 1])
    expected_value = array([0, 0, -1])
    assert_array_almost_equal(ned2ecef(v_ned, lat, lng), expected_value)


def test_pressure():
    """
    Test pressure function
    """
    h = 0
    expected_value = 101325
    assert_array_almost_equal(pressure(h), expected_value)

    h = 20e3
    expected_value = 5529
    assert_array_almost_equal(pressure(h), expected_value)
