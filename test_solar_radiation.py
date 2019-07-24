# coding: utf-8

"""
    Test for the solar radiation model, based on Duffie, J.A., and
    Beckman, W. A., 1974, "Solar energy thermal processes"
"""

import solar_radiation as sr
from solar_radiation import NoSunsetNoSunrise

import numpy as np
from numpy.testing import (assert_equal, assert_almost_equal,
                           assert_array_equal, assert_array_almost_equal)
from scipy.optimize import fmin
import unittest as ut


class Test_Ranges(ut.TestCase):
    """
    Tests ranges checks
    """
    def test_nth_day_range(self):
        self.assertRaises(ValueError, sr.check_nth_day_range, 0)
        self.assertRaises(ValueError, sr.check_nth_day_range, 366)

    def test_latitude_range(self):
        self.assertRaises(ValueError, sr.check_latitude_range, -91)
        self.assertRaises(ValueError, sr.check_latitude_range, 91)

    def test_longitude_range(self):
        self.assertRaises(ValueError, sr.check_longitude_range, -1)
        self.assertRaises(ValueError, sr.check_longitude_range, 360)


class Test_day_of_the_year(ut.TestCase):
    """
    Tests day of the year values.
    """
    def test_Jan1(self):
        month, day = 1, 1
        expected_value = 1  # January 1
        self.assertEqual(sr.day_of_the_year(month, day), expected_value)

    def test_Feb1(self):
        month, day = 2, 1
        expected_value = 32  # February 1
        self.assertEqual(sr.day_of_the_year(month, day), expected_value)

    def test_summerSolstice(self):
        month, day = 6, 20
        expected_value = 171  # summer solstice
        self.assertEqual(sr.day_of_the_year(month, day), expected_value)

    def test_Dec31(self):
        month, day = 12, 31
        expected_value = 365  # December 31
        self.assertEqual(sr.day_of_the_year(month, day), expected_value)


def test_B_nth_day():
    """
    Tests B(n) values.
    """
    n = 1
    expected_value = 0
    assert_equal(sr.B_nth_day(n), expected_value)

    n = 365
    expected_value = 6.2659711
    assert_almost_equal(sr.B_nth_day(n), expected_value, 6)


class Test_Gon(ut.TestCase):
    """
    Tests radiation on a plane normal values.
    """
    def test_min(self):
        a = np.arange(1, 366)  # array with all days
        self.assertTrue((sr.Gon(a) > 1320).all())

    def test_max(self):
        a = np.arange(1, 366)  # array with all days
        self.assertTrue((sr.Gon(a) < 1420).all())

    def test_Jan1(self):
        n = 1  # January 1
        expected_value = 1415
        self.assertAlmostEqual(sr.Gon(n), expected_value, delta=1)

    def test_summerSolstice(self):
        n = 171  # July 21
        expected_value = 1322
        self.assertAlmostEqual(sr.Gon(n), expected_value, delta=1)


class Test_Eq_time(ut.TestCase):
    """
    Tests equation of time values.
    """
    def test_min(self):
        a = np.arange(1, 366)  # array with all days
        self.assertTrue((sr.Eq_time(a) > -15).all())

    def test_max(self):
        a = np.arange(1, 366)  # array with all days
        self.assertTrue((sr.Eq_time(a) < 17).all())

    def test_Jan1(self):
        n = 1  # January 1
        expected_value = -3
        self.assertAlmostEqual(sr.Eq_time(n), expected_value, delta=1)

    def test_summerSolstice(self):
        n = 171  # July 21
        expected_value = -1
        self.assertAlmostEqual(sr.Eq_time(n), expected_value, delta=1)


class Test_declination(ut.TestCase):
    """
    Tests equation of time values.
    """
    def test_min(self):
        a = np.arange(1, 366)  # array with all days
        self.assertTrue((np.rad2deg(sr.declination(a)) > -24).all())

    def test_max(self):
        a = np.arange(1, 366)  # array with all days
        self.assertTrue((np.rad2deg(sr.declination(a)) < 24).all())

    def test_Jan1(self):
        n = 1  # January 1
        expected_value = np.deg2rad(-23)
        self.assertAlmostEqual(sr.declination(n), expected_value, 1)

    def test_summerSolstice(self):
        n = 171  # July 21
        expected_value = np.deg2rad(23)
        self.assertAlmostEqual(sr.declination(n), expected_value, 1)

    def test_Feb13(self):  # Example 1.6.1
        n = 44  # February 13
        expected_value = np.deg2rad(-14)
        self.assertAlmostEqual(sr.declination(n), expected_value, 1)

    def test_Mar16(self):  # Example 1.6.3
        n = sr.day_of_the_year(3, 16)  # March 16
        expected_value = np.deg2rad(-2.4)
        self.assertAlmostEqual(sr.declination(n), expected_value, 1)


class Test_solar_time(ut.TestCase):
    """
    Tests solar time function. Values from Duffie and Beckman example 1.5.1
    """
    def test_errors(self):
        n = 34
        lng = 89.4
        self.assertRaises(ValueError, sr.solar_time, n, -1, 0, lng)
        self.assertRaises(ValueError, sr.solar_time, n, 24, 0, lng)
        self.assertRaises(ValueError, sr.solar_time, n, 0, -1, lng)
        self.assertRaises(ValueError, sr.solar_time, n, 0, 60, lng)

    def test_Feb3(self):
        n = 34
        t_std_h, t_std_min = 10, 30  # standard time (hour, minute)
        lng = 89.4

        expected_value = (10, 18, 54)
        self.assertEqual(sr.solar_time(n, t_std_h, t_std_min, lng),
                         expected_value)


class Test_hour_angle(ut.TestCase):
    """
    Tests hour angle function. Values from Duffie and Beckman
    """
    def test_errors(self):
        self.assertRaises(ValueError, sr.hour_angle, -1, 0)
        self.assertRaises(ValueError, sr.hour_angle, 24, 0)
        self.assertRaises(ValueError, sr.hour_angle, 0, -1)
        self.assertRaises(ValueError, sr.hour_angle, 0, 60)

    def test_examples(self):
        # noon
        hour, minute = 12, 0
        expected_value = np.deg2rad(0)
        self.assertEqual(sr.hour_angle(hour, minute), expected_value)

        # Example 1.6.1
        hour, minute = 10, 30
        expected_value = np.deg2rad(-22.5)
        self.assertEqual(sr.hour_angle(hour, minute), expected_value)

        # Example 1.6.2a
        hour, minute = 9, 30
        expected_value = np.deg2rad(-37.5)
        self.assertEqual(sr.hour_angle(hour, minute), expected_value)

        # Example 1.6.2b
        hour, minute = 18, 30
        expected_value = np.deg2rad(97.5)
        self.assertEqual(sr.hour_angle(hour, minute), expected_value)

        # Example 1.6.3
        hour, minute = 16, 0
        expected_value = np.deg2rad(60)
        self.assertEqual(sr.hour_angle(hour, minute), expected_value)

        # Example 1.7.1
        hour, minute = 14, 0
        expected_value = np.deg2rad(30)
        self.assertEqual(sr.hour_angle(hour, minute), expected_value)


def test_angle_of_incidence():
    """
    Tests angle of incidence function. Values from Duffie and Beckman
    example 1.6.1
    """
    n = 44
    lat = 43
    beta = 45
    surf_az = 15
    hour, minute = 10, 30

    expected_value = np.deg2rad(35)
    assert_almost_equal(sr.theta(n, lat, beta, surf_az, hour, minute),
                        expected_value, decimal=3)


def test_zenith_angle():
    """
    Tests zenith angle function. Values from Duffie and Beckman
    """
    # noon at summer solstice and lat = 23.4º
    n = 171
    lat = 23.45
    hour, minute = 12, 0

    expected_value = np.deg2rad(0)
    assert_almost_equal(sr.theta_z(n, lat, hour, minute),
                        expected_value, decimal=2)

    # Example 1.6.2a
    n = 44
    lat = 43
    hour, minute = 9, 30

    expected_value = np.deg2rad(66.5)
    assert_almost_equal(sr.theta_z(n, lat, hour, minute),
                        expected_value, decimal=2)

    # Example 1.6.2b
    n = sr.day_of_the_year(7, 1)
    lat = 43
    hour, minute = 18, 30

    expected_value = np.deg2rad(79.6)
    assert_almost_equal(sr.theta_z(n, lat, hour, minute),
                        expected_value, decimal=2)

    # Example 1.6.3
    n = sr.day_of_the_year(3, 16)
    lat = 43
    hour, minute = 16, 0

    expected_value = np.deg2rad(70.3)
    assert_almost_equal(sr.theta_z(n, lat, hour, minute),
                        expected_value, decimal=2)


def test_solar_azimuth():
    """
    Tests solar azimuth angle function. Values from Duffie and Beckman
    """
    # different values of (n, lat) at noon
    hour, minute = 12, 0

    expected_value = np.deg2rad(0)
    n, lat = 1, 0
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)
    n, lat = 90, 30
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)
    n, lat = 180, 60
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)
    n, lat = 270, 90
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)

    expected_value = np.deg2rad(180)
    n, lat = 90, -30
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)
    n, lat = 180, -60
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)
    n, lat = 270, -90
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)

    # Example 1.6.2a
    n = 44
    lat = 43
    hour, minute = 9, 30

    expected_value = np.deg2rad(-40.0)
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)

    # Example 1.6.2b
    n = sr.day_of_the_year(7, 1)
    lat = 43
    hour, minute = 18, 30

    expected_value = np.deg2rad(112.0)
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)

    # Example 1.6.3
    n = sr.day_of_the_year(3, 16)
    lat = 43
    hour, minute = 16, 0

    expected_value = np.deg2rad(66.8)
    assert_almost_equal(sr.solar_azimuth(n, lat, hour, minute),
                        expected_value, decimal=2)


def test_solar_altitude():
    """
    Tests solar azimuth angle function. Values from Duffie and Beckman
    """
    # Example 1.6.3
    n = sr.day_of_the_year(3, 16)
    lat = 43
    hour, minute = 16, 0

    expected_value = np.deg2rad(19.7)
    assert_almost_equal(sr.solar_altitude(n, lat, hour, minute),
                        expected_value, decimal=2)


class Test_sunset_hour_angle(ut.TestCase):
    """
    Tests sunset hour angle function. Values from Duffie and Beckman
    """
    def test_errors(self):
        self.assertRaises(NoSunsetNoSunrise, sr.sunset_hour_angle, 1, 80)
        self.assertRaises(NoSunsetNoSunrise, sr.sunset_hour_angle, 171, -75)

    def test_examples(self):
        # Example 1.6.3
        n = sr.day_of_the_year(3, 16)
        lat = 43
        expected_value = np.deg2rad(87.8)
        self.assertAlmostEqual(sr.sunset_hour_angle(n, lat), expected_value, 1)


class Test_sunset_time(ut.TestCase):
    """
    Tests sunset time function. Values from Duffie and Beckman
    """
    def test_examples(self):
        # Example 1.6.3
        n = sr.day_of_the_year(3, 16)
        lat = 43
        expected_value = 17
        self.assertAlmostEqual(sr.sunset_time(n, lat).hour,
                               expected_value, 2)
        expected_value = 52  # diferent year than boook!
        self.assertAlmostEqual(sr.sunset_time(n, lat).minute,
                               expected_value, 2)


class Test_sunrise_hour_angle(ut.TestCase):
    """
    Tests sunrise hour angle function. Values from Duffie and Beckman
    """
    def test_errors(self):
        self.assertRaises(NoSunsetNoSunrise, sr.sunrise_hour_angle, 1, 80)
        self.assertRaises(NoSunsetNoSunrise, sr.sunrise_hour_angle, 171, -75)

    def test_examples(self):
        # Example 1.6.3
        n = sr.day_of_the_year(3, 16)
        lat = 43
        expected_value = np.deg2rad(-87.8)
        self.assertAlmostEqual(sr.sunrise_hour_angle(n, lat),
                               expected_value, 1)


class Test_sunrise_time(ut.TestCase):
    """
    Tests sunrise time function. Values from Duffie and Beckman
    """
    def test_examples(self):
        # Example 1.6.3
        n = sr.day_of_the_year(3, 16)
        lat = 43
        expected_value = 6
        self.assertAlmostEqual(sr.sunrise_time(n, lat).hour,
                               expected_value, 2)
        expected_value = 7  # diferent year than boook!
        self.assertAlmostEqual(sr.sunrise_time(n, lat).minute,
                               expected_value, 2)
