# coding: utf-8

"""
    Test for the solar radiation model, based on Duffie, J.A., and
    Beckman, W. A., 1974, "Solar energy thermal processes"
"""

import solar_radiation as sr

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
        n = fmin(sr.Gon, 150)[0]  # gets "n" for the min value
        min_value = 1320
        self.assertTrue(sr.Gon(n) > min_value)

    def test_max(self):
        n = fmin(lambda x: -sr.Gon(x), 1)  # gets "n" for the max value
        max_value = 1420
        self.assertTrue(sr.Gon(n) < max_value)

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
        self.assertAlmostEqual(sr.declination(n), expected_value, delta=1)

    def test_summerSolstice(self):
        n = 171  # July 21
        expected_value = np.deg2rad(23)
        self.assertAlmostEqual(sr.declination(n), expected_value, delta=1)
