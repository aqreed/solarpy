# coding: utf-8

"""
    Tests of the solar panel class
"""


from solarpy import solar_panel
from numpy import array
from datetime import datetime
import unittest as ut


class Test_pvpanel(ut.TestCase):
    """
    Tests exception
    """
    def test_instance_var_exception(self):
        s, eff = 0, 0
        self.assertRaises(ValueError, solar_panel, -1, eff)
        self.assertRaises(ValueError, solar_panel, s, -0.1)
        self.assertRaises(ValueError, solar_panel, s, 1.1)
        self.assertRaises(TypeError, solar_panel, s, eff, id_name=99)

    def test_position_exception(self):
        s, eff = 0, 0
        lat, lng, h = 0, 0, 0
        v = array([0, 0, -1])
        d = datetime(2019, 1, 21, 12, 0)

        sp = solar_panel(s, eff, id_name=None)
        sp.set_orientation(v)
        sp.set_datetime(d)

        with self.assertRaises(ValueError):
            sp.set_position(-91, lng, h)  # latitude exceeds range
            sp.power()

        with self.assertRaises(ValueError):
            sp.set_position(90.1, lng, h)  # latitude exceeds range
            sp.power()

        with self.assertRaises(TypeError):
            sp.set_position('1', lng, h)  # wrong latitude type
            sp.power()

        with self.assertRaises(ValueError):
            sp.set_position(lat, lng, -1)  # altitude exceeds range
            sp.power()

        with self.assertRaises(TypeError):
            sp.set_position(lat, lng, '21')  # altitude wrong type
            sp.power()

    def test_orientation_exception(self):
        s, eff = 0, 0
        lat, lng, h = 0, 0, 0
        d = datetime(2019, 1, 21, 12, 0)

        sp = solar_panel(s, eff, 'panel1')
        sp.set_position(lat, lng, h)
        sp.set_datetime(d)

        with self.assertRaises(ValueError):
            sp.set_orientation(1)  # wrong type
            sp.power()

    def test_datetime_exception(self):
        s, eff = 0, 0
        lat, lng, h = 0, 0, 0
        v = array([0, 0, -1])

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)

        with self.assertRaises(TypeError):
            sp.set_datetime('a')  # wrong type
            sp.power()

    def test_null_surface(self):
        s, eff = 0, 0.5
        lat, lng, h = 0, 0, 0
        v = array([0, 0, -1])
        d = datetime(2019, 1, 21, 12, 0)

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

    def test_zero_efficiency(self):
        s, eff = 1, 0
        lat, lng, h = 0, 0, 0
        v = array([0, 0, -1])
        d = datetime(2019, 1, 21, 12, 0)

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

    def test_upside_down(self):
        s, eff = 1, 0.3
        lat, lng, h = 0, 0, 0
        v = array([0, 0, 1])  # upside down
        d = datetime(2019, 1, 21, 12, 0)

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

    def test_night(self):
        s, eff = 1, 0.3
        lat, lng, h = 0, 0, 0
        v = array([0, 0, -1])

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)

        d = datetime(2019, 1, 21, 22, 15)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

        d = datetime(2019, 7, 21, 0, 0)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

        d = datetime(2019, 11, 1, 3, 40)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

        d = datetime(2019, 4, 12, 5, 0)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

    def test_permanent_darkness(self):
        # north hemisphere
        s, eff = 1, 0.3
        lat, lng, h = 80, 0, 0
        v = array([0, 0, -1])
        d = datetime(2019, 1, 15, 12, 0)  # Jan 15

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

        # south hemisphere
        s, eff = 1, 0.3
        lat, lng, h = -85, 0, 0
        v = array([0, 0, -1])
        d = datetime(2019, 7, 1, 12, 0)  # Jul 1

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

    def test_pointing_opposite_direction(self):
        # panel points west in the morning
        s, eff = 1, 0.3
        lat, lng, h = 0, 0, 0
        v = array([0, -1, 0])  # west
        d = datetime(2019, 8, 1, 10, 0)

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)

        # panel points east in the afternoon
        s, eff = 1, 0.3
        lat, lng, h = 0, 0, 0
        v = array([0, 1, 0])  # east
        d = datetime(2019, 8, 1, 18, 0)

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)
        sp.set_datetime(d)
        self.assertAlmostEqual(sp.power(), 0)
