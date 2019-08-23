# coding: utf-8

"""
    Tests of the solar panel class
"""


from solarpy.pvpanel import solar_panel
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

    def test_position_exception(self):
        s, eff = 0, 0
        lat, lng, h = 0, 0, 0
        v = array([0, 0, -1])
        d = datetime(2019, 1, 21, 12, 0)

        sp = solar_panel(s, eff)
        sp.set_orientation(v)
        sp.set_datetime(d)

        with self.assertRaises(ValueError):
            sp.set_position(-91, lng, h)
            sp.power()

        with self.assertRaises(ValueError):
            sp.set_position(90.1, lng, h)
            sp.power()

        with self.assertRaises(TypeError):
            sp.set_position('1', lng, h)
            sp.power()

        with self.assertRaises(ValueError):
            sp.set_position(lat, lng, -1)
            sp.power()

        with self.assertRaises(TypeError):
            sp.set_position(lat, lng, '21')
            sp.power()

    def test_orientation_exception(self):
        s, eff = 0, 0
        lat, lng, h = 0, 0, 0
        d = datetime(2019, 1, 21, 12, 0)

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_datetime(d)

        with self.assertRaises(ValueError):
            sp.set_orientation(1)
            sp.power()

    def test_datetime_exception(self):
        s, eff = 0, 0
        lat, lng, h = 0, 0, 0
        v = array([0, 0, -1])

        sp = solar_panel(s, eff)
        sp.set_position(lat, lng, h)
        sp.set_orientation(v)

        with self.assertRaises(TypeError):
            sp.set_datetime('a')
            sp.power()
