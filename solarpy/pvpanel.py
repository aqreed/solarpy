# coding: utf-8

"""
    Photovoltaic panel class
"""
from .radiation import irradiance_on_plane


class solar_panel(object):
    """
    Photovoltaic solar panel class

    Parameters
    ----------
    s : float
        panel surface in m2
    eff : float
        panel efficiency
    """
    def __init__(self, s, eff, id_name=None):
        if (s >= 0):
            self.s = s
        else:
            raise ValueError('surface must be s >= 0')

        if (eff >= 0) and (eff <= 1):
            self.eff = eff
        else:
            raise ValueError('efficiency must be 0 <= eff <= 1')

        if isinstance(id_name, str):
            self.id_name = id_name
        elif id_name is None:
            pass
        else:
            raise TypeError('the id name must be a string')

    def set_position(self, lat, lng, h):
        """
        Sets LLA position (latitude, longitude, altitude)

        Parameters
        ----------
        lat : float
            latitude (-90 to 90) in degrees
        lng : float
            longitude, (-180 to 180) in degrees
        h : float
            altitude above sea level in meters
        """
        self.lat = lat
        self.lng = lng
        self.h = h

    def set_orientation(self, vnorm):
        """
        Sets the orientation of the panel in NED frame

        Parameters
        ----------
        vnorm : array-like
            unit vector normal to plane in NED frame
        """
        self.vnorm = vnorm

    def set_datetime(self, date):
        """
        Sets date and time

        Parameters
        ----------
        date : datetime object
            date and *solar* time
        """
        self.date = date

    def power(self):
        """
        Returns the output power of a solar panel
        """
        return irradiance_on_plane(self.vnorm, self.h,
                                   self.date, self.lat) * self.s * self.eff
