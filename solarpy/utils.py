# coding: utf-8

"""
    Utilities for solar radiation calculations, based on Duffie,
    J.A., and Beckman, W. A., 1974, "Solar energy thermal processes"
"""

import numpy as np
from numpy import sin, cos, deg2rad, array
from datetime import datetime


def check_lat(lat):
    """
    Checks whether the input latitude is within range and correct type

    Parameters
    ----------
    lat : float or int
        latitude (-90 to 90) in degrees

    Returns
    -------
    None. Raises an exception in case
    """
    if isinstance(lat, (int, float)):
        if abs(lat) > 90:
            raise ValueError('latitude should be -90 <= latitude <= 90')
    else:
        raise TypeError('latitude should be "float" or "int"')

    return None


def check_long(lng):
    """
    Checks whether the input longitude is within range and correct type

    Parameters
    ----------
    lng : float or int
        longitude (-179 to 180) in degrees

    Returns
    -------
    None. Raises an exception in case
    """
    if isinstance(lng, (int, float)):
        if abs(lng) > 180:
            raise ValueError('longitude should be -180 <= longitude <= 180')
    else:
        raise TypeError('longitude should be "float" or "int"')

    return None


def check_alt(h):
    """
    Checks whether the input altitude is within range and correct type

    Parameters
    ----------
    h : float or int
        altitude (0 to 24k) in meters

    Returns
    -------
    None. Raises an exception in case
    """
    if isinstance(h, (int, float)):
        if ((h < 0) or (h > 24000)):
            raise ValueError('pressure model is only valid if 0 <= h <= 24000')
    else:
        raise TypeError('altitude should be "float" or "int"')

    return None


def day_of_the_year(date):
    """
    Returns the day of the year

    Parameters
    ----------
    date : datetime object
        date of interest

    Returns
    -------
    day : int
        day of the year (1 to 365)
    """
    if isinstance(date, datetime):
        return date.timetuple().tm_yday
    else:
        msg = "date must be a datetime object or array of datetime objects"
        raise TypeError(msg)


class NoSunsetNoSunrise(Exception):
    """
    Raised when the latitude in question is in permanent night or day
    """
    def __init__(self):
        self.msg = "Permanent night (or day) on this latitude on this day"


def lla2ecef(lat, lng, h):
    """
    Calculates geocentric coordinates (ECEF - Earth Centered, Earth Fixed) for
    a given set of latitude, longitude and altitude inputs.

    Parameters
    ----------
    lat : float
        latitude in degrees
    lng : float
        longitude in degrees
    h : float
        geometric altitude above sea level in meters

    Returns
    -------
    array-like
        ECEF coordinates in meters
    """
    check_lat(lat)
    check_long(lng)
    check_alt(h)

    a = 6378137  # [m] Earth equatorial axis
    b = 6356752.3142  # [m] Earth polar axis
    e = 0.081819190842622  # Earth eccentricity

    lat = deg2rad(lat)  # degrees to radians
    lng = deg2rad(lng)  # degrees to radians

    N = a / (1 - (e * sin(lat))**2)**(.5)

    x = (N + h) * cos(lat) * cos(lng)
    y = (N + h) * cos(lat) * sin(lng)
    z = (((b/a)**2) * N + h) * sin(lat)

    return array([x, y, z])


def ned2ecef(v_ned, lat, lng):
    """
    Converts vector from local geodetic horizon reference frame (NED - North,
    East, Down) at a given latitude and longitude to geocentric coordinates
    (ECEF - Earth Centered, Earth Fixed).

    Parameters
    ----------
    v_ned: array-like
        vector expressed in NED coordinates
    lat : float
        latitude in degrees
    lng : float
        longitude in degrees

    Returns
    -------
    v_ecef : array-like
        vector expressed in ECEF coordinates
    """
    check_lat(lat)
    check_long(lng)

    lat = deg2rad(lat)
    lng = deg2rad(lng)

    Lne = array([[-sin(lat) * cos(lng), -sin(lat) * sin(lng), cos(lat)],
                 [-sin(lng), cos(lng), 0],
                 [-cos(lat) * cos(lng), -cos(lat) * sin(lng), -sin(lat)]])

    Len = Lne.transpose()
    v_ecef = Len.dot(v_ned)

    return v_ecef


def pressure(h):
    """
    Interim function that returns ISA standard day pressure at a desired
    altitud. A new release of https://github.com/AeroPython/scikit-aero
    that includes the COESA atmospheric model (which would allow its
    installation with pip) will substitute this function.

    Parameters
    ----------
    h : float
        altitude above sea level in meters

    Returns
    -------
    p : float
        pressure in Pa

    Notes
    -----
    http://www.pdas.com/atmosTable2SI.html until 20km
    http://www.pdas.com/atmosTable1SI.html until 24km
    """
    check_alt(h)

    alt_ = np.append(np.linspace(0, 20e3, 21),
                     np.linspace(22e3, 24e3, 2))

    p_ = array([101325, 89876, 79501, 70121, 61660, 54048, 47217, 41105,
                35651, 30800, 26499, 22699, 19399, 16579, 14170, 12111,
                10352, 8849, 7565, 6467, 5529, 4047, 2972])

    return np.interp(h, alt_, p_)
