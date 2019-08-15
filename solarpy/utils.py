# coding: utf-8

"""
    Utilities for solar radiation calculations, based on Duffie,
    J.A., and Beckman, W. A., 1974, "Solar energy thermal processes"
"""

import numpy as np
from numpy import sin, cos, tan, deg2rad, array, ndarray
from datetime import datetime, timedelta


def check_lat_range(lat):
    """
    Checks whether the input latitude is within range

    Parameters
    ----------
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    None. Raises an exception in case
    """
    if isinstance(lat, ndarray) and ((abs(lat) > 90).any()):
            raise ValueError('latitude should be -90 <= latitude <= 90')
    elif isinstance(lat, int) and (abs(lat) > 90):
            raise ValueError('latitude should be -90 <= latitude <= 90')
    elif isinstance(lat, float) and (abs(lat) > 90):
            raise ValueError('latitude should be -90 <= latitude <= 90')

    return None


def check_long_range(lng):
    """
    Checks whether the input longitude is within range

    Parameters
    ----------
    lng : float
        longitude (-179 to 180) in degrees

    Returns
    -------
    None. Raises an exception in case
    """
    if isinstance(lng, ndarray) and ((lng < -180).any() or (lng > 180).any()):
            raise ValueError('longitude should be -180 <= longitude <= 180')
    elif isinstance(lng, int) and ((lng < -180) or (lng > 180)):
            raise ValueError('longitude should be -180 <= longitude <= 180')
    elif isinstance(lng, float) and ((lng < -180) or (lng > 180)):
            raise ValueError('longitude should be -180 <= longitude <= 180')

    return None


def check_alt_range(h):
    """
    Checks whether the input altitude is within range

    Parameters
    ----------
    h : float
        altitude (0 to 24k) in meters

    Returns
    -------
    None. Raises an exception in case
    """
    if isinstance(h, ndarray) and ((h < 0).any() or (h > 24000).any()):
            raise ValueError('pressure model is only valid if 0 <= h <= 24000')
    elif isinstance(h, int) and ((h < 0) or (h > 24000)):
            raise ValueError('pressure model is only valid if 0 <= h <= 24000')
    elif isinstance(h, float) and ((h < 0) or (h > 24000)):
            raise ValueError('pressure model is only valid if 0 <= h <= 24000')

    return None


def day_of_the_year(date):
    """
    Returns the day of the year

    Parameters
    ----------
    date : datetime object or array-like (datetime objects inside)
        date of interest

    Returns
    -------
    day : int or array-like (int inside)
        day of the year (1 to 365)
    """
    if (isinstance(date, ndarray) and all(isinstance(i, datetime)
        for i in date)):
        # the parameter is an array of datetime objects
        return array([i.timetuple().tm_yday for i in date])

    elif isinstance(date, datetime):
        # the parameter is a datetime object
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
    check_lat_range(lat)
    check_long_range(lng)
    check_alt_range(h)

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
    check_lat_range(lat)
    check_long_range(lng)

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
    check_alt_range(h)

    alt_ = np.append(np.linspace(0, 20e3, 21),
                     np.linspace(22e3, 24e3, 2))

    p_ = array([101325, 89876, 79501, 70121, 61660, 54048, 47217, 41105,
                35651, 30800, 26499, 22699, 19399, 16579, 14170, 12111,
                10352, 8849, 7565, 6467, 5529, 4047, 2972])

    return np.interp(h, alt_, p_)
