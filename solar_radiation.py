# coding: utf-8

"""
    Solar radiation model, based on Duffie, J.A., and
    Beckman, W. A., 1974, "Solar energy thermal processes"
"""

import numpy as np
from numpy import sin, cos, tan, deg2rad, rad2deg
from datetime import datetime, timedelta


def day_of_the_year(month, day):
    """
    Returns the day of the year

    Parameters
    ----------
    month : integer
        day of the year (1 to 12)
    day : integer
        hour of the day (0 to 31)
    Returns
    -------
    day : integer
        day of the year(1 to 365)
    """
    t = datetime(datetime.now().year, month, day) - \
        datetime(datetime.now().year, 1, 1)

    return t.days


def B_nth_day(n):
    """
    Day of the year angle

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    Returns
    -------
    B : float
        angle of the day of the year in radians
    """
    if isinstance(n, np.ndarray) and ((n < 1).any() or (n > 365).any()):
            raise ValueError('n should be 1 <= n <= 365')
    elif isinstance(n, int) and ((n < 1) or (n > 365)):
            raise ValueError('n should be 1 <= n <= 365')

    return deg2rad((n - 1) * (360 / 365))


def Gon(n):
    """
    Extraterrestrial radiation on a plane normal to
    the radiation on the nth day of the year.

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)

    Returns
    -------
    Gon : float
        extraterrestrial radiation in W/m2
    """
    B = B_nth_day(n)

    return 1367 * (1.00011 + 0.034221 * cos(B) +
                   0.00128 * sin(B) + 0.000719 * cos(2 * B) +
                   0.000077 * sin(2 * B))


def Eq_time(n):
    """
    Equation of time on the nth day of the year.

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)

    Returns
    -------
    E : float
        equation of time in minutes
    """
    B = B_nth_day(n)

    return 229.2 * (0.000075 + 0.001868 * cos(B) -
                    0.032077 * sin(B) - 0.014615 * cos(2 * B) -
                    0.04089 * sin(2 * B))


def declination(n):
    """
    Angular position of the Sun at solar noon. Must comply
    with -23.45º < declination < 23.45º

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)

    Returns
    -------
    declination : float
        declination in radians
    """
    B = B_nth_day(n)

    return 0.006918 - 0.399912 * cos(B) + 0.070257 * sin(B) - \
           0.006758 * cos(2 * B) + 0.000907 * sin(2 * B) - \
           0.002679 * cos(3 * B) + 0.00148 * sin(3 * B)


def solar_time(n, hour, minute, long):
    """
    Solar time (local) for a particular day of the year (nth),
    hour-minute and longitude

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    hour : integer
        hour of the day (0 to 23)
    minute : integer
        minutes (0 to 59)
    longitude : float
        east-west position wrt the Prime Meridian in degrees
    Returns
    -------
    solar time : tuple-like
        local solar time (hour, minute, second)
    """
    if isinstance(n, np.ndarray) and ((n < 1).any() or (n > 365).any()):
            raise ValueError('n should be 1 <= n <= 365')
    elif isinstance(n, int) and ((n < 1) or (n > 365)):
            raise ValueError('n should be 1 <= n <= 365')

    # standard time
    t_std = datetime(datetime.now().year, 1, 1) + \
            timedelta(days=(n-1),
                      hours=hour,
                      minutes=minute)

    # displacement from standard meridian for that longitude
    long_std = round(long / 15) * 15
    delta_std_meridian = timedelta(minutes=(4 * (long_std - long)))

    # eq. of time for that day
    E = timedelta(minutes=Eq_time(n))
    t_solar = t_std + delta_std_meridian + E

    return t_solar.hour, t_solar.minute, t_solar.second


def hour_angle(hour, minute):
    """
    Angular displacement of the sun east-west of the local
    meridian. 15 degrees per hour, morning < 0 < afternoon

    Parameters
    ----------
    hour : integer
        hour (solar time) of the day (0 to 23)
    minute : integer
        solar (solar time) minutes (0 to 59)
    Returns
    -------
    hour angle : float
        local hour angle in radians
    """
    w = ((hour + (minute / 60)) - 12) * 15

    return deg2rad(w)


def theta(n, lat, beta, surf_az, hour, minute):
    """
    Angle of incidence of the sun beam on a surface wrt the normal
    to that surface, for a particular day of the year (nth), latitude,
    surface slope, surface azimuth and hour-minute.

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    beta : float
        slope angle of the surface wrt the local horizon
        in degrees (0 to 180)
    surf_az : float
        azimuth angle of the surface in degrees wrt the local
        meridian (-180 to 180). 0-> south, east negative
    hour : integer
        hour of the day (0 to 23)
    minute : integer
        minutes (0 to 59)
    Returns
    -------
    theta : float
            angle of incidence in radians
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    dec = declination(n)
    lat = deg2rad(lat)
    beta = deg2rad(beta)
    surf_az = deg2rad(surf_az)
    w = hour_angle(hour, minute)

    cos_theta = sin(dec) * sin(lat) * cos(beta) - \
                sin(dec) * cos(lat) * sin(beta) * cos(surf_az) + \
                cos(dec) * cos(lat) * cos(beta) *                cos(w) + \
                cos(dec) * sin(lat) * sin(beta) * cos(surf_az) * cos(w) + \
                cos(dec) *            sin(beta) * sin(surf_az) * sin(w)

    return np.arccos(cos_theta)


def theta_z(n, lat, hour, minute):
    """
    * Zenith angle *

    Angle of incidence of the sun beam on a horizontal surface wrt
    the normal to that surface, for a particular day of the year (nth),
    latitude and hour-minute.

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    hour : integer
        hour of the day (0 to 23)
    minute : integer
        minutes (0 to 59)
    Returns
    -------
    theta_z : float
        angle of incidence in radians
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    dec = declination(n)
    lat = deg2rad(lat)
    w = hour_angle(hour, minute)

    cos_theta_z = sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(w)

    return np.arccos(cos_theta_z)


def solar_azimuth(n, lat, hour, minute):
    """
    * Solar azimuth angle *

    Angle between the projection of the sun beam on a horizontal
    surface wrt N-S, for a particular day of the year (nth),
    latitude and hour-minute. Positive to the West.

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    hour : integer
        hour of the day (0 to 23)
    minute : integer
        minutes (0 to 59)
    Returns
    -------
    solar_az : float
        azimuth angle in radians
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    # to avoid undefined values at lat = 90º or lat = -90º
    # the error incurred is acceptable
    if abs(lat) == 90:
        lat = np.sign(lat) * 89.999

    w = hour_angle(hour, minute)
    dec = declination(n)
    th_z = theta_z(n, lat, hour, minute)
    lat = deg2rad(lat)

    tmp = (cos(th_z) * sin(lat) - sin(dec)) / (sin(th_z) * cos(lat))

    # herculean fight against floating-point errors
    if (abs(tmp) > 1):
        tmp = int(tmp)  # TODO: improve

    # to avoid undefined values at noon (12:00)
    if w == 0:
        s = 1
    else:
        s = np.sign(w)

    return s * np.arccos(tmp)


def solar_altitude(n, lat, hour, minute):
    """
    * Solar altitude angle *

    Angle between the projection of the sun beam on a horizontal
    surface wrt the beam, for a particular day of the year (nth),
    latitude and hour-minute.

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    hour : integer
        hour of the day (0 to 23)
    minute : integer
        minutes (0 to 59)
    Returns
    -------
    solar_altitude : float
        altitude angle in radians
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    th_z = theta_z(n, lat, hour, minute)

    return np.arcsin(cos(th_z))


def sunset_hour_angle(n, lat):
    """
    When theta_z = 90º

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    Returns
    -------
    sunset_hour_angle : float
        hour angle at sunset in radians
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    dec = declination(n)
    lat = deg2rad(lat)
    cos_ws = (-1) * tan(lat) * tan(dec)

    return np.arccos(cos_ws)


def sunset_time(n, lat):
    """
    Calculates the time (hours, minutes) at sunset

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    Returns
    -------
    sunset_hour : datetime-like
        time at sunset
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    ws = sunset_hour_angle(n, lat)  # degrees

    aux = (rad2deg(ws) / 15) * 60 * 60  # seconds
    minutes, seconds = divmod(aux, 60)
    hours, minutes = divmod(minutes, 60)

    st = datetime(datetime.now().year, 1, 1) + \
         timedelta(days=n, hours=(12+hours), minutes=minutes)

    return st


def sunrise_hour_angle(n, lat):
    """
    When theta_z = -90º

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    Returns
    -------
    sunrise_hour_angle : float
        hour angle at sunrise in radians
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    return -sunset_hour_angle(n, lat)


def sunrise_time(n, lat):
    """
    Calculates the time (hours, minutes) at sunrise

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    Returns
    -------
    sunset_hour : datetime-like
        time at sunset
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    ws = sunrise_hour_angle(n, lat)  # degrees

    aux = (rad2deg(ws) / 15) * 60 * 60  # seconds
    minutes, seconds = divmod(aux, 60)
    hours, minutes = divmod(minutes, 60)

    st = datetime(datetime.now().year, 1, 1) + \
         timedelta(days=n, hours=(12+hours), minutes=minutes)

    return st


def daylight_hours(n, lat):
    """
    Nº of hours of light for a particular day

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    Returns
    -------
    day_hours : float
        number of hours of light within the day
    Note
    ----
        http://mathforum.org/library/drmath/view/56478.html
    """
    dec = declination(n)
    lat = deg2rad(lat)

    tmp = -tan(lat) * tan(dec)

    # used mask to allow posterior visualization
    b = np.zeros(tmp.shape)

    b[(tmp < -1.0)] = 1
    b[(abs(tmp) < 1.0)] = (2 * np.arccos(tmp[(abs(tmp) < 1.0)]) / (2 * np.pi))
    b[(tmp > 1.0)] = 0

    return b * 24


def lla2ecef(lat, long, h):
    """
    Calculates geocentric coordinates (ECEF - Earth Centered, Earth Fixed) for
    a given set of latitude, longitude and altitude inputs.

    Parameters
    ----------
    lat : float
        latitude in degrees
    long : float
        longitude in degrees
    h : float
        altitude above sea level in feet

    Returns
    -------
    array-like
        ECEF coordinates in meters
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    a = 6378137  # [m] Earth equatorial axis
    b = 6356752.3142  # [m] Earth polar axis
    e = 0.081819190842622  # Earth eccentricity

    lat = deg2rad(lat)  # degrees to radians
    long = deg2rad(long)  # degrees to radians
    h = h * 0.3048  # feets to meters

    N = a / (1 - (e * sin(lat))**2)**(.5)

    x = (N + h) * cos(lat) * cos(long)
    y = (N + h) * cos(lat) * sin(long)
    z = (((b/a)**2) * N + h) * sin(lat)

    return np.array([x, y, z])


def ned2ecef(v_ned, lat, long):
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
    long : float
        longitude in degrees

    Returns
    -------
    v_ecef : array-like
        vector expressed in ECEF coordinates
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    lat = deg2rad(lat)
    long = deg2rad(long)

    Lne = np.array([[-sin(lat) * cos(long), -sin(lat) * sin(long), cos(lat)],
                    [-sin(long), cos(long), 0],
                    [-cos(lat) * cos(long), -cos(lat) * sin(long), -sin(lat)]])

    Len = Lne.transpose()
    v_ecef = Len.dot(v_ned)

    return v_ecef


def solar_vector_NED(n, lat, hour, minute):
    """
    Calculates solar vector (sun beam) in local geodetic horizon reference
    frame (NED - North, East, Down) of a point on the Earth surface at a
    defined time (day, hour, minute) at a defined latitude.

    Parameters
    ----------
    n : integer
        day of the year (1 to 365)
    lat : float
        latitude (-90 to 90) in degrees
    hour : integer
        hour of the day (0 to 23)
    minute : integer
        minutes (0 to 59)
    Returns
    -------
    array-like
        vector of the solar beam
    """
    if abs(lat) > 90:
        raise ValueError('latitude should be -90 < lat < 90')

    solar_az = solar_azimuth(n, lat, hour, minute)
    solar_alt = solar_altitude(n, lat, hour, minute)

    return np.array([-cos(solar_az) * cos(solar_alt),
                     -sin(solar_az) * cos(solar_alt),
                     -sin(solar_alt)])
