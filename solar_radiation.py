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
    B = deg2rad((n - 1) * (360 / 365))

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
    B = deg2rad((n - 1) * (360 / 365))

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
    B = deg2rad((n - 1) * (360 / 365))

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
    beta = 0
    surf_az = 0

    return theta(n, lat, beta, surf_az, hour, minute)


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
    beta = 0
    surf_az = 0

    w = hour_angle(hour, minute)
    dec = declination(n)
    th_z = theta_z(n, lat, hour, minute)
    lat = deg2rad(lat)

    tmp = (cos(th_z) * sin(lat) - sin(dec)) / (sin(th_z) * cos(lat))

    return np.sign(w) * np.arccos(tmp)


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
    sunset_hour : tuple-like
        time at sunset (hours, minutes)
    """
    ws = sunset_hour_angle(n, lat)  # degrees

    aux = (rad2deg(ws) / 15) * 60 * 60  # seconds
    minutes, seconds = divmod(aux, 60)
    hours, minutes = divmod(minutes, 60)

    st = datetime(datetime.now().year, 1, 1) + \
         timedelta(days=n, hours=(12+hours), minutes=minutes)

    return (st.hour, st.minute)


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
    sunset_hour : tuple-like
        time at sunrise (hours, minutes)
    """
    ws = sunrise_hour_angle(n, lat)  # degrees

    aux = (rad2deg(ws) / 15) * 60 * 60  # seconds
    minutes, seconds = divmod(aux, 60)
    hours, minutes = divmod(minutes, 60)

    st = datetime(datetime.now().year, 1, 1) + \
         timedelta(days=n, hours=(12+hours), minutes=minutes)

    return (st.hour, st.minute)


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
