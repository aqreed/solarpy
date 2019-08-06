# coding: utf-8

"""
    Solar radiation model, based on Duffie, J.A., and
    Beckman, W. A., 1974, "Solar energy thermal processes"
"""

import numpy as np
from numpy import sin, cos, tan, deg2rad, rad2deg
from datetime import datetime, timedelta
from solarpy.utils import *


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
    check_day_range(n)

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


def solar_time(n, hour, minute, lng):
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
    lng : float
        longitude, east-west position wrt the Prime Meridian in degrees

    Returns
    -------
    solar time : tuple-like
        local solar time (hour, minute, second)
    """
    check_day_range(n)
    check_long_range(lng)

    if (hour < 0) or (hour > 23):
        raise ValueError('hour should be 0 <= hour <= 23')
    if (minute < 0) or (minute > 59):
        raise ValueError('minute should be 0 <= minute <= 59')

    # standard time
    t_std = datetime(datetime.now().year, 1, 1) + \
            timedelta(days=(n-1),
                      hours=hour,
                      minutes=minute)

    # displacement from standard meridian for that longitude
    lng_std = round(lng / 15) * 15
    delta_std_meridian = timedelta(minutes=(4 * (lng_std - lng)))

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
    if (hour < 0) or (hour > 23):
        raise ValueError('hour should be 0 <= hour <= 23')
    if (minute < 0) or (minute > 59):
        raise ValueError('minute should be 0 <= minute <= 59')

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
    check_lat_range(lat)

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
        zenith angle of incidence in radians
    """
    check_lat_range(lat)

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
    check_lat_range(lat)

    dec = declination(n)
    lat = deg2rad(lat)
    cos_ws = (-1) * tan(lat) * tan(dec)

    if abs(cos_ws) > 1:
        raise NoSunsetNoSunrise
    else:
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
    try:
        ws = sunset_hour_angle(n, lat)  # degrees

        aux = (rad2deg(ws) / 15) * 60 * 60  # seconds
        minutes, seconds = divmod(aux, 60)
        hours, minutes = divmod(minutes, 60)

        st = datetime(datetime.now().year, 1, 1) + \
             timedelta(days=n, hours=(12+hours), minutes=minutes)
        return st

    except NoSunsetNoSunrise as e:
        print(e.msg)


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
    sunset_hour : datetime-like
        time at sunset
    """
    try:
        ws = sunrise_hour_angle(n, lat)  # degrees

        aux = (rad2deg(ws) / 15) * 60 * 60  # seconds
        minutes, seconds = divmod(aux, 60)
        hours, minutes = divmod(minutes, 60)

        st = datetime(datetime.now().year, 1, 1) + \
             timedelta(days=n, hours=(12+hours), minutes=minutes)
        return st

    except NoSunsetNoSunrise as e:
        print(e.msg)


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
    check_lat_range(lat)

    dec = declination(n)
    lat = deg2rad(lat)

    tmp = -tan(lat) * tan(dec)

    # used mask to allow posterior visualization
    b = np.zeros(tmp.shape)

    b[(tmp < -1.0)] = 1
    b[(abs(tmp) < 1.0)] = (2 * np.arccos(tmp[(abs(tmp) < 1.0)]) / (2 * np.pi))
    b[(tmp > 1.0)] = 0

    return b * 24


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
    solar_az = solar_azimuth(n, lat, hour, minute)
    solar_alt = solar_altitude(n, lat, hour, minute)

    w = hour_angle(hour, minute)
    lh = daylight_hours(n, lat)

    try:
        w_sr = sunrise_hour_angle(n, lat)
        w_ss = sunset_hour_angle(n, lat)

        if (w > w_ss) or (w < w_sr):
            # the point on the earth surface is at night
            return np.array([0, 0, 0])
        else:
            return np.array([-cos(solar_az) * cos(solar_alt),
                             -sin(solar_az) * cos(solar_alt),
                             -sin(solar_alt)])

    except NoSunsetNoSunrise:
        if (lh == 0):
            # the point on the earth surface is in permanent darkness
            return np.array([0, 0, 0])
        else:
            # the point on the earth surface is in permanent light
            return np.array([-cos(solar_az) * cos(solar_alt),
                             -sin(solar_az) * cos(solar_alt),
                             -sin(solar_alt)])


def air_mass_KastenYoung1989(theta_z, h):
    """
    Returns the ratio obetween air mass crossed by a sun beam to the mass
    it would pass if the sun were in the zenith.

    Parameters
    ----------
    theta_z : float
        zenith angle of incidence in degrees
    h : float
        altitude above sea level in meters

    Returns
    -------
    m : float
        vector expressed in ECEF coordinates

    Notes
    -----
    Kasten, F.H., Young, A.T. (1989) "Revised optical air mass tables and
    approximation formula"
    """
    check_alt_range(h)

    # this conditional is needed to avoid KY1989 model limitations beyond 90º
    if theta_z < 91.5:
        theta_z_rad = deg2rad(theta_z)
        m = np.exp(-0.0001184 * h) / (cos(theta_z_rad) +
                                  0.50572 * (96.07995 - theta_z) ** (-1.634))
    else:
        theta_z_rad = deg2rad(91.5)
        m = np.exp(-0.0001184 * h) / (cos(theta_z_rad) +
                                  0.50572 * (96.07995 - 91.5) ** (-1.634))

    return m


def beam_irradiance(h, n, lat, hour, minute):
    """
    Returns the solar beam irradiance on a plane normal to the sun vector (not
    taking into account the diffuse component) at a certain altitude, day,
    latitude and time of the day (hour and minute).

    Parameters
    ----------
    h : float
        altitude above sea level in meters
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
    G : float
        beam irradiance in W/m2

    Notes
    -----
    Aglietti, G.S., Redi, S., Tatnall,A.R., Markvart, T., (2009) "Harnessing
    High-Altitude Solar Power"
    """
    alpha_int = 0.32  # atmospheric extinction. TODO: improve, as it changes
                      # throughout the year. Visible light? 4000-7000A
    prel = pressure(h) / pressure(0)  # pressure relation

    # the maximum zenith angle is the one that points to the horizon
    a = 6378137  # [m] Earth equatorial axis
    theta_lim = (1 / 2) * np.pi + np.arccos(a / (a + h))  # radians

    theta_zenith = theta_z(n, lat, hour, minute)  # radians

    if theta_zenith < theta_lim:
        m = air_mass_KastenYoung1989(rad2deg(theta_zenith), h)
        G = Gon(n) * np.exp(-prel * m * alpha_int)
    else:
        G = 0

    return G


def irradiance_on_plane(vnorm, h, n, lat, hour, minute):
    """
    Returns the solar beam irradiance on a plane (not taking into account the
    diffuse component) defined by its unit normal vector in NED frame at a
    certain altitude, day, latitude and time of the day (hour and minute).

    Parameters
    ----------
    vnorm : array-like
        unit vector normal to plane
    h : float
        altitude above sea level in meters
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
    G : float
        beam irradiance in W/m2
    """
    vsol = solar_vector_NED(n, lat, hour, minute)

    if (vsol == np.array([0, 0, 0])).all():
        # in case there is no sun (night or permanent darkness)
        G = 0
        theta = np.nan
    else:
        vnorm_abs = np.linalg.norm(vnorm)
        vsol_abs = np.linalg.norm(vsol)

        theta = np.arccos(np.dot(vnorm, vsol) / (vnorm_abs * vsol_abs))

        # for future solar panel applications: only one side of it has cells
        if cos(theta) > 0:
            G = beam_irradiance(h, n, lat, hour, minute) * cos(theta)
        else:
            G = 0

    return G
