# coding: utf-8

"""
    Solar radiation model, based on Duffie, J.A., and
    Beckman, W. A., 1974, "Solar energy thermal processes"
"""

import numpy as np
from numpy import sin, cos, tan, deg2rad, rad2deg,\
                  array, arccos, exp, ndarray
from datetime import datetime, timedelta
from solarpy.utils import *


def B_nth_day(date):
    """
    Day-of-the-year angle on a desired date and time.

    Parameters
    ----------
    date : datetime object or array-like (datetime objects inside)
        date of interest

    Returns
    -------
    B : float or array-like (float inside)
        angle of the day of the year in radians
    """
    try:
        n = day_of_the_year(date)

        if (isinstance(date, ndarray) and
            all(isinstance(i, datetime) for i in date)):
            # the parameter is an array of datetime objects
            return array([deg2rad((i - 1) * (360 / 365)) for i in n])

        elif isinstance(date, datetime):
            # the parameter is a datetime object
            return deg2rad((n - 1) * (360 / 365))

    except TypeError as e:
        raise e


def Gon(date):
    """
    Extraterrestrial radiation on a plane normal to
    the radiation on a desired date and time.

    Parameters
    ----------
    date : datetime object or array-like (datetime objects inside)
        date of interest

    Returns
    -------
    Gon : float or array-like (float inside)
        extraterrestrial radiation in W/m2
    """
    try:
        B = B_nth_day(date)

        if (isinstance(date, ndarray) and
            all(isinstance(i, datetime) for i in date)):
            # the parameter is an array of datetime objects
            return array([1367 * (1.00011 + 0.034221 * cos(i) +
                                  0.00128 * sin(i) + 0.000719 * cos(2 * i) +
                                  0.000077 * sin(2 * i)) for i in B])

        elif isinstance(date, datetime):
            # the parameter is a datetime object
            return 1367 * (1.00011 + 0.034221 * cos(B) +
                           0.00128 * sin(B) + 0.000719 * cos(2 * B) +
                           0.000077 * sin(2 * B))

    except TypeError as e:
        raise e


def Eq_time(date):
    """
    Equation of time on a desired date and time.

    Parameters
    ----------
    date : datetime object or array-like (datetime objects inside)
        date of interest

    Returns
    -------
    E : float or array-like (float inside)
        equation of time in minutes
    """
    try:
        B = B_nth_day(date)

        if (isinstance(date, ndarray) and
            all(isinstance(i, datetime) for i in date)):
            # the parameter is an array of datetime objects
            return array([229.2 * (0.000075 + 0.001868 * cos(i) -
                                   0.032077 * sin(i) -
                                   0.014615 * cos(2 * i) -
                                   0.04089 * sin(2 * i)) for i in B])

        elif isinstance(date, datetime):
            # the parameter is a datetime object
            return 229.2 * (0.000075 + 0.001868 * cos(B) -
                            0.032077 * sin(B) - 0.014615 * cos(2 * B) -
                            0.04089 * sin(2 * B))

    except TypeError as e:
        raise e


def declination(date):
    """
    Angular position of the Sun at solar noon on a desired date and time.
    Must comply with -23.45º < declination < 23.45º

    Parameters
    ----------
    date : datetime object or array-like (datetime objects inside)
        date of interest

    Returns
    -------
    declination : float or array-like (float inside)
        declination in radians
    """

    try:
        B = B_nth_day(date)

        if (isinstance(date, ndarray) and
            all(isinstance(i, datetime) for i in date)):
            # the parameter is an array of datetime objects
            return array([0.006918 - 0.399912 * cos(i) + 0.070257 * sin(i) -
                          0.006758 * cos(2 * i) + 0.000907 * sin(2 * i) -
                          0.002679 * cos(3 * i) + 0.00148 * sin(3 * i)
                          for i in B])

        elif isinstance(date, datetime):
            # the parameter is a datetime object
            return 0.006918 - 0.399912 * cos(B) + 0.070257 * sin(B) - \
                   0.006758 * cos(2 * B) + 0.000907 * sin(2 * B) - \
                   0.002679 * cos(3 * B) + 0.00148 * sin(3 * B)

    except TypeError as e:
        raise e


def standard2solar_time(date, lng):
    """
    Solar time for a particular longitude, date and *standard* time.

    Parameters
    ----------
    date : datetime object
        standard (or local) time
    lng : float
        longitude, east-west position wrt the Prime Meridian in degrees

    Returns
    -------
    solar time : datetime object
        solar time
    """
    check_long_range(lng)

    if isinstance(date, datetime):
        # standard time
        t_std = date

        # displacement from standard meridian for that longitude
        lng_std = round(lng / 15) * 15
        delta_std_meridian = timedelta(minutes=(4 * (lng_std - lng)))

        # eq. of time for that day
        E = timedelta(minutes=Eq_time(date))
        t_solar = t_std + delta_std_meridian + E

        return t_solar
    else:
        raise TypeError('date must be a datetime object')


def hour_angle(date):
    """
    Angular displacement of the sun east-west of the local meridian for a
    date, *solar* time and latitude.
    Note: 15 degrees per hour, morning < 0 < afternoon

    Parameters
    ----------
    date : datetime object
        date and *solar* time

    Returns
    -------
    hour angle : float
        local hour angle in radians
    """
    if isinstance(date, datetime):
        w = (date.hour + (date.minute / 60) - 12) * 15
        return deg2rad(w)
    else:
        raise TypeError('date must be a datetime object')


def theta(date, lat, beta, surf_az):
    """
    Angle of incidence of the sun beam on a surface wrt the normal
    to that surface, for a date, *solar* time, latitude, surface
    slope and surface azimuth.

    Parameters
    ----------
    date : datetime object
        date and *solar* time
    lat : float
        latitude (-90 to 90) in degrees
    beta : float
        slope angle of the surface wrt the local horizon
        in degrees (0 to 180)
    surf_az : float
        azimuth angle of the surface in degrees wrt the local
        meridian (-180 to 180). 0-> south, east negative

    Returns
    -------
    theta : float
            angle of incidence in radians
    """
    check_lat_range(lat)

    if isinstance(date, datetime):
        dec = declination(date)
        lat = deg2rad(lat)
        beta = deg2rad(beta)
        surf_az = deg2rad(surf_az)
        w = hour_angle(date)

        cos_theta = sin(dec) * sin(lat) * cos(beta) - \
                    sin(dec) * cos(lat) * sin(beta) * cos(surf_az) + \
                    cos(dec) * cos(lat) * cos(beta) *                cos(w) + \
                    cos(dec) * sin(lat) * sin(beta) * cos(surf_az) * cos(w) + \
                    cos(dec) *            sin(beta) * sin(surf_az) * sin(w)
        return arccos(cos_theta)
    else:
        raise TypeError('date must be a datetime object')


def theta_z(date, lat):
    """
    * Zenith angle *

    Angle of incidence of the sun beam on a horizontal surface wrt the
    normal to that surface, for a date, *solar* time and latitude.

    Parameters
    ----------
    date : datetime object
        date and *solar* time
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    theta_z : float
        zenith angle of incidence in radians
    """
    check_lat_range(lat)

    if isinstance(date, datetime):
        dec = declination(date)
        lat = deg2rad(lat)
        w = hour_angle(date)

        cos_theta_z = sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(w)
        return arccos(cos_theta_z)
    else:
        raise TypeError('date must be a datetime object')


def solar_azimuth(date, lat):
    """
    * Solar azimuth angle *

    Angle between the projection of the sun beam on a horizontal surface
    wrt N-S, for a date, *solar* time and latitude. Positive to the West.

    Parameters
    ----------
    date : datetime object
        date and *solar* time
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    solar_az : float
        azimuth angle in radians
    """
    check_lat_range(lat)

    # to avoid undefined values at lat = 90º or lat = -90º
    # the error incurred is acceptable
    if abs(lat) == 90:
        lat = np.sign(lat) * 89.999

    if isinstance(date, datetime):
        w = hour_angle(date)
        dec = declination(date)
        th_z = theta_z(date, lat)
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

        return s * arccos(tmp)
    else:
        raise TypeError('date must be a datetime object')


def solar_altitude(date, lat):
    """
    * Solar altitude angle *

    Angle between the projection of the sun beam on a horizontal
    surface wrt the beam, for a date, *solar* time and latitude.

    Parameters
    ----------
    date : datetime object
        date and *solar* time
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    solar_altitude : float
        altitude angle in radians
    """
    if isinstance(date, datetime):
        th_z = theta_z(date, lat)
        return np.arcsin(cos(th_z))
    else:
        raise TypeError('date must be a datetime object')


def sunset_hour_angle(date, lat):
    """
    Sunset hour angle for a date and latitude

    Note: theta_z = 90º

    Parameters
    ----------
    date : datetime object
        date (indifferent time)
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    sunset_hour_angle : float
        hour angle at sunset in radians
    """
    check_lat_range(lat)

    if isinstance(date, datetime):
        dec = declination(date)
        lat = deg2rad(lat)
        cos_ws = (-1) * tan(lat) * tan(dec)

        if abs(cos_ws) > 1:
            raise NoSunsetNoSunrise
        else:
            return arccos(cos_ws)
    else:
        raise TypeError('date must be a datetime object')


def sunset_time(date, lat):
    """
    Calculates the *solar* time at sunset for a desired date

    Parameters
    ----------
    date : datetime object
        date (indifferent time)
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    sunset_hour : datetime object
        time at sunset
    """
    if isinstance(date, datetime):
        try:
            ws = sunset_hour_angle(date, lat)  # degrees

            aux = (rad2deg(ws) / 15) * 60 * 60  # seconds
            minutes, seconds = divmod(aux, 60)
            hours, minutes = divmod(minutes, 60)

            st = datetime(date.year, date.month, date.day) + \
                 timedelta(hours=(12+hours), minutes=minutes)
            return st

        except NoSunsetNoSunrise as e:
            raise e
    else:
        raise TypeError('date must be a datetime object')


def sunrise_hour_angle(date, lat):
    """
    Sunset hour angle for a date and latitude

    Note: theta_z = -90º

    Parameters
    ----------
    date : datetime object
        date (indifferent time)
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    sunrise_hour_angle : float
        hour angle at sunrise in radians
    """
    return -sunset_hour_angle(date, lat)


def sunrise_time(date, lat):
    """
    Calculates the *solar* time at sunrise for a desired date

    Parameters
    ----------
    date : datetime object
        date (indifferent time)
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    sunset_hour : datetime object
        time at sunrise
    """
    if isinstance(date, datetime):
        try:
            ws = sunrise_hour_angle(date, lat)  # degrees

            aux = (rad2deg(ws) / 15) * 60 * 60  # seconds
            minutes, seconds = divmod(aux, 60)
            hours, minutes = divmod(minutes, 60)

            st = datetime(date.year, date.month, date.day) + \
                 timedelta(hours=(12+hours), minutes=minutes)
            return st

        except NoSunsetNoSunrise as e:
            raise e
    else:
        raise TypeError('date must be a datetime object')


def daylight_hours(date, lat):
    """
    Nº of hours of light for a particular day

    Parameters
    ----------
    date : datetime object
        date (indifferent time)
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

    if isinstance(date, datetime):
        dec = declination(date)
        lat = deg2rad(lat)

        tmp = -tan(lat) * tan(dec)

        if (tmp < -1.0):
            return 24
        elif (tmp > 1.0):
            return 0
        else:
            return 24 * (2 * arccos(tmp) / (2 * np.pi))
    else:
        raise TypeError('date must be a datetime object')


def solar_vector_NED(date, lat):
    """
    Calculates solar vector (sun beam) in local geodetic horizon reference
    frame (NED - North, East, Down) of a point on the Earth surface at a
    defined date, time and latitude.

    Parameters
    ----------
    date : datetime object
        date and *solar* time
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    array-like
        vector of the solar beam
    """
    if isinstance(date, datetime):
        solar_az = solar_azimuth(date, lat)
        solar_alt = solar_altitude(date, lat)

        w = hour_angle(date)
        lh = daylight_hours(date, lat)

        try:
            w_sr = sunrise_hour_angle(date, lat)
            w_ss = sunset_hour_angle(date, lat)

            if (w > w_ss) or (w < w_sr):
                # the point on the earth surface is at night
                return array([0, 0, 0])
            else:
                return array([-cos(solar_az) * cos(solar_alt),
                              -sin(solar_az) * cos(solar_alt),
                              -sin(solar_alt)])

        except NoSunsetNoSunrise:
            if (lh == 0):
                # the point on the earth surface is in permanent darkness
                return array([0, 0, 0])
            else:
                # the point on the earth surface is in permanent light
                return array([-cos(solar_az) * cos(solar_alt),
                              -sin(solar_az) * cos(solar_alt),
                              -sin(solar_alt)])
    else:
        raise TypeError('date must be a datetime object')


def air_mass_KastenYoung1989(theta_z, h, limit=True):
    """
    Returns the ratio between air mass crossed by a sun beam to the mass
    it would pass if the sun were in the zenith at any altitude.

    Parameters
    ----------
    theta_z : float
        zenith angle of incidence in degrees
    h : float
        altitude above sea level in meters
    limit : boolean
        activates or deaactivates altitude limit

    Returns
    -------
    m : float
        ratio

    Notes
    -----
    Kasten, F.H., Young, A.T. (1989) "Revised optical air mass tables and
    approximation formula"
    """
    # needed until the atmosphere (pressure) model is extended beyond 24km
    if limit:
        check_alt_range(h)
    else:
        pass

    # this saturation is an interim solution needed to avoid KY1989 model
    # limitations beyond 90º. TODO: improve
    if theta_z < 91.5:
        theta_z_rad = deg2rad(theta_z)
        m = exp(-0.0001184 * h) / (cos(theta_z_rad) +
                                   0.50572 * (96.07995 - theta_z) ** (-1.634))
    else:
        theta_z_rad = deg2rad(91.5)
        m = exp(-0.0001184 * h) / (cos(theta_z_rad) +
                                   0.50572 * (96.07995 - 91.5) ** (-1.634))

    return m


def air_mass_Young1994(theta_z):
    """
    Returns the ratio between air mass crossed by a sun beam to the mass
    it would pass if the sun were in the zenith at sea level.

    Parameters
    ----------
    theta_z : float
        zenith angle of incidence in degrees

    Returns
    -------
    m : float
        ratio

    Notes
    -----
    Kasten, F.H. (1994) "Air mass and refraction"
    """
    th_z = deg2rad(theta_z)
    co_thz = cos(th_z)

    a = 1.002432 * co_thz**2 + 0.148386 * co_thz + 0.0096467
    b = co_thz**3 + 0.149864 * co_thz**2 + 0.0102963 * co_thz + 0.000303978
    m = a / b

    return m


def beam_irradiance(h, date, lat):
    """
    Returns the solar beam irradiance on a plane normal to the sun vector (not
    taking into account the diffuse component) at a certain altitude, date,
    time and latitude

    Parameters
    ----------
    h : float
        altitude above sea level in meters
    date : datetime object
        date and *solar* time
    lat : float
        latitude (-90 to 90) in degrees

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
    theta_lim = (1 / 2) * np.pi + arccos(a / (a + h))  # radians

    if isinstance(date, datetime):
        theta_zenith = theta_z(date, lat)  # radians

        if theta_zenith < theta_lim:
            m = air_mass_KastenYoung1989(rad2deg(theta_zenith), h)
            G = Gon(date) * exp(-prel * m * alpha_int)
        else:
            G = 0

        return G
    else:
        raise TypeError('date must be a datetime object')


def irradiance_on_plane(vnorm, h, date, lat):
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
    date : datetime object
        date and *solar* time
    lat : float
        latitude (-90 to 90) in degrees

    Returns
    -------
    G : float
        beam irradiance in W/m2
    """
    try:
        vsol = solar_vector_NED(date, lat)

        if isinstance(date, datetime):
            if (vsol == array([0, 0, 0])).all():
                # in case there is no sun (night or permanent darkness)
                G = 0
                theta = np.nan
            else:
                vnorm_abs = np.linalg.norm(vnorm)
                vsol_abs = np.linalg.norm(vsol)

                theta = arccos(np.dot(vnorm, vsol) / (vnorm_abs * vsol_abs))

                # for future solar panel applications: only one side has cells
                if cos(theta) > 0:
                    G = beam_irradiance(h, date, lat) * cos(theta)
                else:
                    G = 0

            return G
    except TypeError as e:
        raise e
