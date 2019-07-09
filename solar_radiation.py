# coding: utf-8

"""
	Solar radiation model, based on Duffie, J.A., and Beckman, W. A., 1974, 
    "Solar energy thermal processes" 
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