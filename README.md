# solarpy

[![Build Status](https://travis-ci.com/aqreed/solarpy.svg?branch=master)](https://travis-ci.com/aqreed/solarpy)
[![codecov.io](https://codecov.io/gh/aqreed/solarpy/branch/master/graph/badge.svg)](https://codecov.io/gh/aqreed/solarpy/branch/master)
[![license](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://github.com/aqreed/solarpy/raw/master/COPYING)

|  |  |
| ------ | ------ |
| Description | Python Solar Radiation model |
| Author | aqreed <aqreed@protonmail.com> |
| Version | 0.1 |
| Python Version | 3.6 |
| Requires | Numpy, Matplotlib |

This packages aims to provide a reliable solar radiation model, mainly based on the work of Duffie, J.A., and Beckman, W. A., 1974, "Solar energy thermal processes".

The main purpose is to generate a **solar beam irradiance** (W/m2) prediction on:
* **any plane**, thanks to the calculation of the solar vector in NED (North East Down) coordinates, suitable for its use in flight dynamics simulations...
* **any place of the earth**, taking into account the solar time wrt the standard time, geometric altitude, the latitude influence on solar azimuth and solar altitude as well as sunset/sunrise time and hour angle, etc.
* **any day of the year**, taking into account the variations of the extraterrestrial radiation, the equation of time, the declination, etc., throughout the year

Solar [irradiance](https://en.wikipedia.org/wiki/Solar_irradiance) on the southern hemisphere on October 17, at sea-level 13.01UTC (plane pointing upwards)?

```
from solarpy.radiation import irradiance_on_plane
from solarpy.utils import day_of_the_year
import numpy as np

vnorm = np.array([0, 0, -1])  # plane pointing zenith
h = 0  # sea-level
n = day_of_the_year(10, 17)  # October 17
lat = -23.5  # southern hemisphere
hour, minute = 13, 1  # midday

irradiance_on_plane(vnorm, h, n, lat, hour, minute)
```

A dedicated Jupyter Notebook on beam irradiance can be found [here](https://github.com/aqreed/solarpy/blob/master/examples/solar_irradiance.ipynb).

Solar [declination](https://en.wikipedia.org/wiki/Position_of_the_Sun#Declination_of_the_Sun_as_seen_from_Earth) on August 5?

```
from solarpy.radiation import declination
from solarpy.utils import day_of_the_year

n = day_of_the_year(8, 5)  # August 5

declination(n)
```

Please find more notebooks on the 'examples' folder.

---
**NOTE**:
solarpy is under development and might change in the near future.

---

### Dependencies

This package depends on Python, NumPy and Matplotlib and is usually tested on Linux with the following versions:

Python 3.6, NumPy 1.16, Matplotlib 3.0

### Installation

solarpy has been written in Python3

```sh
$ git clone https://github.com/aqreed/solarpy.git
$ cd solarpy
$ pip install -e .
```

### Testing

solarpy recommends py.test for running the test suite. Running from the top directory:

```sh
$ pytest
```

To test coverage (also from the top directory):

```sh
$ pytest --cov
```

### Bug reporting

Please feel free to open an [issue](https://github.com/aqreed/solarpy/issues) on GitHub!

### License

MIT (see `COPYING`)