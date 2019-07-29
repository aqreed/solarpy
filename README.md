# pysolarRadiation
|  |  |
| ------ | ------ |
| Description | Python Solar Radiation model |
| Author | aqreed |
| Version | 0.1.dev0 |
| Python Version | 3.6 |
| Requires | Numpy, Matplotlib, scikit-aero |

This packages aims to provide a reliable solar radiation model, mainly based on the work of Duffie, J.A., and Beckman, W. A., 1974, "Solar energy thermal processes".

The main purpose is to generate a **solar beam irradiance** (W/m2) prediction on:
* **any plane**, thanks to the calculation of the solar vector in NED (North East Down) coordinates, suitable for its use in flight dynamics simulations...
* **any place of the earth**, taking into account the solar time wrt the standard time, the latitude effect on solar azimuth and altitude as well as sunset/sunrise time and hour angle, etc.
* **any day of the year**, taking into account the variations of the extraterrestrial radiation, the equation of time, the declination, etc., throughout the year

### Dependencies

This package depends on Python, NumPy and scikit-aero and is usually tested on Linux with the following versions:

Python 3.6, NumPy 1.16, scikit-aero0.2.dev0

### Installation

pysolarRadiation has been written in Python3

```sh
$ git clone https://github.com/aqreed/pysolarRadiation.git
$ cd pysolarRadiation
$ pip install -e .
```

You will also need skaero:

```sh
$ git clone https://github.com/AeroPython/scikit-aero.git
$ cd scikit-aero
$ pip install -e .
```

### Testing

pysolarRadiation recommends py.test for running the test suite. Running from the top directory:

```sh
$ pytest
```

### Usage

Solar [declination](https://en.wikipedia.org/wiki/Position_of_the_Sun#Declination_of_the_Sun_as_seen_from_Earth) on August 5?

```
import solar_radiation as sr

n = sr.day_of_the_year(8, 5)  # August 5
dec = sr.declination(n)

print(dec)
```

Please find some Jupyter Notebooks on the 'examples' folder.

### License

MIT (see `COPYING`)