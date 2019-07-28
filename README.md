# pysolarRadiation
|  |  |
| ------ | ------ |
| Description | Python Solar Radiation model |
| Author | AeroPython Team <aeropython@groups.io> |
| Version | 0.1.dev0 |
| Python Version | 3.6 |
| Requires | Numpy, Matplotlib, scikit-aero |

This packages aims to provide with a reliable solar radiation model, mainly based on the work of Duffie and Beckman "Solar energy thermal processes" (1974). 

The main objective is generate a **solar beam irradiance** (W/m2) prediction on:
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
$ cd pysolarRadiation
$ pip install -e .
```

### Testing

pysolarRadiation recommends py.test for running the test suite. Running from the top directory:

```sh
$ pytest
```
### License

MIT (see `COPYING`)