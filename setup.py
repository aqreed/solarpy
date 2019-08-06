from setuptools import setup, find_packages

setup(
    name="pysolarradiation",
    author='aqreed',
    description='Solar radiation model based on Duffie & Beckman\
                "Solar energy thermal processes" (1974)',
    version="0.2.dev0",
    url='https://github.com/aqreed/pysolarRadiation',
    packages=['pysolar'],
    install_requires=['numpy', 'matplotlib', 'scipy',
                      'scikit-aero@git+https://git@github.com/AeroPython/scikit-aero@master'],
    tests_requires=['pytest']
    )
