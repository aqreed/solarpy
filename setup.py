from setuptools import setup, find_packages

setup(
    name="pysolarradiation",
    author='aqreed',
    description='Solar radiation model based on Duffie & Beckman\
                "Solar energy thermal processes" (1974)',
    version="0.2.dev0",
    url='https://github.com/aqreed/pysolarRadiation',
    packages=['pysolar'],
    tests_requires=['pytest']
    )
