from setuptools import setup, find_packages

setup(
    name="solarradiation",
    author='aqreed',
    description='Solar radiation model based on Duffie & Beckman\
                "Solar energy thermal processes" (1974)',
    version="0.1.dev0",
    url='https://github.com/aqreed/pysolarRadiation',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    tests_requires=['pytest']
    )
