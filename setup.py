from setuptools import setup, find_packages

setup(
    name="solarpy",
    author='aqreed',
    description='Solar radiation model based on Duffie & Beckman\
                "Solar energy thermal processes" (1974)',
    version="0.1dev1",
    url='https://github.com/aqreed/solarpy',
    packages=['solarpy'],
    install_requires=['numpy', 'matplotlib'],
    tests_requires=['pytest']
    )
