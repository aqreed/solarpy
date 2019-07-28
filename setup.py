from setuptools import setup, find_packages

setup(
    name="solarradiation",
    author='aqreed',
    version="0.1.dev0",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=['numpy', 'matplotlib', 'scikit-aero'],
    tests_requires=['pytest']
    )