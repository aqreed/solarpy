from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="solarpy",
    author='aqreed',
    description='Solar radiation model based on Duffie & Beckman\
                "Solar energy thermal processes" (1974)',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version="0.1",
    url='https://github.com/aqreed/solarpy',
    packages=['solarpy'],
    install_requires=['numpy', 'matplotlib'],
    tests_requires=['pytest']
    )
