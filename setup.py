import sys, os
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='wingstructure',
    version='0.0.4',
    description='A library for structure calculations in airplane wings',
    url='https://github.com/akafliegdarmstadt/wingstructure',
    author='helo',
    author_email='kontakt@akaflieg.tu-darmstadt.de',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    packages=['wingstructure', 
              'wingstructure.data', 
              'wingstructure.aero',
              'wingstructure.structure'],
    install_requires= [
        'numpy',
        'scipy',
        'pandas',
        'sortedcontainers',
        'matplotlib',
        'pytest',
        'strictyaml',
        'Shapely',
        'ipython',
        'pandas'
    ],
    zip_safe=False
    )
