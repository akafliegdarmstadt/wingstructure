import sys, os
from setuptools import setup, find_packages

setup(name='wingstructure',
    version='0.0.2',
    description='A library for structure calculations in airplane wings',
    url='https://github.com/helo9/wingstructure',
    author='helo',
    author_email='kontakt@akaflieg.tu-darmstadt.de',
    license='MIT',
    packages=['wingstructure'],
    install_requires= [
        'numpy',
        'scipy',
        'sortedcontainers',
        'matplotlib',
        'pytest',
        'strictyaml',
        'shapely'
    ],
    zip_safe=False
    )
