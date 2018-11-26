import sys, os
from setuptools import setup, find_packages

setup(name='wingstructure',
    version='0.0.2',
    description='',
    url='None',
    author='helo',
    author_email='-',
    license='GPL 2.0',
    packages=find_packages('.'),
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
