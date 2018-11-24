import sys, os
from setuptools import setup, find_packages

setup(name='wingstructure',
    version='0.0.2',
    description='',
    url='None',
    author='helo',
    author_email='-',
    license='GPL 2.0',
    packages=['wingstructure'],
    install_requires= [
        'numpy',
        'sortedcontainers',
        'matplotlib',
        'pytest'
    ],
    zip_safe=False
    )
