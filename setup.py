import sys, os
from setuptools import setup, find_packages

# taken from http://www.pydanny.com/python-dot-py-tricks.html
if sys.argv[-1] == 'test':
    test_requirements = [
        'pytest',
    ]
    try:
        modules = map(__import__, test_requirements)
    except ImportError as e:
        err_msg = e.message.replace("No module named ", "")
        msg = "%s is not installed. Install your test requirments." % err_msg
        raise ImportError(msg)
    os.system('py.test')
    sys.exit()
# -------------------------------------------

setup(name='wingstructure',
    version='0.0.1',
    description='',
    url='None',
    author='helo',
    author_email='-',
    license='GPL 2.0',
    packages=['wingstructure'],
    install_requires= [
        'numpy',
        'shapely',
        'svgwrite'
    ],
    zip_safe=False
    )
