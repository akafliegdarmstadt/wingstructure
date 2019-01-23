import numpy as np
from wingstructure.data import geometry
import pytest


def test_simple_wing_creation():

    wing = geometry.Wing()

    wing.add_section(geometry.Point(0, 0, 0), 1)
    wing.add_section(geometry.Point(0, 1, 0), 1)

    assert wing.area == 2.
    assert wing.span == 2.

def test_section_sorting():

    wing = geometry.Wing()

    wing.add_section(geometry.Point(0,1,0), 1)
    wing.add_section(geometry.Point(0,0,0), 1)

    assert wing.area == 2.
    assert wing.span == 2.
    assert wing.sections[0].y == 0.


def test_create_from_dict():
    from . import wingdict as wd
    
    wing = geometry.Wing.create_from_dict(wd.wingdict)

    assert wing is not None


def test_mac():
    import numpy as np

    wing = geometry.Wing()

    # F-117 example from getmac
    wing.add_section(geometry.Point(0.0, 0.0, 0.0), 1120)
    wing.add_section(geometry.Point(344, 200, 0), 400)
    wing.add_section(geometry.Point(645, 375, 0), 400)
    wing.add_section(geometry.Point(860, 500, 0), 0)

    mac = wing.mac

    # getmac results
    assert wing.area/2 == 247e3
    assert np.isclose(mac.pos.y, 175.4723)
    assert np.isclose(mac.chord, 643.0229)
    
