import numpy as np
from wingstructure.geometry import Wing, Point
import pytest

# test sorting of sections

def test_simple_wing_creation():

    wing = Wing()

    wing.add_section(Point(0,0,0), 1)
    wing.add_section(Point(0,1,0), 1)

    assert wing.area==2.
    assert wing.span==2.

def test_section_sorting():
    
    wing = Wing()

    # sections should be sorted according y-position
    wing.add_section(Point(0,2,0), 2.)
    wing.add_section(Point(0,3,0), 1.)
    wing.add_section(Point(0,1,0), 3.)

    # -> constant decreasing chord length?
    assert np.all(np.diff(wing.chords) < 0 )