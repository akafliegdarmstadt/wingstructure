from wingstructure.data.geometry import Wing, Point
import pytest

# test sorting of sections


def test_simple_wing_creation():

    wing = Wing()

    wing.add_section(Point(0, 0, 0), 1)
    wing.add_section(Point(0, 1, 0), 1)

    assert wing.area == 2.
    assert wing.span == 2.
