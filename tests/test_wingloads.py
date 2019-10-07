from wingstructure import wingloads
from .wings import d38wing, d43wing

# helper functions 

def isclose(a, b, tol=1e-5, rtol=1e-5):
    diff = abs(a-b)
    maxab = max((abs(a), abs(b)))

    return diff<tol and (diff/maxab)<rtol


# tests

def test_calculationpoints(d38wing):
    
    # calculate base distribution
    ys = wingloads.calculation_points(d38wing, 40)

    # check for existance of all sections, exept last one
    assert all((y in ys for y in d38wing.ys[:-1]))
    assert all((-y in ys for y in d38wing.ys[:-1]))


def test_flatwing(d43wing):
    from wingstructure.wingloads import FlatWing

    flatwing = FlatWing(d43wing)

    # check new spanwidth
    assert isclose(flatwing.span, 2*9.0408708357794)
    # area should have increased
    assert flatwing.area > d43wing.area
    # mac should decrease
    assert flatwing.mac < d43wing.mac