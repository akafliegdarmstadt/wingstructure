from wingstructure import wingloads
from .wings import d38wing


def test_calculationpoints(d38wing):
    
    # calculate base distribution
    ys = wingloads.calculation_points(d38wing, 40)

    # check for existance of all sections, exept last one
    assert all((y in ys for y in d38wing.ys[:-1]))
    assert all((-y in ys for y in d38wing.ys[:-1]))