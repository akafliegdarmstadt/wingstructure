import pytest
import numpy as np

from .wings import d38wing, d43wing

# helper functions 

def isclose(a, b, tol=1e-5, rtol=1e-5):
    diff = abs(a-b)
    maxab = max((abs(a), abs(b)))

    return diff<tol and (diff/maxab)<rtol

# tests


def test_wingbasicprops(d38wing):
    assert isclose(d38wing.area, 11.03516)
    assert d38wing.span==15.0
    assert isclose(d38wing.aspectratio, 20.3893736)
    # compare mac results with getmac
    macpos, maclength = d38wing.get_mac()
    assert isclose(maclength, 0.7706271)
    assert isclose(macpos[0], 0.2109074-0.7706271/4)


def test_properties(d43wing):
    # check array values with input data
    assert (d43wing.chords == [1.12, 1.028, 0.673, 0.454537, 0.36]).all()
    assert (d43wing.ys == [0.0, 4.0, 7.223, 8.5, 9.0]).all()
    assert (d43wing.twists == 0.0).all()
    assert (d43wing.airfoils == '').all()


def test_heperfunctions(d43wing):
    d43wing.add_controlsurface('aileron1', 4.0, 8.5, 0.8, 0.8, 'aileron')
    d43wing.add_controlsurface('airbrake1', 2.4, 3.83, 0.5, 0.5, 'airbrake')

    assert d43wing.within_control('airbrake1', 2.5)
    assert d43wing.within_control('aileron1', 4.5)

    assert (d43wing.within_airbrake(np.array([1.0, 2.5, 4.0])) == np.array([False, True, False])).all()

    # except Exception with wrong control name
    with pytest.raises(KeyError):
        d43wing.within_control('airbrake15', 3.3)


def test_serialization(d43wing):
    from wingstructure.data.wing import Wing

    # add control surfaces
    d43wing.add_controlsurface('aileron1', 4.0, 8.5, 0.8, 0.8, 'aileron')
    d43wing.add_controlsurface('airbrake1', 2.4, 3.83, 0.5, 0.5, 'airbrake')

    # serialize
    data = d43wing.serialize()
   
    # rebuild wing
    d43wing2 = Wing.deserialize(data)

    # check that wings are equal
    assert (d43wing.ys == d43wing2.ys).all()
    assert (d43wing.chords == d43wing2.chords).all()
    assert (d43wing.twists == d43wing2.twists).all()
    assert (d43wing.airfoils == d43wing2.airfoils).all()
    
    assert d43wing.controlsurfaces['airbrake1'].pos1 == d43wing2.controlsurfaces['airbrake1'].pos1

    # control surfaces are optional, should not raise any Exception if missing
    del data['controlsurfaces']
    Wing.deserialize(data)


def test_plot(d43wing):
    d43wing.plot()