import numpy as np
from collections import defaultdict, namedtuple

from wingstructure.aero import multhop, AirfoilData


defaultairfoil_db = defaultdict(AirfoilData)


# test low level functions

def test_multhop_solve():
    """test multhop_solve with numeric data from schlichting book    
    """

    Λ = 6 # aspect ratio of wing
    
    b = 15 # m span width
    S = b**2/Λ

    cs = np.array([b/Λ]*2) # depth of wing
    ys = np.array([0,b/2]) # section positions

    αs = np.array([1]*2) # angle of attack
    dcls = np.array([2*np.pi]*2)

    # reference result from Schlichting
    ηs_ref = np.array([0,0.3827,0.7071,0.9239,1])
    γs_ref = np.array([0.4320,0.4192,0.3710,0.2485,0])

    # coarse calculation
    M = 7

    solverinput = multhop._prepare_multhop(ys, αs, cs, dcls, S, b, M)

    calc_ys = solverinput[0]

    γs = multhop._multhop_solve(*solverinput[1:], b)

    assert np.isclose(calc_ys[M//2:]/b*2, ηs_ref[:-1], atol=1e-4).all()
    assert np.isclose(γs[M//2:], γs_ref[:-1], atol=1e-4).all()


# test helper functions 

def test_calc_gridpoints():
    """basic tests for multhop's grid point generation"""

    from wingstructure.aero.multhop import _calc_gridpoints
    
    Wing = namedtuple('Wing', ['span', 'aspectratio'])
    wing = Wing(span=9, aspectratio=15)
    
    # span
    b = wing.span
    
    ys = _calc_gridpoints(wing)

    # grid points should be sorted
    assert (np.sort(ys) == ys).all()

    # check if values are within span
    assert (ys<b/2.0).all()
    assert (ys>-b/2.0).all()


def test_calc_base_α():
    """check base α for geometric and aerodynamic twist"""

    from wingstructure.aero.multhop import _calc_base_α

    twists = np.zeros(5)

    Wing = namedtuple('Wing', ['airfoils', 'twists', 'ys'])

    wing = Wing(airfoils=5*[''], twists=twists, ys=[0.0, 1.0, 2.5, 5.0, 7.5])

    res = _calc_base_α(wing, [-7.2, -5, -3, 0.0, 3, 5, 7.2], defaultairfoil_db)

    assert res.shape[0]==7

    assert (res == 0).all()


def test_control_surface_Δα():
    """check α of control surface"""
    
    from wingstructure.aero.multhop import _calc_flap_Δα

    ControlSurface = namedtuple('ControlSurface', ['pos1', 'pos2', 'depth1', 'depth2'])

    length = lambda span_pos: np.interp(span_pos, [3.5, 5.5], [0.8, 0.9],
                right=0.0, left=0.0)

    cs = ControlSurface(3.5, 5.5, 0.8, 0.9)

    res = _calc_flap_Δα(cs, np.array([1.0, 2.5, 3.5, 5.0, 7.5]), 5)

    # outside of controlsurface no change in angle of attack
    assert (res[0:2] == 0).all()
    assert (res[-1] == 0).all()



