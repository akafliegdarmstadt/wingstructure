import numpy as np

from wingstructure.aero import multhop


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
