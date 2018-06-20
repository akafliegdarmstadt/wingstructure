# -*- coding: utf-8 -*-
"""
A module providing multhop method for solving lifting line problem
"""
import numpy as np
from numba import jit

π = np.pi

@jit
def multhopp(αs, θs, chords, b, dcls):
    """Calculates lift distribution with multhopp method.
    
    :param αs: angle of attack in radians at grid points
    :param θs: gridpoints in transformed spanwise coordinates y=b/2 cos(θ)
    :param chords: chord lenghtes at grid points
    :param dcls: lift slope of wing sections at grid points
    :rtype: dict
    """
    M = len(θs)
    
    # create empyt matrix (N_MxN_M) for multhoppcoefficients
    B = np.zeros((M, M))

    # calculation of multhopp coefficients
    for v, (theta_v, c, dcl_v) in enumerate(zip(θs, chords, dcls)):
        for n, theta_n in enumerate(θs):

            # diagonal elements
            if v == n:
                B[v, v] = (M+1)/(4*np.sin(theta_v))+2*b/(dcl_v*c)
            # non diagonal elements
            else:
                B[v, n] = - ((1-(-1.)**(v-n))/2*(np.sin(theta_n)/ \
                            ((M+1)*(np.cos(theta_n)-np.cos(theta_v))**2)))

    # calculation of local circulation
    γs = np.dot(np.linalg.inv(B), αs)

    return γs

def calculate_liftcoefficients(Θs, γs, chords, Λ, b, M):
    """
    Calculates lift coefficients with multhopp results
    """

    c_l = 2*b/(np.array(chords)) * np.array(γs)
    C_L = π*Λ / (M+1) * np.sum(γs * np.sin(Θs))

    return C_L, c_l
