# -*- coding: utf-8 -*-
"""
A module providing multhop method for solving lifting line problem
"""
import numpy as np
from numba import jit
from collections import namedtuple


π = np.pi


_multhopp_result = namedtuple('multhopp_result', ('ys','c_ls','C_L'))


@jit
def _solve_multhopp(αs, θs, chords, b, dcls):
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

def _calculate_liftcoefficients(Θs, γs, chords, Λ, b, M):
    """
    Calculates lift coefficients with multhopp results
    """

    # calculate lift coefficient distritbution
    c_l = 2*b/(np.array(chords)) * np.array(γs)

    # calculate overall lift coefficient (whole wing)
    C_L = π*Λ / (M+1) * np.sum(γs * np.sin(Θs))

    return c_l, C_L

def multhopp(αs: np.array, chords: np.array, ys: np.array, dcls: np.array=None, M:int=None,
             return_γ = False ):

    # calculate wingspan
    b = 2*max(ys)

    # calculate wing area
    S = 2 * np.trapz(y=chords, x=ys)
    
    # calculate aspect ratio
    Λ = b**2 / S

    # calculate number of gridpoints
    if M is None:
        M = int(round(Λ)*4-1) # has to be uneven, not more than 4*aspect ratio
    elif M%2 == 0:
        M += 1 # has to be uneven

    if not dcls:
        dcls = [2*np.pi]*len(chords)

    # grid points as angle
    θs = np.linspace(np.pi/(M+1), M/(M+1)*np.pi, M)

    # calculate grid points
    calc_ys = -b/2 * np.cos(θs)

    # interpolate
    calc_αs = np.interp(np.abs(calc_ys), ys, αs)
    calc_chords = np.interp(np.abs(calc_ys), ys, chords)
    calc_dcls = np.interp(np.abs(calc_ys), ys, dcls)

    # calculate circulation distribution
    γs = _solve_multhopp(calc_αs, θs, calc_chords, b, calc_dcls)

    if return_γ:
        return calc_ys, γs
    
    # calculate coefficients
    c_ls, C_L = _calculate_liftcoefficients(θs, γs, calc_chords, Λ, b, M)

    return _multhopp_result(calc_ys, c_ls, C_L)
