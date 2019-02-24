# -*- coding: utf-8 -*-
"""
Module providing the multhop quadrature method for solving prandtl's lifting line problem.
"""
import numpy as np
#from numba import jit
from collections import namedtuple


π = np.pi


_multhopp_result = namedtuple('multhopp_result', ('ys','c_ls','C_L'))
_ext_multhopp_result = namedtuple('ext_multhopp_result', ('ys', 'c_ls', 'C_L', 'αᵢs','C_Wi'))


#@jit
def _solve_multhopp(αs, θs, chords, b, dcls, return_B = False):

    M = len(θs)

    # create empyt matrix (N_MxN_M) for multhoppcoefficients
    Bb = np.zeros((M, M))
    Bd = np.zeros(M)

    # calculation of multhopp coefficients
    for v, (θv, c, dcl_v) in enumerate(zip(θs, chords, dcls)):
        for n, θn in enumerate(θs):

            # diagonal elements
            if v == n:
                # prevent division throught zero
                if np.isclose(np.sin(θv), 0.0):
                    θv = 1e-15
                
                Bb[v, v] = (M+1)/(4*np.sin(θv))
                Bd[v] = 2*b/(dcl_v*c)
            # non diagonal elements
            else:
                Bb[v, n] = - ((1-(-1.)**(v-n))/2*(np.sin(θn)/ \
                            ((M+1)*(np.cos(θn)-np.cos(θv))**2)))

    B = Bb + np.diag(Bd)

    # calculation of local circulation
    γs = np.dot(np.linalg.inv(B), αs)

    if return_B:
        return γs, Bb, Bd
    else:
        return γs

def _calculate_liftcoefficients(Θs, γs, chords, Λ, b, M):
    """Calculates lift coefficients from circulation
    
    Parameters
    ----------
    chords : np.ndarray
        chord lengthes
    b : float
        span width
    M : int
        number of grid points

    Returns
    -------
    tuple
        c_l (ndarray) - local lift coefficents, C_L (float) - lift coefficient wing
    """


    # calculate lift coefficient distritbution
    c_l = 2*b/(np.array(chords)) * np.array(γs) # TODO: Check this formula

    # calculate overall lift coefficient (whole wing)
    C_L = π*Λ / (M+1) * np.sum(γs * np.sin(Θs))

    return c_l, C_L

def _calcgridpoints(b:float, Λ:float, M:int = None):

    # calculate number of gridpoints
    if M is None:
        M = int(round(Λ)*4-1) # has to be uneven, not more than 4*aspect ratio
    elif M%2 == 0:
        M += 1 # has to be uneven

    # grid points as angle
    θs = np.linspace(np.pi/(M+1), M/(M+1)*np.pi, M)

    # calculate grid points
    calc_ys = -b/2 * np.cos(θs)

    return θs, calc_ys

def multhopp(αs: np.ndarray, chords: np.ndarray, ys: np.ndarray, dcls: np.ndarray or float=np.nan, M:int=None,
             mode = 'c_l', interp = True ):
    """Use multhopp's quadrature to solve prandtl's lifting line problem
    
    Parameters
    ----------
    αs : np.array
        angle of attack regarding zero lift angle
    chords : np.array
        chord lengths of wing
    ys : np.array
        corresponding span position
    dcls : np.array, optional
        lift coefficient slope (the default is np.nan, which means 2pi is used)
    M : int, optional
        number of evaluation points (the default is None, which let function use recommendation number)
    mode : str, optional
         (the default is 'c_l', which [default_description])
    interp : bool, optional
        calculate grid points or use given span positions for calculation (the default is True)
    
    Returns
    -------
    namedtuple
        local lift coefficients, optional circulation distribution
    """

    if np.isnan(dcls).all():
        dcls = np.array([2*π]*len(ys))

    # calculate wingspan
    b = 2*max(ys)

    # calculate wing area
    if min(ys) >= 0.0:
        S = 2 * np.trapz(y=chords, x=ys)
    else:
        S = np.trapz(y=chords, x=ys)
    
    # calculate aspect ratio
    Λ = b**2 / S

    # interpolate
    if interp: 
        θs, calc_ys = _calcgridpoints(b, Λ, M)
        
        calc_αs = np.interp(np.abs(calc_ys), ys, αs)
        calc_chords = np.interp(np.abs(calc_ys), ys, chords)
        calc_dcls = np.interp(np.abs(calc_ys), ys, dcls)
    else:
        calc_ys = ys
        calc_chords = chords
        calc_αs = αs
        calc_dcls = dcls

        θs = np.arccos(-2*calc_ys/b)

    # calculate circulation distribution
    γs, Bb, Bd = _solve_multhopp(calc_αs, θs, calc_chords, b, calc_dcls, return_B=True)

    B = Bb+np.diag(Bd)

    # return only ys and gamma
    if mode in ('gamma', 'γ'):
        return calc_ys, γs
    
    # calculate coefficients
    c_ls, C_L = _calculate_liftcoefficients(θs, γs, calc_chords, Λ, b, M)

    if mode == 'c_l':
        # return lift coefficients
        return _multhopp_result(calc_ys, c_ls, C_L)
    elif mode == 'combined':
        # return lift coefficients, induced angle and induced drag

        αᵢs = Bb@γs
        C_Di = π*Λ/(M+1) * np.sum( γs * αᵢs * np.sin(θs))
        
        return _ext_multhopp_result(calc_ys, c_ls, C_L, αᵢs, C_Di)
