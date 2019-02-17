"""Module providing multhop quadrature method for solving 
prandtl' lifting line problem
"""

import numpy as np
from collections import namedtuple


# Define Pi
π = np.pi


# Definition of low level functions 

_multhop_result = namedtuple('ext_multhopp_result', 
                              ('c_ls', 'α_is', 'C_L', 'C_Di'))


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


def _multhop_solve(αs, θs, chords, dcls, b, return_αi=False):

    # number of grid points
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

    if not return_αi:
        return γs
    else:
        α_is = Bb@γs
        return γs, α_is


def multhop(ys: np.ndarray, αs: np.ndarray, chords: np.ndarray,
             dcls: np.ndarray, S:float, b:float):
    """Low level for multhop quadrature calculation

    The parameters except for S and b have to be numpy arrays of the
    same length and determine the wing geometry as well calculation
    grid points.

    
    
    Parameters
    ----------
    
    ys : np.ndarray
        span positions at which lift is calculated
    chords : np.ndarray
        chord lengthes at span positions
    αs: np.ndarray
        array of angles of attack for chord positions
    dcls : np.ndarray
        lift coefficient slope regarding angle of attack
    S : float
        wing area
    b : float
        span width of wing
    
    Returns
    -------
    tuple
        multhop results - (c_ls, α_is, C_L, C_Di))
         lift coefficients, induced angle of attack, 
         wing's lift coefficient, induced drag coefficient
    """

    # number of grid points
    M = len(αs)

    # calculate aspect ratio
    Λ = b**2 / S

    # multhop quadrature
    θs = np.arccos(-2*ys/b)

    γs, α_is = _multhop_solve(αs, θs, chords, b, dcls, return_αi=True)

    # calculate lift coefficient distritbution
    c_ls = 2*b/(np.array(chords)) * np.array(γs) # TODO: Check this formula

    # calculate overall lift coefficient (whole wing)
    C_L = π*Λ / (M+1) * np.sum(γs * np.sin(θs))

    # calculate induced drag
    C_Di = π*Λ/(M+1) * np.sum( γs * α_is * np.sin(θs))

    return _multhop_result(c_ls, α_is, C_L, C_Di)


# Definition of high level functions and object

def calc_multhoplift()


