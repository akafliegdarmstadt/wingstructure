"""Module providing multhop quadrature method for solving 
prandtl' lifting line problem
"""

import numpy as np
from numpy import pi as π, sqrt, arctan
from collections import namedtuple, defaultdict


# Definition of low level functions 
_multhop_result = namedtuple('ext_multhopp_result', 
                              ('c_ls', 'α_is', 'C_L', 'C_Di'))


def _prepare_multhop(ys: np.ndarray, αs: np.ndarray, chords: np.ndarray,
             dcls: np.ndarray, S:float, b:float, M=None):
    
    # calculate aspect ratio
    Λ = b**2 / S

    # calculate grid points
    if M is None:
        M = int(round(Λ)*4-1) # has to be uneven, not more than 4*aspect ratio
    elif M%2 == 0:
        M += 1

    θs = np.linspace(np.pi/(M+1), M/(M+1)*np.pi, M)
    calc_ys = -b/2 * np.cos(θs)

    # interpolate input values
    
    αs_int = np.interp(np.abs(calc_ys), ys, αs)
    chords_int = np.interp(np.abs(calc_ys), ys, chords)
    dcls_int = np.interp(np.abs(calc_ys), ys, dcls)

    return calc_ys, θs, αs_int, chords_int, dcls_int


def _multhop_solve(θs, αs, chords, dcls, b, return_αi=False):

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

    #print(B)

    # calculation of local circulation
    γs = np.dot(np.linalg.inv(B), αs)

    if not return_αi:
        return γs
    else:
        α_is = Bb@γs
        return γs, α_is


def multhop(ys: np.ndarray, αs: np.ndarray, chords: np.ndarray,
             dcls: np.ndarray, S:float, b:float, M:int=None, do_prep=True):
    """Low level function for multhop quadrature calculation

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
    M: int
        number of grid points
    
    Returns
    -------
    tuple
        multhop results - (c_ls, α_is, C_L, C_Di))
         lift coefficients, induced angle of attack, 
         wing's lift coefficient, induced drag coefficient
    """

    # calculate aspect ratio
    Λ = b**2 / S

    if do_prep:
        solverinput = _prepare_multhop(ys, αs, chords, dcls, S, b, M)

        γs, α_is = _multhop_solve(*solverinput[1:], b, return_αi=True)

        θs = solverinput[1]
    else:
        M = len(ys)

        θs = np.arccos(-2* ys/b)
        #print(θs)
        γs, α_is = _multhop_solve(θs, αs, chords, dcls, b, return_αi=True)
    
    # calculate lift coefficient distritbution
    c_ls = 2*b/(np.array(chords)) * np.array(γs) # TODO: Check this formula

    # calculate overall lift coefficient (whole wing)
    C_L = π*Λ / (M+1) * np.sum(γs * np.sin(θs))

    # calculate induced drag
    C_Di = π*Λ/(M+1) * np.sum( γs * α_is * np.sin(θs))

    return _multhop_result(c_ls, α_is, C_L, C_Di)


# Helper functions for high level interface

class AirfoilData(object):
    def __init__(self, alpha0: float = 0, dif_ca_alpha: float = 2*np.pi, cm0: float = 0):
        self.alpha0 = alpha0
        self.dif_ca_alpha = dif_ca_alpha
        self.cm0 = cm0

def _calc_gridpoints(wing, M:int):

    b = wing.span

    θs = np.linspace(np.pi/(M+1), M/(M+1)*np.pi, M)
    calc_ys = -b/2 * np.cos(θs)

    return calc_ys

def _calc_base_α(wing, ys, airfoil_db):
    """calculates of geometric and aerodynamic twist"""

    def get_α0(airfoil):
        return np.radians(airfoil_db[airfoil].alpha0)

    αs = wing.twists # geometric twist

    αs -= np.array([get_α0(airfoil) for airfoil in wing.airfoils])

    return np.interp(np.abs(ys), wing.ys, αs)

def _calc_flap_Δα(controlsurface, ys, η):

    λ_k = 1 - controlsurface.chordpos_at(ys)

    λ_k[ys<0.0] = 0.0 # only consider one of two (symmetric wing)

    η_k = np.full_like(ys, η)

    η_keff = 22.743 * arctan(0.04715 * η_k)

    k = -2/π * (sqrt(λ_k * (1-λ_k)) + arctan(sqrt(λ_k)))

    αs = -(0.75 * k - 0.25 * λ_k) * η_keff

    return αs


# Definition of high level functions and analysis object

def calc_multhoplift(wing, α, controls:dict={}, airbrake:bool=False, airfoil_db:dict=None, options:dict=None):
    
    # use default airfoil if no database is set
    if airfoil_db is None:
        #warn('No airfoil database defined, using default airfoil.')
        airfoil_db = defaultdict(AirfoilData)

    # calculate grid points
    ys = _calc_gridpoints(wing, 87)

    # sum up aoa parts

    # geometric and aerodynamic twist
    αs = _calc_base_α(wing, ys, airfoil_db) 

    # control surface parts
    for name, (η_r, η_l) in controls.items():
        try:
            control_surf = wing.flaps[name]
        except:
            raise Exception('control surface "{}" not defined in wing!'.format(name))
            
        #αs += _calc_flap_Δα(control_surf, ys, η_r)
        #αs += _calc_flap_Δα(control_surf, ys, η_l)[::-1]

    # - airbrake
    if airbrake:
        #αs += np.where(wing.within_airbrake(ys), -12.0, 0.0)
        pass

    # add aoa
    αs += α

    # interpolate chord lengthes
    chords = np.interp(np.abs(ys), wing.ys, wing.chords)

    # determine lift slope
    dcls = np.full_like(chords, 2*π) #TODO: make adaptable

    res = multhop(ys, αs, chords, dcls, wing.area, wing.span, do_prep=False)

    return {'ys': ys, 'c_ls':res.c_ls, 'C_L': res.C_L, 'C_Di': res.C_Di}