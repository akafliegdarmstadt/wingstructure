"""Module for structural section analysis based on polygon geometry definition.

Inside this module basis formulas for calculation of area, static moments and 
moments of inertia based on polygon geometry are defined. Furthermore functions
for calculation of the derived quanities neutral center and principal axis angle
are included. For a simplified application to sections defined with *section*
module the class *StructuralAnalysis* is also included. 

.. HINT::
   All coordinate sequences passed directly to the low level functions (not the class)
   have to be defined in a clockwise manner (negative mathematical direction).
"""

import numpy as np


def calcarea(outline):
    """Calculate area within polygon
    
    Parameters
    ----------
    outline : np.array
        2D coordinates of a polygon area
    
    Returns
    -------
    float
        area within outline
    """

    
    x_i, y_i = outline.T
    x_ip1, y_ip1 = np.roll(outline.T, 1, axis=1)
    
    A = 0.5 * np.sum(y_ip1*x_i-y_i*x_ip1)
    
    return A


def calcstaticmoments(outline):
    """Calculate the static moment of a polygon
    
    Parameters
    ----------
    outline : np.array
        2D coordinates of a polygon area
    
    Returns
    -------
    tuple(float)
        static moments
    """

    x_i, y_i = outline.T
    x_ip1, y_ip1 = np.roll(outline.T, 1, axis=1)
    
    S_x = 1/6 * np.sum((y_i+y_ip1)*(y_ip1*x_i-y_i*x_ip1))
    S_y = 1/6 * np.sum((x_i+x_ip1)*(y_ip1*x_i-y_i*x_ip1))
    
    return S_x, S_y


def calcinertiamoments(outline):
    """Calculate moment of inertia for given polygon
    
    Parameters
    ----------
    outline : np.array
        2D coordinates of a polygon area
    
    Returns
    -------
    tuple(float)
        moments of inertia (I_xx, I_yy, I_xy)
    """

    x_i, y_i = outline.T
    x_ip1, y_ip1 = np.roll(outline.T, 1, axis=1)
    
    I_xx = 1/12 * np.sum((y_ip1**2 + (y_i+y_ip1)*y_i)\
                         *(y_ip1*x_i-y_i*x_ip1))
    I_yy = 1/12 * np.sum((x_ip1**2 + (x_i+x_ip1)*x_i)\
                         *(y_ip1*x_i-y_i*x_ip1))
    I_xy = 1/12 * np.sum(0.5*x_ip1**2*y_i**2-0.5*x_i**2*y_ip1**2\
                          -(y_ip1*x_i-y_i*x_ip1)*(x_i*y_i+x_ip1*y_ip1))
    
    return I_xx, I_yy, I_xy


def calcneutralcenter(outline):
    """Calculate the neutral center of polygon
    
    Parameters
    ----------
    outline : np.array
        2D coordinates of a polygon area
    
    Returns
    -------
    tuple(float)
        (x_n, y_n)
    """

    
    area_ = calcarea(outline)
    
    S_x, S_y = calcstaticmoments(outline)
    
    x_nc = S_y/area_
    y_nc = S_x/area_
    
    return x_nc, y_nc


def calcprincipalaxis(I_xx, I_yy, I_xy):
    """Calculates angle of first principal axis
    
    Parameters
    ----------
    I_xx : float
        area moment of inertia
    I_yy : float
        area moment of inertia
    I_xy : float
        area moment of inertia
    
    Returns
    -------
    float
        principal axis angle
    """

    return np.arctan2(2*I_xy, I_yy-I_xx)


def transform_intertiamoments(I_xx, I_yy, I_xy, φ):
    """Transform moments of inertia from into rotated coordinate system
    
    Parameters
    ----------
    I_xx : float
        area moment of inertia
    I_yy : float
        area moment of inertia
    I_xy : float
        area moment of inertia
    φ: float
        angle to roatet
    
    Returns
    -------
    tuple of floats
        moment of inertia in rotated coordinate system
    """

    I_ξξ = (I_xx+I_yy)/2 + (I_yy-I_xx)/2 * np.cos(2*φ) - I_xy * np.sin(2*φ)
    I_ηη = (I_xx+I_yy)/2 - (I_yy-I_xx)/2 * np.cos(2*φ) + I_xy * np.sin(2*φ)
    I_ξη = (I_yy-I_xx)/2 * np.cos(2*φ) + I_xy * np.sin(2*φ)

    return I_ξξ, I_ηη, I_ξη


def _calcgeom(geom, property_functions):
    """helper function calculating asked properties for geometry
    """

    properties = {propfun:None for propfun in property_functions}
  
    isexterior = True
    for outline in (geom.exterior, *geom.interiors):
        if isexterior:
            for propfun in property_functions:
                properties[propfun] = np.array(propfun(outline))
            isexterior = False
        else:
            for propfun in property_functions:
                properties[propfun] -= np.array(propfun(outline))

    return [properties[propfun] for propfun in property_functions]


class StructuralAnalysis:
    """Class for structural analysis of Sections or Features

    Parameters
    ----------
    structure: SectionBase or Feature
        structure element to be analysed

    Attributes
    ----------
    nc : np.array
        Neutralcenter of structure, available after calling update
    area : float
        Overall Area of structure, available after calling update
    bendingstiffness : np.array
        bending stiffness of structure, available after calling update
    """

    def __init__(self, structure):
        self._structure = structure
        self.nc = None
        self.area = None
        self.bendingstiffness = None

    def update(self):
        # calculate normal center

        self.area = 0

        weighted_area_ = 0
        weighted_staticmoment_ = np.zeros((2))

        for geom in self._structure.exportgeometry():
            area, staticmoment = _calcgeom(geom, [calcarea, calcstaticmoments])

            self.area += area
            weighted_area_ += geom.material.E*area
            weighted_staticmoment_ += geom.material.E*staticmoment

        self.nc = (weighted_staticmoment_/weighted_area_)[::-1]  

        # calculate properties regarding normal center

        self.bendingstiffness = np.zeros(3)

        for geom in self._structure.exportgeometry(self.nc):
            inertiamoment = _calcgeom(geom, [calcinertiamoments])[0]

            self.bendingstiffness += inertiamoment*geom.material.E
