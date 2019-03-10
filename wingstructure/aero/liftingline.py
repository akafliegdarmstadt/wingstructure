# -*- coding: utf-8 -*-
"""
Module providing objects for calculation of prandtl's lifting line problem.
"""
import numpy as np
from collections import namedtuple, defaultdict
from .multhop import _calc_gridpoints, Multhop


π = np.pi


class LiftAnalysis:
    @classmethod
    def generate(cls, wing, airfoil_db=defaultdict(AirfoilData), method='multhop'):
        
        ys = _calc_gridpoints(wing, 81)

        if method == 'multhop':
            calculator = Multhop(wing, ys, airfoil_db)

        analysis = cls()

        analysis.ys = ys
        analysis._base = calculator.baselift()
        analysis._airbrake = calculator.airbrakelift()
        analysis._coontrol_surfaces = [calculator.controlsurfacelift(name, 1.0) \
                                                     for name in wing.controls.keys()]

        return analysis


    def __init__(self):
        self.ys = None
        self._base = None
        self._airbrake = None
        self._coontrol_surfaces = []
    
    
    def calculate(self, C_L, controls:dict={}, airbrake:bool=False, airfoil_db:dict=None, options:dict=None):
        
        c_ls = np.zeros_like(self.ys)
        C_Di = 0.0

        # collect contributions to lift distribution
        contributions = [self._base, self._airbrake]
        contributions += [self._calc_controlsurface(*cs) for cs in controls]

        for contrib in contributions:
            C_Di += contrib.C_Di
            c_ls += contrib.c_ls
            C_L -= contrib.C_L

    def _calc_controlsurface(self, name, deflections):
        controllift = self._coontrol_surfaces[name]
        return controllift * deflections[0] + controllift.flip() * deflections[1]
