# -*- coding: utf-8 -*-
"""
Module providing objects for calculation of prandtl's lifting line problem.
"""
import numpy as np
from collections import namedtuple, defaultdict
from .multhop import _calc_gridpoints, Multhop, AirfoilData, _calc_eta_eff


π = np.pi

_calculator_dict = {
    'multhop': Multhop,
}


class LiftAnalysis:
    @classmethod
    def generate(cls, wing, airfoil_db=defaultdict(AirfoilData), M=None, method='multhop'):
        """Build a LiftAnalysis object
        
        Parameters
        ----------
        wing : Wing
            a wing object
        airfoil_db : [type], optional
            , by default defaultdict(AirfoilData)
        M : int, optional
            number of grid points for calculation, must be uneven, by default None
        method : str, optional
            calculation method, by default 'multhop' (only option atm)
        
        Returns
        -------
        LiftAnalysis
            lift analysis object
        
        Raises
        ------
        ValueError
            M has to be uneven if chosen!
        NotImplementedError
            is raised when chosen cacluation method is available
        """
        
        if M is not None and M%2 != 1:
            raise ValueError('Number of grid points (M) has to be uneven!')

        ys = _calc_gridpoints(wing, M)

        ys[len(ys)//2] = 0.0

        try:
            calculator = _calculator_dict[method](wing, ys, airfoil_db)
        except:
            raise NotImplementedError('{} is not implemented yet!'.format(method))

        analysis = cls()

        analysis.ys = ys
        analysis._base = calculator.baselift()
        analysis._airbrake = calculator.airbrakelift()
        analysis._aoa = calculator.aoa(1)
        analysis._control_surfaces = {name : calculator._controlsurfacelift_n(name) \
                                                     for name in wing.controlsurfaces.keys()}
        analysis.chords = calculator.chords

        return analysis

    def __init__(self):
        self.ys = None
        self.chords = None
        self._base = None
        self._airbrake = None
        self._aoa = None
        self._control_surfaces = {}
    
    def calculate(self, C_L, controls:dict={}, airbrake:bool=False):
        
        α, res = self.__call__(C_L, 'C_L', controls, airbrake)

        return np.rad2deg(α), res.c_ls, res.C_Di

    def _calc_controlsurface(self, name, deflections):
        # control surface lift distribution is proportional
        # to effective deflection not to deflection itself
        fac = lambda η: _calc_eta_eff(np.radians(η))
        controllift = self._control_surfaces[name]
        return controllift * fac(deflections[0]) \
                + controllift.flip() * fac(deflections[1])

    def __call__(self, target=0.0, target_type='C_L', controls={}, airbrake=False):
        
        from .multhop import MulthopResult

        C_L = -target if target_type=='C_L' else 0.0

        def zly(): return np.zeros_like(self.ys)

        res = MulthopResult(self.ys, zly(), zly(), C_L, 0.0)

        # collect contributions to lift distribution
        res += self._base

        if airbrake: 
            res += self._airbrake

        for cs in controls.items():
            res += self._calc_controlsurface(*cs)

        # choose angle of attack
        if target_type == 'C_L':
            α = -res.C_L / self._aoa.C_L
        else:
            α = target

        res += α * self._aoa

        res.C_L += target if target_type=='C_L' else 0.0

        return α, res


def calculate(wing, target=0.0, target_type='C_L', controls={}, airbrake=False, M=None, 
        method='multhop', airfoil_db:dict=defaultdict(AirfoilData), calc_cmx=False):
        
    la = LiftAnalysis.generate(wing, airfoil_db, M, method)

    α, res = la(target, target_type, controls, airbrake)

    additional = {}

    if calc_cmx:
        # Moment coefficient in flight direction
        A = wing.area
        b = wing.span
        C_Mx = np.trapz(la.ys*res.c_ls*la.chords, la.ys)/ (A*b)

        additional['C_Mx'] = C_Mx
    
    return {'c_ls': res.c_ls, 'a_is': res.α_is, 'C_L': res.C_L, 'alpha': np.rad2deg(α), **additional}
