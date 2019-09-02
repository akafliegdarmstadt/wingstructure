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
    
    def calculate(self, C_L, controls:dict={}, airbrake:bool=False, airfoil_db:dict=None, options:dict=None):
        
        c_ls = np.zeros_like(self.ys)
        C_Di = 0.0

        # collect contributions to lift distribution
        contributions = [self._base]
        if airbrake: contributions.append(self._airbrake)
        
        contributions += [self._calc_controlsurface(*cs) for cs in controls.items()]
        
        for contrib in contributions:
            C_Di += contrib.C_Di
            c_ls += contrib.c_ls
            C_L -= contrib.C_L

        # choose angle of attack to reach C_L
        α = C_L / self._aoa.C_L

        c_ls += α * self._aoa.c_ls
        C_Di += abs(α) * self._aoa.C_Di

        return np.rad2deg(α), c_ls, C_Di

    def _calc_controlsurface(self, name, deflections):
        # control surface lift distribution is proportional
        # to effective deflection not to deflection itself
        fac = lambda η: _calc_eta_eff(np.radians(η))
        controllift = self._control_surfaces[name]
        return controllift * fac(deflections[0]) \
                + controllift.flip() * fac(deflections[1])

def calculate(wing, alpha=None, C_L=None, controls={}, airbrake=False, M=None, method='multhop', airfoil_db:dict=defaultdict(AirfoilData),
    calc_cmx=False):
        
        ys = _calc_gridpoints(wing, M)

        try:
            calculator = _calculator_dict[method](wing, ys, airfoil_db)
        except:
            raise Exception('{} is not implemented yet!'.format(method))

        base_res = calculator.baselift()

        if airbrake:
            base_res += calculator.airbrakelift()

        for name, (η_r, η_l) in controls.items():
            base_res += calculator.controlsurfacelift(name, η_r)
            base_res += calculator.controlsurfacelift(name, η_l).flip()

        if alpha is not None:
            base_res += calculator.aoa(alpha)

            if C_L is not None:
                raise Warning('C_L value is ignored because aoa is defined!')

        elif C_L is not None:
            res_aoa = calculator.aoa(1)

            C_L -= base_res.C_L

            alpha = C_L/res_aoa.C_L

            base_res += alpha * res_aoa

        additional = {}

        if calc_cmx:
            # Moment coefficient in flight direction
            A = wing.area
            b = wing.span
            C_Mx = np.trapz(ys*base_res.c_ls*calculator.chords, ys)/ (A*b)

            additional['C_Mx'] = C_Mx
        
        return {'c_ls': base_res.c_ls, 'a_is': base_res.α_is, 'C_L': base_res.C_L, 'alpha': np.rad2deg(alpha), **additional}
