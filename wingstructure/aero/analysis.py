from collections import defaultdict
from warnings import warn

import numpy as np
from scipy import interpolate

from ..data import geometry
from . import liftingline as ll

class AirfoilData(object):
    def __init__(self, alpha0: float = 0, dif_ca_alpha: float = 2*np.pi, cm0: float = 0):
        self.alpha0 = alpha0
        self.dif_ca_alpha = dif_ca_alpha
        self.cm0 = cm0


class LiftAnalysis(object):
    """
    Class for simplified handling of Multhopp

    Parameters
    ----------
    wing : geometry.Wing
        Definition of wing to be analysed.
    airfoil_db : dict
        Database of airfoils used in wing.

    """
    def __init__(self, wing: geometry.Wing, airfoil_db: dict = None):

        # use default airfoil if no database is set
        if airfoil_db is None:
            warn('No airfoil database defined, using default airfoil.')
            airfoil_db = defaultdict(AirfoilData)

        # calculate grid points
        θs, self.calc_ys = ll._calcgridpoints(wing.span, wing.aspect_ratio)

        # number of gridpoints
        M = len(θs)
        self.n = M

        # interpolate 
        self.calc_chords = np.interp(np.abs(self.calc_ys), wing.ys, wing.chords)
        calc_dcls = np.array([2*np.pi]*len(self.calc_ys))

        # shorter multhopp call
        def multhopp_(αs):
            res = ll.multhopp(αs, self.calc_chords, self.calc_ys, calc_dcls, M=M,
                              mode='combined', interp=False)
            return res.c_ls, res.C_L, res.C_Wi

        # calculate distributions

        ## base distribution
        def get_α0(airfoil, airfoil_db):
            return np.radians(airfoil_db[airfoil].alpha0)
        
        aerotwists = np.array([get_α0(airfoil, airfoil_db) for airfoil in wing.airfoils])

        αs_base1 = np.interp(np.abs(self.calc_ys), wing.ys, wing.twists) # geometric twist
        αs_base2 = np.interp(np.abs(self.calc_ys), wing.ys, aerotwists) # aerodynamic twist

        self.base_liftdist, self.base_lift, self.base_drag = multhopp_(αs_base1-αs_base2)

        ## flap distributions
        self.flap_liftdist = {}
        self.flap_lift = {}
        self.flap_drag = {}
        for name, flap in wing.flaps.items():
            αs_flap = self._calculate_flap_α(flap)
            self.flap_liftdist[name], self.flap_lift[name], self.flap_drag[name] = multhopp_(αs_flap)

        ## air brake distribution
        if wing.airbrake:
            from numpy import abs
            α_ab = np.zeros(M)
            α_ab[ (abs(self.calc_ys)>wing.airbrake['start']) & (abs(self.calc_ys)<wing.airbrake['end'])] = np.radians(-12)
            self.airbrake_distribution, self.airbrake_lift, self.airbrake_drag = multhopp_(α_ab)

        ## lift due to angle of attack
        #print('\n\naoa')
        α_aoa = np.ones(M)
        self.aoa_c_ls, aoa_C_L, self.aoa_C_Di = multhopp_(α_aoa)
        self.aoa_c_ls /= aoa_C_L
        self.aoa_α = np.rad2deg(1/aoa_C_L)

    def _calculate_flap_α(self, flap, angle=np.radians(1)):

        alphas = np.zeros(len(self.calc_ys))
            
        for ii, span_pos in enumerate(self.calc_ys):
            
            if flap.length(span_pos) > 0.0 and span_pos > 0:
                
                lambda_k = 1-flap.chordpos_at(span_pos)
                
                eta_k = angle
                
                k = -2 / np.pi * (np.sqrt(lambda_k * (1-lambda_k)) + np.arcsin(np.sqrt(lambda_k)))
        
                eta_keff = self._calculate_eta_keff(eta_k)
                
                alphas[ii] = -(0.75 * k - 0.25 * lambda_k ) * eta_keff

        return alphas

    def _calculate(self, C_L: float, air_brake=False, flap_deflections={})->tuple:

        # create float variable for drag
        C_Di = self.base_drag

        # create empty array for lift distribution
        distribution = np.copy(self.base_liftdist)

        # subtract lift resulting from aerodynamic and geometric twist
        # from demanded lift coefficient
        C_L -= float(self.base_lift)

        # take air brake into account
        if air_brake:
            C_L -= self.airbrake_lift
            C_Di += self.airbrake_drag
            distribution += self.airbrake_distribution

        # iterate over given flap deflections
        for flap_name in flap_deflections:
            # when flap exists take flap's impact on lift into account
            if flap_name in self.flap_lift:
                factor = self._calculate_eta_keff(np.array(flap_deflections[flap_name]))
                C_L -= self.flap_lift[flap_name] * factor
                C_Di += self.flap_drag[flap_name] * np.abs(factor)
                flap_distribution = self.flap_liftdist[flap_name]
                distribution += flap_distribution * factor[0] + flap_distribution[::-1] * factor[1]
            # otherwise warn user
            else:
                from warnings import warn
                warn('flap {} does not exist'.format(flap_name))

        return np.mean(C_L), self.aoa_c_ls * np.mean(C_L) + distribution, C_Di + self.aoa_C_Di*np.mean(C_L)

    def calculate(self, C_L: float, air_brake=False, flap_deflections={}, return_C_Di=False)->tuple:
        """Calculates lift distribution with defined control surface setting
      
        Parameters
        ----------
        C_L : float
            wing lift coefficient
        air_brake : bool, optional
            Is air brake active (the default is False, which mean not active)
        flap_deflections : dict, optional
            deflections of control surfaces as dictionary {cs-name: (deflection-left, deflection-right)}
        return_C_Di : bool, optional
            Is induced drag also returned? (the default is False, which means it will not)
        
        Returns
        -------
        tuple
            (angle of attack, local lift coefficients, [induced drag])
        """

        C_L_, c_ls, C_Di = self._calculate(C_L, air_brake, flap_deflections)

        if not return_C_Di:
            return self.aoa_α * C_L_, c_ls
        else:
            return self.aoa_α * C_L_, c_ls, C_Di
    
    @staticmethod
    def _calculate_eta_keff(eta_k: float or np.ndarray) -> float:
        return 22.743 * np.arctan(0.04715 * eta_k)
        
    
class LiftAndMomentAnalysis(LiftAnalysis):
    """
    Extended Analysis, calculates lift and aerodynamic moments

    Parameters
    ----------
    wing : geometry.Wing
        Definition of wing to analyse
    airfoil_db : dict
        Database of airfoil coefficients {airfoilname: AirfoilData}
    """

    def __init__(self, wing: geometry.Wing, airfoil_db: dict = {}):

        super().__init__(wing, airfoil_db)  
        
        cm_0 = [airfoil_db[airfoil].cm0 for airfoil in wing.airfoils]
        
        self.moment_basic_distribution = np.interp(np.abs(self.calc_ys), wing.ys, cm_0)

        self.moment_aileron_distributions = {}

        # TODO: Fix Flap Moments
        for name, flap in wing.flaps.items():
            
            moment_temp = np.zeros(self.n)
            
            for ii, span_pos in enumerate(self.calc_ys):
            
                if flap.y_start <= span_pos <= flap.y_end:
                
                    lambda_k = flap.chordpos_at(span_pos)
                    
                    moment_temp[ii] = -1.5 * np.sqrt(lambda_k * (1 - lambda_k)**3)
           
            self.moment_aileron_distributions[name] = moment_temp
           
    def calculate(self, lift: float, air_brake=False, flap_deflections={}):
        """Calculates lift and moment distribution  
        
        Parameters
        ----------
        lift : float
            wing lift coefficient
        air_brake : bool, optional
            Is airbrake active?
        flap_deflections : dict, optional
            Dictionary containing flap deflections, e.g. {'flap1':[np.radians(5),np.radians(-5)]}
        
        Returns
        -------
        tuple
            (wing angle of attack, distribution of local lift coefficients, 
                distribution of local moment coefficients)
        """

        # TODO: Moment airbrake - Zottel 4.45
        
        alpha, lift_distribution = super().calculate(lift, air_brake, flap_deflections)
        
        moment_distribution = np.zeros(self.n)
        
        moment_distribution += self.moment_basic_distribution
        
        for flap_name in flap_deflections:
            
            if flap_name in self.moment_aileron_distributions:
                
                eta_k = np.array(flap_deflections[flap_name])
                
                eta_keff = self._calculate_eta_keff(eta_k)
                
                temp_distribution = self.moment_aileron_distributions[flap_name]
                
                moment_distribution += (temp_distribution * eta_keff[0] + temp_distribution[::-1] * eta_keff[1])/2
                
        return alpha, lift_distribution, moment_distribution
