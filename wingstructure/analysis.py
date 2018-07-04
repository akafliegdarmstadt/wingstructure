from collections import defaultdict
from warnings import warn

import numpy as np
from scipy import interpolate

from .airfoil import Airfoil
from .geometry import BaseWing
from .liftingline import _calculate_liftcoefficients, _solve_multhopp, multhopp


class LiftAnalysis(object):
    """
    Class for simplified handling of Multhopp
    """
    def __init__(self, wing: BaseWing, airfoil_db: dict = None):
        """
        Initialise Analysis object, calculate distributions
        """

        # use default airfoil if no database is set
        if airfoil_db is None:
            warn('No airfoil database defined, using default airfoil.')
            airfoil_db = defaultdict(Airfoil)

        # calculate the distributions
        # 1) lift distribution of wing without any flaps set
        self._calculate_base_distribution(wing, airfoil_db)
        # 2) lift distribution of individual flaps
        self.flaps_distribution, self.flaps_lift = self._calculate_aileron_distribution(wing)
        # 3) lift distribution of airbrakes
        self.airbrake_distribution, self.airbrake_lift = self._calculate_airbrake_distribution(wing)
        # 4) lift distribution due to angle of attack
        self._calculate_aoa_distribution(wing)

    def _calculate(self, C_L: float, air_brake=False, flap_deflections={})->tuple:
        """
        Calculates lift distribution and resultant coefficients for defined flap settings and lift
        :param c_l: demanded lift coefficient for whole wing
        :param air_brake: enables or disables airbrake [True/False]
        :param flap_deflections: dictionary containing flap deflections, e.g. {'flap1':[np.radians(5),np.radians(-5)]}
        :return: (angle of attack, lift distribution)
        """

        # create empty array for lift distribution
        distribution = np.copy(self.base_distribution)

        C_L -= float(self.base_lift)

        # take air brake impact into account
        if air_brake:
            C_L -= self.airbrake_lift
            distribution += self.airbrake_distribution

        # iterate over given flap deflections
        for flap_name in flap_deflections:
            # when flap exists take flap's impact on lift into account
            if flap_name in self.flaps_lift:
                factor = self._calculate_eta_keff(np.array(flap_deflections[flap_name]))
                C_L -= self.flaps_lift[flap_name] * factor
                flap_distribution = self.flaps_distribution[flap_name]
                distribution += flap_distribution * factor[0] + flap_distribution[::-1] * factor[1]
            # otherwise warn user
            else:
                from warnings import warn
                warn('flap {} does not exist'.format(flap_name))

        return C_L, self.aoa_distribution * np.mean(C_L) + distribution

    def calculate(self, C_L: float, air_brake=False, flap_deflections={})->tuple:
        """
        Calculates the lift distribution for specific lift an defined flap settings
        :param c_l: demanded lift coefficient for whole wing
        :param air_brake: enables or disables airbrake [True/False]
        :param flap_deflections: dictionary containing flap deflections, e.g. {'flap1':[np.radians(5),np.radians(-5)]}
        :return: (angle of attack, lift distribution)
        """

        C_L_, c_ls = self._calculate(C_L, air_brake, flap_deflections)
                
        return self.aoa_base_angle * np.mean(C_L_), c_ls

    def calculate_resultant(self, C_L, air_brake=False, flap_deflections={})->tuple:
        """
        Determines resulting lift coefficients for left and right wing half
        :param c_l: demanded lift coefficient for whole wing
        :param air_brake: enables or disables airbrake [True/False]
        :param flap_deflections: dictionary containing flap deflections, e.g. {'flap1':[np.radians(5),np.radians(-5)]}
        :return: (lift coefficient left, lift coefficient right)
        """

        C_L_, c_l_dis = self._calculate(C_L, air_brake, flap_deflections)

        # calculate lever for lift

        y_lever = np.zeros(2)

        # shorter names for used values
        c_i = self.calc_chords
        y_i = self.calc_ys

        # span width of calculation wing segments
        diff_y_i = -np.diff(y_i)

        # area of calculated wing segments
        area_i = (c_i[:-1] + c_i[1:]) / 2 * diff_y_i

        # redistribute area, target-> area relevant for given lift coeffients
        area_ip = np.pad(area_i, 1, 'edge')

        area_im = (area_ip[1:] + area_ip[:-1]) / 2

        # calculate leaver arms for local lift
        y_l = np.copy(y_i)
        y_l[0] = np.mean(y_l[:2])
        y_l[-1] = np.mean(y_l[-2:])

        mid = len(y_l)//2

        y_l_left = y_l[mid:]
        y_l_left[0] = np.mean(y_l_left[:2])
        y_l_right = y_l[:mid+1]
        y_l_right[-1] = np.mean(y_l_right[-2:])

        y_lever[0] = sum(c_l_dis[mid:] * area_im[mid:] * y_l_left) /\
                    sum(c_l_dis[mid:] * area_im[mid:])
        y_lever[1] = sum(c_l_dis[:mid+1] * area_im[:mid+1] * y_l_right) /\
                    sum(c_l_dis[:mid+1] * area_im[:mid+1])

        # add base and flap lift
        C_L_ += C_L - np.mean(C_L_)

        return (*y_lever, *C_L_)

    def _calculate_base_distribution(self, wing, airfoil_db):
        """
        calculate basic lift distribution for aoa 1°
        """

        alphas = np.radians(wing.alphas)

        for ii, section in enumerate(wing.sections):
            alpha0 = airfoil_db[section.airfoil].alpha0
            alphas[ii] -= np.radians(alpha0)

        # TODO: refine calculation grid
        self.n = int(round(wing.aspect_ratio) * 4 - 1)

        self.thetas = np.linspace(np.pi/(self.n+1), self.n/(self.n+1)*np.pi, self.n)
        self.calc_ys = wing.span / 2 * np.cos(self.thetas)
        self.calc_chords = np.interp(np.abs(self.calc_ys), wing.ys, wing.chords)
        calc_alphas = np.interp(np.abs(self.calc_ys), wing.ys, alphas)

        # no lift for fuselage
        calc_alphas[np.abs(self.calc_ys) < wing.root_pos] = 0.0

        # TODO: use airfoils lift coefficient slope
        self.dcl = np.array([2 * np.pi] * self.n)

        # use multhopp method for calculation
        result = _solve_multhopp(calc_alphas, self.thetas, self.calc_chords, wing.span, self.dcl)
        
        c_ls, C_L = _calculate_liftcoefficients(self.thetas, result, self.calc_chords,
                                         wing.aspect_ratio, wing.span, self.n)

        self.base_lift = C_L
        self.base_distribution = c_ls
        
    def _calculate_aoa_distribution(self, wing):
        """
        calculates lift distribution due to angle of attack
        """

        res = multhopp([1]*len(wing.chords), wing.chords, wing.ys)

        self.aoa_distribution = res.c_ls/res.C_L
        self.aoa_base_angle = np.rad2deg(1/res.C_L)

    def _calculate_aileron_distribution(self, wing, angle=np.radians(1)):
        """
        calculates aileron lift distributions and lift coefficients
        """

        distributions = {}
        lift = {}
        
        for name, flap in wing.flaps.items():
            
            alphas = np.zeros(len(self.calc_ys))
            
            for ii, span_pos in enumerate(self.calc_ys):
                
                if flap.depth_at(span_pos) > 0 and span_pos > 0:
                    
                    lambda_k = flap.depth_at(span_pos)
                    
                    eta_k = angle
                    
                    k = -2 / np.pi * (np.sqrt(lambda_k * (1-lambda_k)) + np.arcsin(np.sqrt(lambda_k)))
            
                    eta_keff = self._calculate_eta_keff(eta_k)
                    
                    alphas[ii] = -(0.75 * k - 0.25 * lambda_k ) * eta_keff
                    
            result = _solve_multhopp(alphas, self.thetas, self.calc_chords, wing.span, self.dcl)
            
            c_ls, C_L = _calculate_liftcoefficients(self.thetas, result, self.calc_chords,
                                         wing.aspect_ratio, wing.span, self.n)

            distributions[name] = c_ls/eta_keff
            lift[name] = C_L
    
        return distributions, lift
        
    def _calculate_airbrake_distribution(self, wing):
        """
        calculates airbrake lift distribution and lift coefficient
        """
        alphas = np.zeros(self.n)
        
        for ii, span_pos in enumerate(self.calc_ys):
            
            if wing.is_airbrake_pos( span_pos ):
                #TODO: let user choose value
                alphas[ii] = np.radians(-12)
            else:
                alphas[ii] = 0

        result = _solve_multhopp(alphas, self.thetas, self.calc_chords, wing.span, self.dcl)
        
        return _calculate_liftcoefficients(self.thetas, result, self.calc_chords,
                                         wing.aspect_ratio, wing.span, self.n)

    @staticmethod
    def _calculate_eta_keff(eta_k: float or np.ndarray) -> float:
        return 22.743 * np.arctan(0.04715 * eta_k)
        
    
class LiftAndMomentAnalysis(LiftAnalysis):
    """
    Extended Analysis, calculates lift and moments
    """
    def __init__(self, wing: BaseWing, airfoil_db: dict = None):

        super().__init__(wing, airfoil_db)  
        
        cm_0 = [airfoil_db[airfoil].cm0 for airfoil in wing.airfoils]
        
        self.moment_basic_distribution = np.interp(np.abs(self.calc_ys), wing.ys, cm_0)

        self.moment_aileron_distributions = {}

        # TODO: Fix Flap Moments
        for name, flap in wing.flaps.items():
            
            moment_temp = np.zeros(self.n)
            
            for ii, span_pos in enumerate(self.calc_ys):
            
                if flap.y_start <= span_pos <= flap.y_end:
                
                    lambda_k = flap.depth_at(span_pos)
                    
                    moment_temp[ii] = -1.5 * np.sqrt(lambda_k * (1 - lambda_k)**3)
           
            self.moment_aileron_distributions[name] = moment_temp
           
    def calculate(self, lift: float, air_brake=False, flap_deflections={}):
        """
        calculates lift and moment distribution
        :param lift: demanded lift coefficient
        :param air_brake: enable or disable [True/False]
        :param flap_deflections: flap_deflections: dictionary containing flap deflections, e.g. {'flap1':[np.radians(5),np.radians(-5)]}
        :return: (aoa, lift-distribution, moment-distribution)
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
