from .airfoil import Airfoil
from .geometry import Wing
from . import solve_multhopp
import numpy as np
from scipy import interpolate


class LiftAnalysis(object):
    
    def __init__(self, wing: Wing, airfoil_db: dict = None):
        
        if (airfoil_db == None):
            airfoils = set([sec.airfoil for sec in wing.sections])
            airfoil_db = self.airfoil_db = dict()
            
            for airfoil in airfoils:
                self.airfoil_db[airfoil] = Airfoil()

        self.basic_distribution = self._calculate_basic_distribution(wing, airfoil_db)

        self.flaps_distribution, self.flaps_lift = self._calculate_aileron_distribution(wing)

        self.airbrake_distribution, self.airbrake_lift = self._calculate_airbrake_distribution(wing)

    def calculate(self, lift, airbrake = False, flap_deflections = {}):
        """calculate the lift distribution for specific lift an defined flap settings"""
        distributions = np.zeros(self.n)
        
        if airbrake:
            lift -= self.airbrake_lift
            distributions += self.airbrake_distribution
            
        for flap_name in flap_deflections:
            if flap_name in self.flaps_lift:
                factor = 22.743 * np.arctan( 0.04715 * np.array(flap_deflections[flap_name]) )
                lift -= self.flaps_lift[flap_name] * np.sum(factor)
                flap_distribution = self.flaps_distribution[flap_name]
                distributions += flap_distribution*factor[0]+flap_distribution[::-1]*factor[1] 
                
        return self.basic_distribution*lift + distributions

    def _calculate_basic_distribution(self, wing, airfoil_db):
        """calculate basic lift distribution angle 1°"""
        self.span_positions = list()
        self.chord_lenghtes = list()
        alphas = list()

        for section in wing.sections:
            self.span_positions.append(section.pos.y)
            self.chord_lenghtes.append(section.chord)
            alpha0 = airfoil_db[section.airfoil].alpha0
            alphas.append(section.alpha - alpha0)
        
        #TODO: refine calculation grid
        
        self.n = int(round(wing.aspect_ratio()) * 4 - 1)
        vs = np.arange(1, self.n + 1)
        thetas = vs * np.pi / (self.n + 1)
        self.calculation_positions = wing.span_width() / 2 * np.cos(np.arange(1, self.n + 1) * np.pi / (self.n + 1))
        self.calculation_chord_lengths = interpolate.interp1d(self.span_positions, self.chord_lenghtes)(np.abs(self.calculation_positions))
        calculation_alphas = interpolate.interp1d(self.span_positions, alphas)(np.abs(self.calculation_positions))

        alphas = 1/180*np.pi + calculation_alphas
        
        #TODO: solution for fuselages lift
        dcl = np.array([2 * np.pi] * self.n)
        
        #TODO: use airfoils lift coefficient slope
        result = solve_multhopp(alphas, self.calculation_positions, self.calculation_chord_lengths, dcl,
                                wing.span_width(), wing.aspect_ratio())
        
        return result['c_l']/result['C_L']
        
    def _calculate_aileron_distribution(self, wing, angle = 1):
        """calculate aileron lift distribution"""
        distributions = {}
        lift = {}
        
        for name, flap in wing.flaps.items():
            
            alphas = np.zeros(len(self.calculation_positions))
            
            for ii, span_pos in enumerate(self.calculation_positions):
                
                if flap.depth_at(span_pos) > 0 and span_pos > 0:
                    
                    lambda_k = flap.depth_at(span_pos)
                    
                    eta_k = angle
                    
                    k = -2 / np.pi * ( np.sqrt( lambda_k * (1-lambda_k) )  + np.arcsin ( np.sqrt(lambda_k) ) )
            
                    eta_keff = self._calculate_eta_keff( eta_k )
                    
                    alphas[ii] = -(0.75 * k - 0.25 * lambda_k ) * eta_keff
                    
            result = solve_multhopp(alphas, self.calculation_positions, self.calculation_chord_lengths, np.array( [2*np.pi] *self.n),
                                            wing.span_width(), wing.aspect_ratio())
            
            distributions[name] = result['c_l']
            lift[name] = result['C_L']
    
        return distributions, lift
        
    def _calculate_airbrake_distribution(self, wing):
        """calculate airbrake distribution"""
        alphas = np.zeros(self.n)
        
        for ii, span_pos in enumerate(self.calculation_positions):
            
            if wing.is_airbrake_pos( span_pos ):
                #TODO: let user choose value
                alphas[ii] = np.deg2rad(-12)
            else:
                alphas[ii] = 0
        
       
        result = solve_multhopp(alphas, self.calculation_positions, self.calculation_chord_lengths, np.array( [2*np.pi] *self.n),
                                wing.span_width(), wing.aspect_ratio())
        
        return result['c_l'], result['C_L']

    def _calculate_eta_keff(self, eta_k):
    
        return 22.743 * np.arctan( 0.04715 * eta_k )
        
    
class LiftAndMomentAnalysis(LiftAnalysis):
    def __init__(self, wing: Wing, airfoil_db: dict = None):
        super().__init__(wing, airfoil_db)  
        
        cm_0 = [airfoil_db[sec.airfoil].cm0 for sec in wing.sections]
        
        self.moment_basic_distribution = np.interp(np.abs(self.calculation_positions), self.span_positions, cm_0)

        self.moment_aileron_distributions = {}
        
        for name, flap in wing.flaps.items():
            
            moment_temp = np.zeros(self.n)
            
            for ii, span_pos in enumerate(self.calculation_positions):
            
                if flap.y_start <= span_pos <= flap.y_end:
                
                    lambda_k = flap.depth_at(span_pos)
                    
                    moment_temp[ii] = -1.5 * np.sqrt( lambda_k * ( 1 - lambda_k )**3 )
           
            self.moment_aileron_distributions[name] = moment_temp
           
    def calculate(self, lift: float, airbrake = False, flap_deflections = {}):
        """"""
        #TODO: Moment airbrake - Zottel 4.45
        
        lift_distribution = super().calculate(lift, airbrake, flap_deflections)
        
        moment_distribution = np.zeros(self.n)
        
        moment_distribution += self.moment_basic_distribution
        
        for flap_name in flap_deflections:
            
            if flap_name in self.moment_aileron_distributions:
                
                eta_k = np.array( flap_deflections[flap_name] )
                
                eta_keff = self._calculate_eta_keff( eta_k )
                
                temp_distribution = self.moment_aileron_distributions[flap_name]
                
                moment_distribution += (temp_distribution * eta_keff[0] + temp_distribution[::-1] * eta_keff[1])/2
                
        return lift_distribution, moment_distribution

