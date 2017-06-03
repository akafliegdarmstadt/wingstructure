from .airfoil import Airfoil
from .geometry import Wing
from . import solve_multhopp
import numpy as np
from scipy import interpolate


def wing_lift(alpha: float, wing: Wing, airfoil_db: dict):
    """Calculates wing lift with multhopp method"""

    span_positions = list()
    chord_lenghtes = list()
    alphas = list()

    for section in wing.sections:
        span_positions.append(section.pos.y)
        chord_lenghtes.append(section.chord)
        alpha0 = airfoil_db[section.airfoil].alpha0
        alphas.append(section.alpha-alpha0)

    n = int(round(wing.aspect_ratio())*4-1)
    vs = np.arange(1, n+1)
    thetas= vs*np.pi/(n+1)
    calculation_positions = wing.span_width()/2 * np.cos(np.arange(1, n+1)*np.pi/(n+1))
    calculation_chord_lengths = interpolate.interp1d(span_positions, chord_lenghtes)(np.abs(calculation_positions))
    calculation_alphas = interpolate.interp1d(span_positions, alphas)(np.abs(calculation_positions))

    print(chord_lenghtes, span_positions, calculation_chord_lengths)

    alphas = alpha + calculation_alphas

    result = solve_multhopp(alphas, calculation_positions, calculation_chord_lengths, np.array([2*np.pi]*n),
                            wing.span_width(), wing.aspect_ratio())

    return result


class LiftAnalysis(object):

    def __init__(self, wing: Wing, airfoil_db: dict = None):
        
        if (airfoil_db == None):
            airfoils = set([sec.airfoil for sec in wing.sections])
            airfoil_db = dict()
            
            for airfoil in airfoils:
                airfoil_db[airfoil] = Airfoil()

        self.basic_distribution = self._calculate_basic_distribution(wing, airfoil_db)

        self.aileron_distribution = self._calculate_aileron_distribution(wing)

        self.flap_distribuiton = None

        self.airbrake_distribution = self._calculate_airbrake_distribution(wing)

    def calculate(self, lift):

        return self.basic_distribution*lift

    def _calculate_basic_distribution(self, wing, airfoil_db):

        self.span_positions = list()
        self.chord_lenghtes = list()
        alphas = list()

        for section in wing.sections:
            self.span_positions.append(section.pos.y)
            self.chord_lenghtes.append(section.chord)
            alpha0 = airfoil_db[section.airfoil].alpha0
            alphas.append(section.alpha - alpha0)

        self.n = int(round(wing.aspect_ratio()) * 4 - 1)
        vs = np.arange(1, self.n + 1)
        thetas = vs * np.pi / (self.n + 1)
        self.calculation_positions = wing.span_width() / 2 * np.cos(np.arange(1, self.n + 1) * np.pi / (self.n + 1))
        self.calculation_chord_lengths = interpolate.interp1d(self.span_positions, self.chord_lenghtes)(np.abs(self.calculation_positions))
        calculation_alphas = interpolate.interp1d(self.span_positions, alphas)(np.abs(self.calculation_positions))

        #print(chord_lenghtes, span_positions, calculation_chord_lengths)

        alphas = 1/180*np.pi + calculation_alphas
        
        # TODO: use airfoils lift coefficient slope
        result = solve_multhopp(alphas, self.calculation_positions, self.calculation_chord_lengths, np.array([2 * np.pi] * self.n),
                                wing.span_width(), wing.aspect_ratio())

        return result['c_a_li']/result['C_A']
        
    def _calculate_aileron_distribution(self, wing):
        
        results = {}
        
        for name, flap in wing.flaps.items():
            
            alphas = np.zeros(len(self.calculation_positions))
            
            for ii, span_pos in enumerate(self.calculation_positions):
                
                if flap.depth_at(span_pos) > 0:
                    
                    lambda_k = flap.depth_at(span_pos)
                    
                    eta_k = np.deg2rad(1)
                    
                    k = -2 / np.pi * ( np.sqrt( lambda_k * (1-lambda_k) )  + np.arcsin ( np.sqrt(lambda_k) ) )
            
                    eta_keff = 22.743 * np.arctan( 0.4715 * eta_k )
                    
                    alphas[ii] = -(0.75 * k - 0.25 * lambda_k ) * eta_keff
                    
            result = solve_multhopp(alphas, self.calculation_positions, self.calculation_chord_lengths, np.array( [2*np.pi] *self.n),
                                            wing.span_width(), wing.aspect_ratio())
            
            results[name] = result['c_a_li']
    
        return results
        
    def _calculate_airbrake_distribution(self, wing):
        
        alphas = np.zeros(self.n)
        
        for ii, span_pos in enumerate(self.calculation_positions):
            
            if wing.is_airbrake_pos( span_pos ):
                #TODO: let user choose value
                alphas[ii] = np.deg2rad(-12)
            else:
                alphas[ii] = 0
        
       
        result = solve_multhopp(alphas, self.calculation_positions, self.calculation_chord_lengths, np.array( [2*np.pi] *self.n),
                                wing.span_width(), wing.aspect_ratio())
        
        return result['c_a_li']

