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

    def __init__(self, wing: Wing, airfoil_db: dict):

        self.basic_distribution = self._calculate_basic_distribution(wing, airfoil_db)

        self.aileron_distribution = None

        self.flap_distribuiton = None

        self.airbrake_distribution = None

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

        n = int(round(wing.aspect_ratio()) * 4 - 1)
        vs = np.arange(1, n + 1)
        thetas = vs * np.pi / (n + 1)
        self.calculation_positions = wing.span_width() / 2 * np.cos(np.arange(1, n + 1) * np.pi / (n + 1))
        self.calculation_chord_lengths = interpolate.interp1d(self.span_positions, self.chord_lenghtes)(np.abs(self.calculation_positions))
        calculation_alphas = interpolate.interp1d(self.span_positions, alphas)(np.abs(self.calculation_positions))

        #print(chord_lenghtes, span_positions, calculation_chord_lengths)

        alphas = 1/180*np.pi + calculation_alphas

        result = solve_multhopp(alphas, self.calculation_positions, self.calculation_chord_lengths, np.array([2 * np.pi] * n),
                                wing.span_width(), wing.aspect_ratio())

        return result['c_a_li']/result['C_A']
