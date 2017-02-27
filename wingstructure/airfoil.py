import numpy as np


class Airfoil(object):
    def __init__(self, alpha0: float = 0, dif_ca_alpha: float = 2*np.pi, cm0: float = 0):
        self.alpha0 = alpha0
        self.dif_ca_alpha = dif_ca_alpha
        self.cm0 = cm0
