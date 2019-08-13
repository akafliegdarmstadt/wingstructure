import numpy as np

def mean_momentcoefficient(wing, airfoil_db):
    """calculate mean coefficient of moment for wing
    
    Parameters
    ----------
    wing : Wing
        object describing wing
    airfoil_db : dict
        dictionary containing airfoil data
    """

    try:
        c_m0s = [airfoil_db[sec.airfoil].c_m0 for sec in wing.sections]
    except KeyError:
        raise KeyError('Not all airfoils used in wing are defined in airfoil_db!')

    ys = wing.ys
    cs = wing.chords

    C_m0 = 2 / (wing.area * wing.mac) * np.trapz(c_m0s * cs**2, ys)

    return C_m0