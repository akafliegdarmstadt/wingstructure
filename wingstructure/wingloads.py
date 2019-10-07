import numpy as np

def calculation_points(flatwing, M):
    """create calculation points for multhopp and structural calculations
    
    Assuming a symmetric wing with first section at span position 0.

    Parameters
    ----------
    flatwing : Wing
        wing representation of flat wing or wing with negligible dihedral
    M : int
        number of points positioned within one half of the wing
    
    Returns
    -------
    array
        span wise calculation positions
    """
    from .aero.multhop import _calc_gridpoints

    b = flatwing.span

    # collect points of interest
    ys_int = []

    for y in flatwing.ys[:-1]:
        ys_int.append(y)

    ys_int.append(b/2 * np.cos(np.pi/(2*(M+1))))

    # sort points of interest

    ys_int = sorted(ys_int)

    # create list for calculation points
    ys_calc = []

    # iterate over segments between points of interest

    for i in range(1, len(ys_int)):
        y1 = ys_int[i-1]
        y2 = ys_int[i]

        seg_b = y2-y1

        seg_M = int(round(M * seg_b/b))

        θs = np.linspace(0, np.pi, seg_M)
        temp_ys = (y1+y2)/2 - seg_b/2 * np.cos(θs)

        ys_calc.extend(temp_ys)

    return [-y for y in ys_calc[:1:-1]] + ys_calc