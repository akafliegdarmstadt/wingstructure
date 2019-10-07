import numpy as np

from wingstructure.data import Wing


class FlatWing(Wing):
    """A class representing the flattend version of a wing

    Can be useful in calculation of aerodynamic loads with methods not suporting
    dihedral. Flat means parallel to the x-y-plane.
    
    Parameters
    ----------
    wing : Wing
        the wing to be flattend
    """
    
    def __init__(self, wing):
        """create flat wing instance
        
        Returns
        -------
        Wing
            a flat wing instance
        """

        super().__init__((wing.x, wing.y, wing.z))

        self.basewing = wing

        lastsec = None
        lastpos = 0.0

        for section in wing.sections:

            curpos = section.pos

            if lastsec is not None:
                # calculate section distance in y-z-plane
                pos1 = np.array(section.pos)
                pos2 = np.array(lastsec.pos)

                # ignore x -> 2D vector
                distvec = (pos1-pos2)[(1,2),] 
                
                dist = np.linalg.norm(distvec)

                newpos = curpos._replace(y=lastpos+dist, z=0.0)
                newsec = section._replace(pos=newpos)
                
            else:
                # set z to zero
                newpos = curpos._replace(z=0.0)
                newsec = section._replace(pos=newpos)
            
            self.sections.append(newsec)
            lastsec = section
            lastpos = newsec.pos[1]

        for cs_name, controlsurface in wing.controlsurfaces.items():
            pos1, pos2 = controlsurface.pos1, controlsurface.pos2

            newpos1, newpos2 = np.interp((pos1, pos2), wing.ys, self.ys)

            self.controlsurfaces[cs_name] = controlsurface._replace(pos1=newpos1, pos2=newpos2)

def calculation_points(flatwing, M):
    """create calculation points for multhopp and structural calculations
    
    Assuming a symmetric wing with first section at span position 0.

    Parameters
    ----------
    flatwing : FlatWing
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

