from collections import namedtuple
from copy import deepcopy

import numpy as np


# data type for 3D-Coordinates
Point = namedtuple('Point', 'x y z')


class _Wing:
    """A data structue for multi trapez wing definitions.

    Not inteded for usage. Have a look at the Wing class.
    """

    _Section = namedtuple('Section', ['pos', 'chord', 'twist', 'airfoil'])

    def __init__(self,pos=(0.0, 0.0, 0.0), symmetric=True):
        self.x, self.y, self.z = pos
        self.symmetric = symmetric
        self.sections = []
    
    def append(self, pos=(0.0, 0.0, 0.0), chord=1.0, twist=0.0, airfoil=''):
        self.sections.append(
            self._Section(Point(*pos), chord, twist, airfoil)
        )

    def get_mac(self):
        """Calculate mean aerodynamic chord.
        
        Returns
        -------
        pos: arraylike
           leading edge position of mean aerodynamic chord
        mac: float
           mac length
        
        Notes
        -----
        Implements formulas reported in http://dx.doi.org/10.1063/1.4951901
        """
        
        pos = np.zeros(3)
        area = 0.0
        mac = 0.0

        lastsec = None

        for sec in self.sections:
            if lastsec is None:
                lastsec = sec
                continue
            
            # short aliases for usage in formulas
            x1, x2 = lastsec.pos.x, sec.pos.x
            y1, y2 = lastsec.pos.y, sec.pos.y
            c1, c2 = lastsec.chord, sec.chord

            # segment properties
            S = (c1+c2)/2 * (y2-y1)
            λ = c2 / c1

            segmac = 2/3 * c1 * (λ**2 + λ + 1) / (λ + 1)         
            segx = x1 + (x2-x1) * (1+2*λ)/(3+3*λ)
            segy = y1 + (y2-y1) * (1+2*λ)/(3+3*λ)

            # sum up values weighted by segment area
            pos += np.array([segx, segy, 0]) * S
            mac += segmac * S

            area += S

            lastsec = sec
        
        pos /= area
        mac /= area

        return pos, mac

    @property
    def span(self):
        """Get span of wing."""
        if self.symmetric:
            return 2*max((sec.pos.y for sec in self.sections))

    @property
    def area(self):
        """Get wing area."""

        span_positions = [sec.pos.y for sec in self.sections]
        chord_lengths = [sec.chord for sec in self.sections]

        area = np.trapz(chord_lengths, span_positions)

        return 2*area if self.symmetric else area

    @property
    def aspectratio(self):
        """Get aspect ratio."""
        return self.span**2/self.area
    
    @property
    def mac(self):
        """Get mac length"""
        return self.get_mac()[1] 


class Wing(_Wing):
    """A object representing lift generating airplane parts.
    
    Parameters
    ----------
    pos: float
       coordinate system offset
    rot: float
       
    symmetric: bool
       wing symmetry regarding y=0.0
    """

    _ControlSurface = namedtuple('ControlSurface', 
            ['pos1', 'pos2', 'depth1', 'depth2', 'cstype'])

    def __init__(self, pos=(0.0, 0.0, 0.0), symmetric=True):
        super().__init__(pos, symmetric)
        self.controlsurfaces = {}

    def add_controlsurface(self, name, pos1, pos2, depth1, depth2, cstype):
        self.controlsurfaces[name] = self._ControlSurface(
                pos1, pos2, depth1, depth2, cstype)

    def _repr_svg_(self):
        pass

    @property
    def chords(self):
        return np.array([sec.chord for sec in self.sections])
    
    @property
    def ys(self):
        return np.array([sec.pos.y for sec in self.sections])

    @property
    def twists(self):
        return np.array([sec.twist for sec in self.sections])

    @property
    def airfoils(self):
        return np.array([sec.airfoil for sec in self.sections])

    def within_airbrake(self, y):
        within_ab = np.full_like(y, False)
        for cs in self.controlsurfaces.values():
            if cs.cstype in ('airbrake', 'spoiler'):
                within_tmp = (cs.pos1 <= y) & (cs.pos2 >= y)
                np.where(within_tmp, True, within_ab)
        return within_ab

    def within_control(self, csname, y):
        try:
            cs = self.controlsurfaces[csname]
            return (cs.pos1 <= y) & (y <= cs.pos2)
        except KeyError:
            return np.full_like(y, False)

    @classmethod
    def from_dict(cls, adict):
        
        # create Wing instance
        wing = cls(pos=Point(**adict['pos']))

        # generate sections
        for secdict in adict['sections']:
            secdict_ = deepcopy(secdict)
            secdict_['pos'] = Point(**secdict_['pos'])
            wing.append(**secdict_)

        # add control surfaces
        try:
            for name, csdict in adict['control-surfaces'].items():
                wing.add_controlsurface(name, **csdict)
        except KeyError:
            pass
        
        return wing

