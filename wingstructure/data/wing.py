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
    
    def serializesection(self):
        data = dict(self._asdict())
        data['pos'] = dict(self.pos._asdict())
        return data

    _Section.serialize = serializesection

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
        """Add controlsurface to Wing instance
        
        Parameters
        ----------
        name : str
            identifier for control surface
        pos1 : float
            starting position (spanwise)
        pos2 : float
            end position (spanwise)
        depth1 : float
            start depth or chordwise position (depends on type)
        depth2 : float
            end depth or chordwise position (depends on type)
        cstype : str
            use one of the following type strings: flap, spoiler, airbrake
        """
        self.controlsurfaces[name] = self._ControlSurface(
                pos1, pos2, depth1, depth2, cstype)

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

    def within_control(self, csname, y):
        try:
            cs = self.controlsurfaces[csname]
            return (cs.pos1 <= y) & (y <= cs.pos2)
        except KeyError:
            raise KeyError('{} is not a control surface'.format(csname))
    
    def within_airbrake(self, ys):
        within_ab = np.full_like(ys, False, dtype=bool)
        for cs in self.controlsurfaces.values():
            if cs.cstype in ('airbrake', 'spoiler'):
                within_tmp = (cs.pos1 <= ys) & (cs.pos2 >= ys)
                within_ab = np.where(within_tmp, True, within_ab)
        return within_ab

    def serialize(self):
        data = {
            'pos': {'x': self.x, 'y': self.y, 'z': self.z},
            'symmetric': self.symmetric,
            'sections': [deepcopy(sec.serialize()) for sec in self.sections],
            'controlsurfaces': {name: dict(cs._asdict()) for name, cs in self.controlsurfaces.items()}
        }

        return data

    @classmethod
    def deserialize(cls, adict):
        """Create new Wing instance from dict
        
        Parameters
        ----------
        adict : dict
            dictionary containing wing data
        
        Returns
        -------
        Wing
            instance object
        """
        # create Wing instance
        wing = cls(pos=Point(**adict['pos']))

        # generate sections
        for secdict in adict['sections']:
            secdict_ = deepcopy(secdict)
            secdict_['pos'] = Point(**secdict_['pos'])
            wing.append(**secdict_)

        # add control surfaces
        try:
            for name, csdict in adict['controlsurfaces'].items():
                wing.add_controlsurface(name, **csdict)
        except KeyError:
            pass
        
        return wing
    
    def flatten(self):
        """create flat wing instance
        
        Returns
        -------
        Wing
            a flat wing instance
        """

        newwing = Wing()

        lastsec = None
        lastpos = 0.0

        for section in self.sections:

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
            
            newwing.sections.append(newsec)
            lastsec = section
            lastpos = newsec.pos[1]

        return newwing