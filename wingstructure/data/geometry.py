# -*- coding: utf-8 -*-
"""
Geometry Module

Contains classes for storing plane geometry
"""

from collections import namedtuple
from sortedcontainers import SortedList
import numpy as np

# data type for 3D-Coordinates
Point = namedtuple('Point', 'x y z')


origin = Point(0, 0, 0)


class Airfoil(object):
    """
    A storage class for Airfoil Coordinates
    """
    def __init__(self, coords: np.array):
        """Creates Airfoil instance from coordinates

        Parameters
        ----------
        coords : np.array
            coordinates of airfoil (dat-file)
        """

        self.coords = coords
        self._n = np.argmin(np.abs(self.coords[:, 0]))

    def interpolate(self, airfoil, beta):
        """interpolation with other airfoil

        Parameters
        ----------
        airfoil: Airfoil
            another airfoil
        beta: float
            interpolation weighting factor

        Returns
        -------
        Airfoil
            interpolated airfoil
        """
        base_upper1 = self._get_upper()
        base_lower1 = self._get_lower()
        
        base_upper2 = airfoil._get_upper()
        base_lower2 = airfoil._get_lower()
        
        upper_coords = np.zeros(base_upper1.shape)
        upper_coords[:, 0] = base_upper1[:, 0]
        upper_coords[:, 1] = np.interp(upper_coords[:, 0], 
                                       base_upper2[::-1, 0],
                                        base_upper2[::-1,1])*(1-beta) +(base_upper1[:,1])*beta
        
        lower_coords = np.zeros(base_lower1.shape)
        lower_coords[:,0] = base_lower1[:,0]
        lower_coords[:,1] = (base_lower1[:,1]*beta + np.interp(lower_coords[:,0], 
                                base_lower2[:,0], base_lower2[:,1])*(1-beta))
        
        return Airfoil(np.vstack((upper_coords, lower_coords)))
        
    def _get_upper(self):
        return self.coords[:self._n,:]
        
    def _get_lower(self):
        return self.coords[self._n:, :]
    
    def plot(self, *args):
        """show airfoil with matplotlib

        Parameters:
        args
            Arguments passed to matplotlib.pyplot.plot function
        """


        from matplotlib import pyplot as plt
        
        plt.plot(self.coords[:,0], self.coords[:,1], *args)


class Section(object):
    """
    A storage class for wing sections
    """
    def __init__(self, pos: Point, chord: float, twist: float = 0.0, airfoil: str = ''):
        """section class constructor
        
        Parameters
        ----------
        pos : Point
            leading edge position
        chord : float
            chord length
        twist : float, optional
            rotation angle around pos (the default is 0.0, which means no rotation)
        airfoil : str, optional
            name of airfoil's section (the default is '', which indicates no airfoil specified)
        
        """

        self.pos = pos
        self.chord = chord
        self.twist = twist
        self.airfoil = airfoil

    @property
    def x(self):
        return self.pos.x

    @property
    def y(self):
        return self.pos.y

    @property
    def z(self):
        return self.pos.z

    def __lt__(self, other) -> bool:
        return self.pos.y < other.pos.y

    def __eq__(self, other):
        return self.pos.y == other.pos.y

    def __repr__(self):
        return 'sec: {{leading edge: {}, chord: {}}}'.format(self.pos,
                                                             self.chord)


class _BaseWing(object):

    def __init__(self, pos=origin, rot=origin, scale=1.0,
                 sections=[], control_surfaces=[]):
        self.sections = SortedList()
        self.pos = pos
        self.rot = rot
        self.root_pos = 0.0

    @classmethod
    def create_from_dict(cls, adict):
        """build _BaseWing from dictionary definition
        """

        pos = Point(**adict['pos'])
        rot = Point(**adict['rot'])

        sections = cls._generate_sections(adict)

        wing = cls(pos=pos, rot=rot, sections=sections)

        return wing
    
    @classmethod
    def _generate_sections(cls, adict: dict)->list:
        """build section list from dictionary
        
        Parameters
        ----------
        adict : dict
            dictionary containing wing definition
        
        Returns
        -------
        list
            list of sections
        """

        sections = []

        for sectiondict in adict['sections']:
            sectiondict['pos'] = Point(**sectiondict['pos'])
            sectiondict['twist'] = np.deg2rad(sectiondict['twist'])
            sections.append(Section(**sectiondict))

        return sections

    def add_section(self, pos: Point, chord: float, twist:float=0., airfoil:str=''):

        tmpsec = Section(pos, chord, twist, airfoil)

        self._add_section(tmpsec)

    def _add_section(self, section: Section):

        self.sections.add(section)

    def set_root_pos(self, span_position: float):
        """Define span position of wing root
        
        Parameters
        ----------
        span_position : float
            wing root span position
        
        """

        self.root_pos = span_position

    @property
    def area(self) -> float:
        """calculate wing area
        
        Returns
        -------
        float
            wing area
        """

        span_positions = [section.pos.y for section in self.sections]
        chord_lengths = [section.chord for section in self.sections]

        area = np.trapz(chord_lengths, span_positions)

        return 2*area # for both sides of wing

    @property
    def aspect_ratio(self) -> float:
        """calculate aspect ratio
        
        Returns
        -------
        float
            aspect ratio
        """

        return self.span**2/self.area

    @property
    def mac(self) -> Section:
        """calculates mean aerodynamic chord.
        
        Returns
        -------
        Section
            mean aerodynamic chord
        """

        wingarea = 0.0
        mac = 0.0
        pos = np.zeros(3)

        for ii in range(len(self.sections)):

            if ii == 0 or self.sections[ii].pos.y < self.root_pos:
                continue

            section = self.sections[ii]
            lastsec = self.sections[ii-1]

            # calculate segments characteristica
            segspan = section.y - lastsec.y
            segarea = (section.chord + lastsec.chord) * segspan / 2
            taper = section.chord / lastsec.chord
            taper1 = taper + 1
            fraction = (taper + taper1) / (3 * taper1)

            segmac = lastsec.chord * (taper ** 2 + taper + 1) / (1.5 * taper1)
            segx = lastsec.x + fraction * (section.x - lastsec.x)
            segy = lastsec.y + fraction * segspan

            # sum values up weighted by segments area
            pos += np.array([segx, segy, 0])* segarea
            mac += segmac * segarea

            wingarea += segarea

        pos /= wingarea
        mac /= wingarea

        return Section(Point(*pos.tolist()),  mac, 0, airfoil=None)

    @property
    def span(self) -> float:
        """Calculate the span width of wing."""
        return 2 * max((section.pos.y for section in self.sections))
          
    def chord_at(self, y: float or np.array) -> float:
        """Calculates the chord depth at given span position"""
        
        y_positions = [section.pos.y for section in self.sections]
        chord_lengths = [section.chord for section in self.sections]

        return np.interp(np.abs(y), y_positions, chord_lengths)

    def le_at(self, y: float) -> np.array:
        """Calculates Leading edge position (x,z) at given span position"""
        
        y_positions = [section.pos.y for section in self.sections]
        x_positions = [section.pos.x for section in self.sections]
        z_positions = [section.pos.z for section in self.sections]
        
        x = np.interp(np.abs(y), y_positions, x_positions)
        z = np.interp(np.abs(y), y_positions, z_positions)
        
        return np.array((x,z))

    def next_airfoil(self, y: float) -> str:
        """Lookup next airfoil at given span position"""

        ys = [section.pos.y for section in self.sections]

        airfoil_index = np.argmin(np.array(ys)-y)

        return self.sections[airfoil_index].airfoil
            
    def airfoils_at(self, y:float) -> tuple:
        """Lookup airfoils around given span position"""
        
        for ii, sec in enumerate(self.sections):
            if sec.pos.y > y:
                break
        
        if ii == 0:
            return {0: self.sections[ii].airfoil}
        
        airfoil1 = self.sections[ii-1].airfoil
        airfoil2 = self.sections[ii].airfoil
        
        if airfoil1 == airfoil2:
            return {0: airfoil1}
            
        else:
            beta = (y-self.sections[ii-1].pos.y)/(self.sections[ii].pos.y-self.sections[ii-1].pos.y)
            return {0: airfoil1, 1: airfoil2, 'beta':beta}

    @property
    def ys(self) -> np.array:
        return np.array([section.pos.y for section in self.sections])
    
    @property
    def chords(self) -> np.array:
        return np.array([section.chord for section in self.sections])

    @property
    def twists(self) -> np.array:
        return np.array([section.twist for section in self.sections])

    @property
    def airfoils(self) -> list:
        return [section.airfoil for section in self.sections]


class ControlSurface(object):
    """Data object for flap definition

    Instances of the class store the span position (start and end) and 
    the chord position of a control surface. The chord position is defined
    relative to the chord length (0.0-1.0) and can either be constant
    or linear. 

    Examples
    --------
    >>> constantflap = ws.ControlSurface(4.5, 6.5, 0.8)
    >>> linearflap = ws.ControlSurface(4.5, 6.5, [0.8, 0.75])
    """
    def __init__(self, span_start, span_end, chord):
        if np.isscalar(chord):
            chord = np.ones(2,1) * chord
        self.y_start = span_start
        self.y_end = span_end
        self.chord_start = chord[0]
        self.chord_end = chord[1]
        
    def chordpos_at(self, span_pos: float):
        """interpolate chord position
        """

        return np.interp(span_pos, [self.y_start, self.y_end], [self.chord_start, self.chord_end],
                         right=0.0, left=0.0)

    def __lt__(self, other) -> bool:
        return self.y_start < other.y_start

    def __eq__(self, other):
        return self.y_start == other.y_start


class Wing(_BaseWing):
    """Storage class for Wing definition

    Instances of the Wing class can be used to store geometry
    information for wing elements. Many useful functions are available
    as well through it's objects. 
    """
    def __init__(self, pos:Point=origin, rot:Point=origin):
        super(Wing, self).__init__(pos, rot)
        self.flaps = dict()
        self.airbrake = None

    @classmethod
    def create_from_dict(cls, adict):

        # create Wing instance
        wing = cls(pos=Point(**adict['pos']),
                   rot=Point(**adict['rot']))

        # generate sections
        wing.sections = cls._generate_sections(adict)

        # add control surfaces
        if 'control-surfaces' in adict.keys() and isinstance(adict['control-surfaces'], list):
            for data in adict['control-surfaces']:
                if data['type'] == 'spoiler':
                    wing.set_airbrake(data['span-start'], data['span-end'])
                elif data['type'] in ('aileron', 'flap', 'flaperon'):
                    wing.set_flap(data['name'], data['span-start'], 
                                  data['span-end'], data['chord-pos'])
                else:
                    raise Exception('not recognized control-surface type: {}'.format(data['type']))

        return wing

    def set_flap(self, name, span_pos_start, span_pos_end, depth):

        if np.isscalar(depth):
            depth = [depth]*2

        flap = ControlSurface(span_pos_start, span_pos_end, depth)
        
        self.flaps[name] = flap
    
    def get_flap_depth(self, span_pos: float)->float:
        for flap in self.flaps.values():
            if flap.y_start <= span_pos <= flap.y_end:
                return flap.depth_at(span_pos)/self.chord_at(span_pos)
        return 0.0
        
    def set_airbrake(self, span_pos_start, span_pos_end):
        if span_pos_end > span_pos_start:
            self.airbrake = {'start': span_pos_start, 'end': span_pos_end}
            
    def is_airbrake_pos(self, span_pos: float) -> bool:
        if self.airbrake is None:
            return False
        else:
            if self.airbrake['start'] <= abs(span_pos) <= self.airbrake['end']:
                return True
            else:
                return False
                
    def plot(self):
        import matplotlib.pyplot as plt
        
        # draw centerline 
        plt.axvline(x=0, linestyle='-.')
        
        # draw sections
        x_positions = []
        y_positions = []
        chord_lengths = []
        
        for section in self.sections:
            x = section.pos.x+self.pos.x
            y = section.pos.y
            chord = section.chord
            
            plt.plot((y, y), (-x, -x-chord), 'r')
            x_positions.append(x)
            y_positions.append(y)
            chord_lengths.append(chord)
        
        y_positions = np.array(y_positions)
        
        # draw leading edge
        plt.plot(y_positions, -1*np.array(x_positions), 'b' )
        # draw trailing edge
        plt.plot(y_positions, -1*np.array(x_positions)-np.array(chord_lengths), 'b')
        
        # draw flaps
        for name, aflap in self.flaps.items():
            y_pos = np.array([aflap.y_start, aflap.y_end])
            y_pos = np.concatenate([y_pos,
                                    y_positions[(y_positions>aflap.y_start) &
                                                (y_positions<aflap.y_end)]])
            y_pos.sort()
            
            chords = np.interp(y_pos, y_positions, chord_lengths)
            offsets = np.interp(y_pos, y_positions, x_positions)
            te = -(chords+offsets)
            depth = (1-aflap.chordpos_at(y_pos))*chords
            
            plt.plot(y_pos[[0,*range(len(y_pos)),-1]], [te[0],*(te+depth),te[-1]], 'g--')

        # draw air brake
        if self.airbrake:
            start = self.airbrake['start']
            end = self.airbrake['end']

            cs = np.interp((start, end), y_positions, chord_lengths)
            xs = np.interp((start, end), y_positions, x_positions)

            plt.plot((start, end), -(xs+cs/2), 'k')

        # format 
        plt.axis('equal')
        plt.axis('off')
        plt.xlim(-1, max(y_positions)+1)
