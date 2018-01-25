# -*- coding: utf-8 -*-
"""
Geometry Module

Contains classes for storing plane geometry
"""

from collections import namedtuple
from sortedcontainers import SortedList
import numpy as np


Point = namedtuple('Point','x y z')


class Airfoil(object):
    """
    A storage class for Airfoil Coordinates
    """
    def __init__(self, coords: np.array):
        self.coords = coords
        self._n = np.argmin(np.abs(self.coords[:,0]))
        
    def interpolate(self, airfoil, beta):
        """
        Calculates interpolation with other airfoil
        :param airfoil: another Airfoil object
        :param beta: interpolation weighting factor
        :return: interpolated Airfoil object
        """
        base_upper1 = self._get_upper()
        base_lower1 = self._get_lower()
        
        base_upper2 = airfoil._get_upper()
        base_lower2 = airfoil._get_lower()
        
        upper_coords = np.zeros(base_upper1.shape)
        upper_coords[:,0] = base_upper1[:,0]
        upper_coords[:,1] = np.interp(upper_coords[:,0], base_upper2[::-1,0], base_upper2[::-1,1])*(1-beta) +(base_upper1[:,1])*beta
        
        lower_coords = np.zeros(base_lower1.shape)
        lower_coords[:,0] = base_lower1[:,0]
        lower_coords[:,1] = (base_lower1[:,1]*beta + np.interp(lower_coords[:,0], base_lower2[:,0], base_lower2[:,1])*(1-beta))
        
        return Airfoil(np.vstack((upper_coords, lower_coords)))
        
    def _get_upper(self):
        return self.coords[:self._n,:]
        
    def _get_lower(self):
        return self.coords[self._n:, :]
    
    def plot(self, *args):
        from matplotlib import pyplot as plt
        
        plt.plot(self.coords[:,0], self.coords[:,1], *args)


class Section(object):
    """
    A storage class for wing sections
    """
    def __init__(self, pos: Point, chord: float, alpha: float, airfoil: str):
        self.pos = pos
        self.alpha = alpha
        self.chord = chord
        self.airfoil = airfoil

    def __lt__(self, other) -> bool:
        return self.pos.x < other.pos.x

    def __eq__(self, other):
        return self.pos.x == other.pos.x


class BaseWing(object):
    """
    A basic storage ClassF
    """

    def __init__(self, pos):
        self.sections = SortedList()
        self.pos = pos
        self.root_pos = 0.0

    @classmethod
    def create_from_planform(cls, span_positions: list, chord_lengths: list, offsets: list, twists: list, airfoils: list):
        """Generates a Wing object faster than through adding sections manually"""

        pos1 = Point(0, 0, 0)
        wing1 = cls(pos1)

        data = zip(span_positions, chord_lengths, offsets, twists, airfoils)

        for span_position, chord_length, offset, twist, airfoil in data:
            pos = Point(offset, span_position, 0)
            section = Section(pos, chord_length, twist, airfoil)

            wing1.add_section(section)

        return wing1

    def add_section(self, section: Section) -> None:

        self.sections.add(section)

    def set_root_pos(self, span_position: float) -> None:
        self.root_pos = span_position

    def area(self) -> float:
        """Calculate area of the wing."""

        span_positions = [section.pos.y for section in self.sections]
        chord_lengths = [section.chord for section in self.sections]
        
        while self.root_pos > span_positions[0]:
            last = (span_positions[0], chord_lengths[0])
            
            del(span_positions[0])
            del(chord_lengths[0])
        
        if self.root_pos < span_positions[0]:
            chord0 = self.root_pos-last[0]/(span_positions[0]-last[0])*(chord_lengths[0]-last[1])+last[1]
        
        area = np.trapz(chord_lengths, span_positions)

        return 2*area # for both sides of wing

    def aspect_ratio(self) -> float:
        return self.span_width()**2/self.area()

    def mac(self) -> Section:
        """Calculate the mean aerodynamic chord of wing."""

        def generate_segment_macs(x_positions, y_positions, chord_lengths):
            """Generate the segments mean aerodynamic chord values."""

            last = None

            for x, y, chord in zip(x_positions, y_positions, chord_lengths):

                if last is not None and y > self.root_pos:
                    x_old, y_old, chord_old = last
                    span = y - y_old
                    area = (chord_old+chord) * span / 2
                    taper = chord / chord_old
                    taper1 = taper + 1
                    fraction = (taper + taper1) / (3 * taper1)

                    mac = chord_old * (taper**2 + taper + 1) / (1.5 * taper1)
                    mac_x = x_old + fraction * (x - x_old)
                    mac_y = y_old + fraction * span

                    yield (mac_x*area, mac_y*area, mac*area, area)

                last = (x, y, chord)

        x_positions = [section.pos.x for section in self.sections]
        y_positions = [section.pos.y for section in self.sections]
        chord_lengths = [section.chord for section in self.sections]

        res = np.array(list(generate_segment_macs(x_positions, y_positions, chord_lengths)))

        area = self.area()

        mac_x = sum(res[:, 0])/area
        mac_y = sum(res[:, 1])/area
        mac = sum(res[:, 2])/area

        return Section(Point(mac_x, mac_y, 0),  mac, 0, airfoil=None)

    def span_width(self) -> float:
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

        y_positions = [section.pos.y for section in self.sections]

        airfoil_index = np.argmin(np.array(y_positions)-y)

        return self.sections[airfoil_index].airfoil
            
    def airfoils_at(self, y:float) -> tuple:
        """ Lookup airfoils around given span position"""
        
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


class Flap(object):
    """
    A storage class for flaps.
    """
    def __init__(self, span_start, span_end, depth):
        """ Create a Flap object

        Arguments
        span_start -- span at which flap begins
        span_end -- span at which flap ends
        depth -- relative depth at span_start and _end as list
        """
        self.y_start = span_start
        self.y_end = span_end
        self.depth_start = depth[0]
        self.depth_end = depth[1]
        
    def depth_at(self, span_pos):
        return np.interp(span_pos, [self.y_start, self.y_end], [self.depth_start, self.depth_end],
                         right=0.0, left=0.0)

    def __lt__(self, other) -> bool:
        return self.y_start < other.y_start

    def __eq__(self, other):
        return self.y_start == other.y_start


class Wing(BaseWing):
    """
    A storage class for Wing with airbrake and flaps.
    """
    def __init__(self, pos):
        super(Wing, self).__init__(pos)
        self.flaps = dict()
        self.airbrake = None
        
    def set_flap(self, name, span_pos_start, span_pos_end, depth):
        
        flap = Flap(span_pos_start, span_pos_end, depth)
        
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
            x = section.pos.x
            y = section.pos.y
            chord = section.chord
            
            plt.plot((y, y), (-x, -x-chord), 'r')
            x_positions.append(x)
            y_positions.append(y)
            chord_lengths.append(chord)
        
        y_positions = np.array(y_positions)
        
        # draw leading edge
        plt.plot(y_positions, -np.array(x_positions), 'b' )
        # draw trailing edge
        plt.plot(y_positions, -np.array(x_positions)-np.array(chord_lengths), 'b')
        
        # draw flaps
        for name, aflap in self.flaps.items():
            y_pos = np.array([aflap.y_start, aflap.y_end])
            y_pos = np.concatenate([y_pos, y_positions[(y_positions>aflap.y_start)&(y_positions<aflap.y_end)]])
            y_pos.sort()
            
            chords = np.interp(y_pos, y_positions, chord_lengths)
            offsets = np.interp(y_pos, y_positions, x_positions)
            te = -(chords+offsets)
            depth = aflap.depth_at(y_pos)*chords
            
            plt.plot(y_pos[[0,*range(len(y_pos)),-1]], [te[0],*(te+depth),te[-1]], 'g--')
        # format 
        plt.axis('equal')
        plt.axis('off')
        plt.xlim(-1, max(y_positions)+1)
        
class Plane:
    def __init__(self, name, wing, hlw):
        self.name = name
        self.wing = wing
        self.hlw = hlw

    @classmethod
    def load_from_file(cls, filename):
        import yaml, io

        # Lade Daten
        with io.open(filename, encoding='utf8') as f:
            content = f.read()
            data = yaml.load(content)

        # Erstelle Flugzeug Objekt
        plane = cls.generate(data)

        return plane

    @classmethod
    def generate(cls, data):

        # Name einlesen
        name = data['Name']

        # Flügel einlesen
        fluegel = BaseWing.generate(data[u'Flügel'])

        # Höhenleitwerk einlesen
        Hoehenleitwerk = namedtuple('Hoehenleitwerk','pos c4tel ')

        hlw_data = data[u'Höhenleitwerk']

        hlw_pos = Point(x=hlw_data['pos']['x'], y=0, z=hlw_data['pos']['z'])
        hlw = Hoehenleitwerk(pos=hlw_pos, c4tel=hlw_data['c4tel'])

        # Flugzeug Objekt erstellen
        plane = Plane(name, wing=fluegel, hlw=hlw)

        return plane

    def plot(self):
        from matplotlib import pyplot as plt

        fig = plt.figure(figsize=(10,6))

        # 3D Plot
        ax = fig.add_subplot(221, projection='3d')

        xle = self.wing.xle
        yle = self.wing.yle
        zle = self.wing.zle
        c = self.wing.c
        kt = self.wing.kt

        ax.plot(xle, yle, zle)

        # top view
        plt.subplot(2,2,2)

        plt.plot(yle, xle, 'black')      # leading edge
        plt.plot(yle, xle+c, 'black')    # trailing edge
        plt.plot(yle, xle+(1-kt)*c, '*') # flap edge
        plt.axis('equal')

        # side view
        plt.subplot(2,2,3)
        
def createDummy():
    return None
