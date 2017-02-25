# -*- coding: utf-8 -*-

from collections import namedtuple
from sortedcontainers import SortedList
import numpy as np

Point = namedtuple('Point','x y z')


class Section(object):

    def __init__(self, pos: Point, chord: float, alpha: float, airfoil: str):
        self.pos = pos
        self.chord = chord
        self.airfoil = airfoil

    def __lt__(self, other) -> bool:
        return self.pos.x < other.pos.x

    def __eq__(self, other):
        return self.pos.x == other.pos.x


class Wing(object):
    """describes symmetric wing"""

    def __init__(self, pos):
        self.sections = SortedList()
        self.pos = pos

    @classmethod
    def fast_creation(cls, span_positions: list, chord_lengths: list, offsets:list, twists:list):
        """Generates a Wing object faster than through adding sections manually"""

        pos1 = Point(0, 0, 0)
        wing1 = cls(pos1)

        data = zip(span_positions, chord_lengths, offsets, twists)

        for span_position, chord_length, offset, twist in data:
            pos = Point(offset, span_position, 0)
            section = Section(pos, chord_length, twist, None)

            wing1.add_section(section)

        return wing1

    def add_section(self, section: Section):

        self.sections.add(section)

    def area(self) -> float:
        """Calculate area of the wing."""

        span_positions = [section.pos.y for section in self.sections]
        chord_lengths = [section.chord for section in self.sections]

        area = np.trapz(chord_lengths, span_positions)

        return area

    def mac(self) -> Section:
        """Calculate the mean aerodynamic chord of wing."""

        def generate_segment_macs(x_positions, y_positions, chord_lenghts):
            """Generate the segments mean aerodynamic chord values."""

            last = None

            for x, y, chord in zip(x_positions, y_positions, chord_lenghts):

                if last is not None:
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

        print(res)

        area = self.area()

        mac_x = sum(res[:, 0])/area
        mac_y = sum(res[:, 1])/area
        mac = sum(res[:, 2])/area

        return Section(Point(mac_x, mac_y, 0),  mac, 0, airfoil=None)

    def span_width(self) -> float:
        """Calculate the span width of wing."""
        return 2 * max((section.pos.y for section in self.sections))

    def chord_at(self, y:float) -> float:
        """Calculates the chord depth at given span position"""

        y_positions = [section.pos.y for section in self.sections]
        chord_lengths = [section.chord for section in self.sections]

        return np.interp(y, y_positions, chord_lengths)


    def airfoil_at(self, y: float) -> str:
        """Lookup airfoil at given span position"""

        y_positions = [section.pos.y for section in self.sections]

        airfoil_index = np.argmin(np.array(y_positions)-y)

        return self.sections[airfoil_index].airfoil

    def _svg_(self):
        import svgwrite

        dwg = svgwrite.Drawing()

        # draw contourline
        x_positions = [section.pos.x for section in self.sections]
        y_positions = [section.pos.y for section in self.sections]
        chord_lengths = [section.chord for section in self.sections]

        leading_edge = [(x*100, y*100) for x, y in zip(x_positions, y_positions)]
        trailing_edge = [((x+c)*100, y*100) for x, y, c in zip(x_positions, y_positions, chord_lengths)]

        all_pts = leading_edge+trailing_edge[::-1]

        contour = dwg.add(dwg.g(id='contour', stroke='black'))

        contour.add(dwg.polyline(all_pts))

        return dwg.tostring()

        # TODO: fix this function


    def _repr_html_(self):
        return """
        <h3>Wing</h3>
        {3}
        <p><b>spanwidth:</b>{0}</p>
        <p><b>wing area:</b>{1}</p>
        <p><b>mac</b>{2}</p>""".format(self.span_width(), self.area(), self.mac().chord, self._svg_())


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
        fluegel = Wing.generate(data[u'Flügel'])

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