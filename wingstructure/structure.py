from shapely.geometry import LineString, Polygon
from shapely.ops import cascaded_union
import numpy as np


class SectionFeature(object):
    def __init__(self):
        self._shell = None

    def update(self):
        if self._shell:

            self._update(self)

        else:

            raise Exception('To update feature it has to be part of structure!')


class SparI(SectionFeature):
    def __init__(self, relpos: float, flangethickness: float, webthickness: float, width: float, density: float):
        super().__init__()
        self._flangethickness = flangethickness
        self._webthickness = webthickness
        self._width = width
        self._pos = relpos
        self._density = density

        self.flangetop = None
        self.flangebottom = None
        self.web = None

        self.geometry = None

    def _update(self):
        shell = self._shell

        xmin, xmax = shell.geometry.bounds[0], shell.geometry.bounds[2]

        pos = xmin + (xmax - xmin) * self._pos

        print((xmin,xmax, pos, self._width))

        top, bottom = shell.interiorheight(pos - self._width/2, pos + self._width/2)

        self.flangetop = LineString([(pos - self._width/2, top), (pos + self._width/2, top)])
        self.flangebottom = LineString([(pos - self._width/2, bottom), (pos + self._width/2, bottom)])

        self.web = LineString([(pos, bottom), (pos, top)])

        self.geometry = cascaded_union([self.flangebottom, self.flangetop, self.web])

    def _sec_mass(self) -> float:
        flangeweight = (self.flangebottom.length + self.flangetop.length) * self._flangethickness * self._density
        webweight = self.web.length * self._webthickness * self._density

        return flangeweight+webweight

    def _sec_cg(self):
        flangetop_weighted_cg = np.array(self.flangetop.centroid) * self.flangetop.length * self._flangethickness * self._density
        flangebottom_weighted_cg = np.array(self.flangebottom.centroid) * self.flangebottom.length * self._flangethickness* self._density
        web_cg = np.array(self.web.centroid) * self.web.length * self._webthickness * self._density

        return (flangebottom_weighted_cg + flangetop_weighted_cg + web_cg) / self._sec_mass()


class Shell(SectionFeature):
    def __init__(self, winggeometry, airfoildb, thickness: float, density: float):

        self.winggeometry = winggeometry
        self.airfoildb = airfoildb
        self.thickness = thickness
        self.density = density

    def mass(self, y: float) -> float:
        return self.geometry(y).length*self.thickness*self.density

    def _sec_cg(self, y) -> np.array:
        return np.array(self.geometry(y).centroid)

    def interiorheight(self, y: float, x_min: float, x_max: float) -> np.array:

        bounds = self.geometry(y).bounds

        poly = Polygon([(x_min, bounds[1] - 0.1), (x_max, bounds[1] - 0.1), (x_max, bounds[3] + 0.1), (x_min, bounds[3] + 0.1)])

        intersection = self.geometry(y).intersection(poly)

        y = np.zeros((2, 2))

        for ii, geo in enumerate(intersection.geoms):
            y[ii, :] = [geo.bounds[1], geo.bounds[3]]

        if y[0, 0] > y[1, 0]:
            top = min(y[0, :])
            bottom = max(y[1, :])
        else:
            top = min(y[1, :])
            bottom = max(y[0, :])

        return top, bottom

    def geometry(self, y):
        try:
            airfoil_str = self.winggeometry.airfoil_at(y)
            airfoil = self.airfoildb[airfoil_str]*self.winggeometry.chord_at(y)
            return LineString(airfoil)
        except:
            raise ValueError('Error while creating geometry.')


class Structure(object):
    def __init__(self, shell):

        self._shell = shell
        self._features = []

        self.geometry = None

    def addfeature(self, feature: SectionFeature) -> None:
        self._features.append(feature)
        self._features[-1].shell = self._shell

    def _sec_mass(self, y) -> float:
        return self._shell.mass() + sum([feature.mass() for feature in self._features])

    def _sec_cg(self, y) -> np.array:
        shell_weighted_cg = self._shell.mass() * self._shell.centerofgravity()
        features_weighted_cg = sum([feature.mass() * feature.centerofgravity() for feature in self._features])

        return (shell_weighted_cg + features_weighted_cg) / self.mass()

    def update(self):
        for feature in self._features:
            feature.update(self._shell)

        self.geometry = self._shell.geometry.union(cascaded_union([feature.geometry for feature in self._features]))
