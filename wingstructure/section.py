from shapely.geometry import LineString, Polygon
from shapely.ops import cascaded_union
import numpy as np


class SectionFeature:
    pass


class SparI(SectionFeature):
    def __init__(self, relpos: float, flangethickness: float, webthickness: float, width: float, density: float):
        self._flangethickness = flangethickness
        self._webthickness = webthickness
        self._width = width
        self._pos = relpos
        self._density = density

        self.flangetop = None
        self.flangebottom = None
        self.web = None

        self.geometry = None

    def update(self, shell):
        xmin, xmax = shell.geometry.bounds[0], shell.geometry.bounds[2]

        pos = xmin + (xmax - xmin) * self._pos

        print((xmin,xmax, pos, self._width))

        top, bottom = shell.interiorheight(pos - self._width/2, pos + self._width/2)

        self.flangetop = LineString([(pos - self._width/2, top), (pos + self._width/2, top)])
        self.flangebottom = LineString([(pos - self._width/2, bottom), (pos + self._width/2, bottom)])

        self.web = LineString([(pos, bottom), (pos, top)])

        self.geometry = cascaded_union([self.flangebottom, self.flangetop, self.web])

    def mass(self) -> float:
        flangeweight = (self.flangebottom.length + self.flangetop.length) * self._flangethickness * self._density
        webweight = self.web.length * self._webthickness * self._density

        return flangeweight+webweight

    def centerofgravity(self):
        flangetop_weighted_cg = np.array(self.flangetop.centroid) * self.flangetop.length * self._flangethickness * self._density
        flangebottom_weighted_cg = np.array(self.flangebottom.centroid) * self.flangebottom.length * self._flangethickness* self._density
        web_cg = np.array(self.web.centroid) * self.web.length * self._webthickness * self._density

        return (flangebottom_weighted_cg + flangetop_weighted_cg + web_cg) / self.mass()


class Shell(SectionFeature):
    def __init__(self, airfoil: np.array, thickness: float,  density: float):

        try:
            self.geometry = LineString(airfoil)
        except:
            raise ValueError('Parameter airfoil does not match.')

        self.thickness = thickness
        self.density = density

    def mass(self) -> float:
        return self.geometry.length*self.thickness*self.density

    def centerofgravity(self) -> np.array:
        return np.array(self.geometry.centroid)

    def interiorheight(self, xmin, xmax) -> np.array:

        bounds = self.geometry.bounds

        poly = Polygon([(xmin, bounds[1]-0.1), (xmax, bounds[1]-0.1), (xmax, bounds[3]+0.1), (xmin, bounds[3]+0.1)])

        intersection = self.geometry.intersection(poly)

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


class Section:
    def __init__(self, shell):

        self._shell = shell
        self._features = []

        self.geometry = None

    def addfeature(self, feature: SectionFeature) -> None:
        self._features.append(feature)

    def mass(self) -> float:
        return self._shell.mass() + sum([feature.mass() for feature in self._features])

    def centerofgravity(self) -> np.array:
        shell_weighted_cg = self._shell.mass() * self._shell.centerofgravity()
        features_weighted_cg = sum([feature.mass() * feature.centerofgravity() for feature in self._features])

        return (shell_weighted_cg + features_weighted_cg) / self.mass()

    def update(self):
        for feature in self._features:
            feature.update(self._shell)

        self.geometry = self._shell.geometry.union(cascaded_union([feature.geometry for feature in self._features]))
