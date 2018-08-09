from abc import ABC, abstractmethod

import numpy as np

from shapely import geometry as geom, ops
from shapely.algorithms import cga


class _AbstractBase(ABC):
    def __init__(self):
            self.children = []
            self.geometry = None
    
    def _update(self, interior):
        for child in self.children:
            child._update()
    
    def _addchild(self, child):
        self.children.append(child)

class _AbstractBaseStructure(_AbstractBase):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        parent._addchild(self)

    @property
    def svgdata(self):
        return self.parent.svgdata

    def _update(self, interior):
        self._update_geometry(interior)

        for child in self.children:
            child._update()

    @abstractmethod
    def _update_geometry(self, exterior):
        pass

    def _get_inside_direction(self, linearring):
        """Gets the inside direction for parallel offset (left or right)
        from signed area of geometry"""
        
        if cga.signed_area(linearring) > 0:
            return 'left'
        else:
            return 'right'


class BaseAirfoil(_AbstractBase):
    def __init__(self, airfoil_coordinates):
        super().__init__()
        self.geometry = geom.LinearRing(airfoil_coordinates)

    @property
    def interior(self):
        return self.geometry

    @property
    def svgdata(self):
        bounds = self.geometry.bounds
        height = bounds[3]-bounds[1]
        width = bounds[2]-bounds[0]

        string = """<svg viewBox=\"{}\" xmlns="http://www.w3.org/2000/svg">
                  {}</svg>""".format("0, 0, {}, {}".format(width, height),"{}")

        heigth_offset = -bounds[3]
        width_offset = bounds[0] 

        return string, (width_offset, heigth_offset)

    def _repr_svg_(self):
        string, offset = self.svgdata
        pts = np.array(self.geometry.coords)
        return string.format(svgpolyline(pts, offset))


class Layer(_AbstractBaseStructure):
    def __init__(self, parent, thickness=0.0):
        super().__init__(parent)
        self.thickness = thickness

        self.interior = None

        self._update_geometry(parent.interior)

    def _update_geometry(self, exterior):
        
        inside_direction = self._get_inside_direction(exterior)

        self.interior = exterior.parallel_offset(self.thickness, side=inside_direction)

        if self.interior.type == 'MultiLineString':
                    self.interior = geom.LinearRing(self.interior.geoms[0])
        else:
            self.interior = geom.LinearRing(self.interior)

        self.geometry = geom.Polygon(exterior)-geom.Polygon(self.interior)

    def _repr_svg_(self):
        string, offset = self.svgdata
        return string.format(svgpolyline(np.array(self.geometry.interiors[0].coords),
                                                    offset)+\
                            svgpolyline(np.array(self.geometry.exterior.coords),
                                                    offset))

class Reinforcement(Layer):
    def __init__(self, parent, thickness=0.0, limits=None):
        _AbstractBaseStructure.__init__(self, parent)#TODO fix
        self.thickness = thickness
        self.limits = limits

        self.interior = None

        self._update_geometry(parent.interior)
    
    def _update_geometry(self, exterior):
        limited_box = geom.box(self.limits[0], exterior.bounds[1]*1.1, 
                        self.limits[1], exterior.bounds[3]*1.1)

        intersection = limited_box.intersection(exterior)

        side = self._get_inside_direction(exterior)

        tmp_geometries = []

        self.interior = geom.LinearRing(exterior)

        for ageo in intersection.geoms:
            tmp_geometry = self._create_offset_box(ageo, self.thickness, side)
            self.interior -= tmp_geometry
            tmp_geometries.append(tmp_geometry)

        self.geometry = tmp_geometries

        self._refine_interior()

    def _refine_interior(self):
        if self.interior.type == 'MultiLineString':
           self.interior = geom.LinearRing(self.interior.geoms[1])
        else:
            self.interior = geom.LinearRing(self.interior)

    def _repr_svg_(self):
         string, offset = self.svgdata
         return string.format(svgpolyline(np.array(self.geometry[0].exterior.coords),
                                                    offset)+\
                            svgpolyline(np.array(self.geometry[1].exterior.coords),
                                                    offset))
    def _create_offset_box(self, line, thickness, side, bevel=0.0):
        
        offsetline = line.parallel_offset(thickness, side=side)
        
        if bevel > 0.0:
            
            raise Exception('Beveling not yet implemented')
        
        if side == 'left':
            connect1 = geom.LineString((line.coords[-1],offsetline.coords[-1]))
            connect2 = geom.LineString((line.coords[0], offsetline.coords[0]))
            return geom.Polygon(ops.linemerge((line,connect1,offsetline,connect2)))
        else:
            connect1 = geom.LineString((line.coords[-1],offsetline.coords[0]))
            connect2 = geom.LineString((line.coords[0], offsetline.coords[-1]))
            return geom.Polygon(ops.linemerge((offsetline,connect1,line,connect2)))

def svgpolyline(pts, offset):
    string = "<polyline points=\"{}\" fill=\"none\" stroke=\"black\" />"

    return string.format(svgpts(pts, offset))

def svgpts(pts, offset=(0,0)):
    pts[:,1] *= -1
    pts[:,1] -= offset[1]
    pts[:,0] -= offset[0]
    
    return " ".join(("{},{}".format(*row) for row in pts))
