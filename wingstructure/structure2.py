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

    @property
    def cut_elements(self):
        return []

class _AbstractBaseStructure(_AbstractBase):
    def __init__(self, parent, material):
        super().__init__()
        self.parent = parent
        self.material = material
        self._cut_elements = []
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

    @property
    def massproperties(self):
        return self.geometry.centroid, self.geometry.area*self.material.ρ

    def _sum_massproperties(self):
        mass = 0.0
        cg = np.zeros(2)
        for child in self.children:
            childcg, childmass = child.massproperties

            mass += childmass
            cg += childmass * np.array(childcg)
        
        return geom.Point(cg/mass), mass

    def _get_inside_direction(self, linearring):
        """Gets the inside direction for parallel offset (left or right)
        from signed area of geometry"""
        
        if cga.signed_area(linearring) > 0:
            return 'left'
        else:
            return 'right'
    
    @property
    def cut_elements(self):
        return [*self.parent.cut_elements,*self._cut_elements]


class SectionBase(_AbstractBase):
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

        string = '<svg viewBox=\"{}\" xmlns="http://www.w3.org/2000/svg">\
                    <defs>\
                        <!-- simple dot marker definition -->\
                        <marker id="dot" viewBox="0 0 10 10" refX="5" refY="5"\
                            markerWidth="3" markerHeight="3">\
                        <circle cx="5" cy="5" r="5" fill="red" />\
                        </marker>\
                    </defs>\
                  {}</svg>'.format("0, 0, {}, {}".format(width, height),"{}")

        heigth_offset = -bounds[3]
        width_offset = bounds[0] 

        return string, (width_offset, heigth_offset)

    def _repr_svg_(self):
        def collectgeometry(part):
            collection = []
            for child in part.children:
                collection.extend(collectgeometry(child))
            collection.append(part.geometry)
            return collection

        collection = collectgeometry(self)

        shply_collection = geom.GeometryCollection(collection)

        return shply_collection._repr_svg_()
    
    def massproperties(self):
        return self.children[0]._sum_massproperties()


class Layer(_AbstractBaseStructure):
    def __init__(self, parent, material, thickness=0.0):
        super().__init__(parent, material)
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

        for elem in self.cut_elements:
            self.geometry -= elem

    def _repr_svg_(self):
        string, offset = self.svgdata
        return string.format(svgpolyline(np.array(self.geometry.interiors[0].coords),
                                                    offset)+\
                            svgpolyline(np.array(self.geometry.exterior.coords),
                                                    offset))


class Reinforcement(Layer):
    def __init__(self, parent, material, thickness=0.0, limits=None):
        _AbstractBaseStructure.__init__(self, parent, material)#TODO fix
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

        tmp_interior = geom.Polygon(exterior)

        for ageo in intersection.geoms:
            tmp_geometry = self._create_offset_box(ageo, self.thickness, side, symmetric=False)
            tmp_interior -= tmp_geometry
            tmp_geometries.append(tmp_geometry)

        self.geometry = geom.GeometryCollection(tmp_geometries)

        self.interior = geom.LinearRing(tmp_interior.exterior)
        self._refine_interior()

    def _refine_interior(self):
        
        if self.interior.type == 'MultiLineString':
           self.interior = geom.LinearRing(self.interior.geoms[0])
        else:
            self.interior = geom.LinearRing(self.interior)

        # find bugs in interior and remove them
        pts = np.array(self.interior)[:-1,:]
        del_pts = []

        for ii in range(len(pts)):
            
            vec1 = pts[ii]-pts[ii-1]
            vec2 = pts[(ii+1)%len(pts)]-pts[ii]
            
            if(np.linalg.norm(vec1/np.linalg.norm(vec1)+vec2/np.linalg.norm(vec2))<1e-3):
                del_pts.append(ii)

        pts2 = np.delete(pts, del_pts, axis=0)
        
        self.interior = geom.LinearRing(pts2)

    def _repr_svg_(self):
         string, offset = self.svgdata
         return string.format(svgpolyline(np.array(self.geometry[0].exterior.coords),
                                                    offset)+\
                            svgpolyline(np.array(self.geometry[1].exterior.coords),
                                                    offset))

    def _create_offset_box(self, line, thickness, side, bevel=0.0, symmetric=False ):
        
        offsetline = line.parallel_offset(thickness, side=side)
        
        if symmetric:
            offsetline2 = line.parallel_offset(thickness, side=otherside(side))
            line = offsetline2

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

class Display(object):
    def __init__(self, geometry):
        self.geometry = geometry
    
    def _repr_svg_(self):
        bounds = self.geometry.bounds
        height = bounds[3]-bounds[1]
        width = bounds[2]-bounds[0]

        string = '<svg viewBox=\"{}\" xmlns="http://www.w3.org/2000/svg">\
                    <defs>\
                        <!-- simple dot marker definition -->\
                        <marker id="dot" viewBox="0 0 10 10" refX="5" refY="5"\
                            markerWidth="5" markerHeight="5">\
                        <circle cx="5" cy="5" r="5" fill="red" />\
                        </marker>\
                    </defs>\
                  {}</svg>'.format("0, 0, {}, {}".format(width, height),"{}")

        offset = [bounds[0], -bounds[3]]

        try:
            pts = np.array(self.geometry)
        except:
            try:
                pts = np.array(self.geometry.exterior)
            except:
                return ''

        return string.format(svgpolyline(pts, offset))


class BoxSpar(_AbstractBaseStructure):
    def __init__(self, parent, material, midpos: float, width: float, flangeheigt, webwidth):
        super().__init__(parent, material)
        self.midpos = midpos
        self.width = width
        self.interior = None
        self.flangeheigt = flangeheigt
        self.webwidth = webwidth

        self._update_geometry(parent.interior)

    def _update_geometry(self, exterior):
        self.interior = exterior
        start = self.midpos-self.width/2
        end = self.midpos+self.width/2

        abox = geom.box(start, exterior.bounds[1]-3, end, exterior.bounds[3]+3)

        intersection = abox.intersection(geom.Polygon(exterior))

        for elem in self.cut_elements:
            if not intersection.intersects(elem):
                continue
            raise Exception('Cannot use cut element in cut element')

        anotherbox = geom.box(start+self.webwidth,
                              exterior.bounds[1]-3,
                              end-self.webwidth,
                              exterior.bounds[3]+3)

        web = intersection-anotherbox
        flangebox = intersection.intersection(anotherbox)

        flangcutout = geom.box(start, flangebox.bounds[1]+self.flangeheigt,
                               end, flangebox.bounds[3]-self.flangeheigt)

        flange = flangebox-flangcutout

        self.geometry = geom.GeometryCollection([*web.geoms, *flange.geoms])       

        self._cut_elements = [abox]

    @property
    def massproperties(self):
        if type(self.material) is dict:
            try:
                webmaterial = self.material['web']
                flangematerial = self.material['flange']
            except:
                raise Exception('wrong material defintion')

            mass = 0.0
            cg = np.zeros(2)

            for ii, geom in enumerate(self.geometry.geoms):
                if ii<2:
                    material = webmaterial
                else:
                    material = flangematerial
                
                tmpmass = geom.area*material.ρ 
                mass += tmpmass
                cg += np.array(geom.centroid) * tmpmass
            mass = (self.geometry.geoms[0].area+self.geometry.geoms[1])*webmaterial.ρ \
                    + (self.geometry.geoms[2].area+self.geometry.geoms[3])*flangematerial.ρ 
            cg /= mass

            return geom.Point(cg), mass
        else:
            return super().massproperties


def svgpolyline(pts, offset):
    string = "<polyline points=\"{}\" fill=\"none\" stroke=\"black\" marker-start=\"url(#dot)\" marker-mid=\"url(#dot)\"  marker-end=\"url(#dot)\"  />"

    return string.format(svgpts(pts, offset))

def svgpts(pts, offset=(0,0)):
    pts[:,1] *= -1
    pts[:,1] -= offset[1]
    pts[:,0] -= offset[0]
    
    return " ".join(("{},{}".format(*row) for row in pts))

def otherside(side):
    if side=='left':
        return 'right'
    else:
        return 'left'
