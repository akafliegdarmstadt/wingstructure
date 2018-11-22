from abc import ABC, abstractmethod

import numpy as np
from numpy.linalg import norm
from shapely import geometry as shpl_geom
from shapely import ops
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
        return (self.geometry.centroid, self.geometry.area*self.material.ρ)

    def _sum_massproperties(self):
        mass = 0.0
        cg = np.zeros(2)
        for child in self.children:
            childcg, childmass = child.massproperties

            mass += childmass
            cg += childmass * np.array(childcg)
        
        return shpl_geom.Point(cg/mass), mass

    def _get_inside_direction(self, linearring):
        """Gets the inside direction for parallel offset (left or right)
        from signed area of geometry"""
        
        if cga.signed_area(linearring) > 0:
            return 'left'
        else:
            return 'right'
    
    def _create_offset_box(self, line, thickness, side, bevel=0.0,
                           symmetric=False):
        
        offsetline = line.parallel_offset(thickness, side=side)
        
        if symmetric:
            offsetline2 = line.parallel_offset(thickness, side=otherside(side))
            line = offsetline2

        if bevel > 0.0:
            
            raise Exception('Beveling not yet implemented')
        
        if side == 'left':
            connect1 = shpl_geom.LineString((line.coords[-1], offsetline.coords[-1]))
            connect2 = shpl_geom.LineString((line.coords[0], offsetline.coords[0]))
            return shpl_geom.Polygon(ops.linemerge((line, connect1, offsetline, connect2)))
        else:
            connect1 = shpl_geom.LineString((line.coords[-1],offsetline.coords[0]))
            connect2 = shpl_geom.LineString((line.coords[0], offsetline.coords[-1]))
            return shpl_geom.Polygon(ops.linemerge((offsetline,connect1,line,connect2)))

    @property
    def cut_elements(self):
        return [*self.parent.cut_elements, *self._cut_elements]
    
    def _repr_svg_(self):
        def collectgeometry(part):
            
            current = self
            while not isinstance(current, SectionBase):
                oldcur = current
                current = current.parent

                yield oldcur.geometry
 
        collection = list(collectgeometry(self))

        shply_collection = shpl_geom.GeometryCollection(collection)

        return shply_collection._repr_svg_()


class SectionBase(_AbstractBase):
    def __init__(self, airfoil_coordinates):
        super().__init__()
        self.geometry = shpl_geom.LinearRing(airfoil_coordinates)

    @property
    def interior(self):
        return self.geometry
    
    def massproperties(self):
        return self.children[0]._sum_massproperties()

    def _repr_svg_(self):
        return self.geometry._repr_svg_()


class Layer(_AbstractBaseStructure):
    def __init__(self, parent, material, thickness=0.0):
        super().__init__(parent, material)
        self.thickness = thickness

        self.interior = None

        self._update_geometry(parent.interior)

    def _update_geometry(self, exterior):
        
        inside_direction = self._get_inside_direction(exterior)

        self.interior = exterior.parallel_offset(self.thickness,
                                                 side=inside_direction)

        if self.interior.type == 'MultiLineString':
                    self.interior = shpl_geom.LinearRing(self.interior.geoms[0])
        else:
            self.interior = shpl_geom.LinearRing(self.interior)

        self.geometry = shpl_geom.Polygon(exterior)-shpl_geom.Polygon(self.interior)


class Reinforcement(Layer):
    def __init__(self, parent, material, thickness=0.0, limits=None):
        _AbstractBaseStructure.__init__(self, parent, material)  #TODO fix
        self.thickness = thickness
        self.limits = limits

        self.interior = None

        self._update_geometry(parent.interior)
    
    def _update_geometry(self, exterior):
        limited_box = shpl_geom.box(self.limits[0], exterior.bounds[1]*1.1,
                                    self.limits[1], exterior.bounds[3]*1.1)

        intersection = limited_box.intersection(exterior)

        side = self._get_inside_direction(exterior)

        tmp_geometries = []

        tmp_interior = shpl_geom.Polygon(exterior)

        for ageo in intersection.geoms:
            tmp_geometry = self._create_offset_box(ageo, self.thickness, side, 
                                                   symmetric=False)
            tmp_interior -= tmp_geometry
            tmp_geometries.append(tmp_geometry)

        self.geometry = shpl_geom.GeometryCollection(tmp_geometries)

        self.interior = shpl_geom.LinearRing(tmp_interior.exterior)
        self._refine_interior()

    def _refine_interior(self):
        
        if self.interior.type == 'MultiLineString':
            self.interior = shpl_geom.LinearRing(self.interior.geoms[0])
        else:
            self.interior = shpl_geom.LinearRing(self.interior)

        # find bugs in interior and remove them
        pts = np.array(self.interior)[:-1, :]
        del_pts = []

        for ii in range(len(pts)):
            
            vec1 = pts[ii]-pts[ii-1]
            vec2 = pts[(ii+1) % len(pts)]-pts[ii]
            
            if(np.linalg.norm(vec1/norm(vec1)+vec2/norm(vec2)) < 1e-3):
                del_pts.append(ii)

        pts2 = np.delete(pts, del_pts, axis=0)
        
        self.interior = shpl_geom.LinearRing(pts2)


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


class ISpar(_AbstractBaseStructure):
    def __init__(self, parent, material, midpos: float, flangewidth: float,
                 flangethickness: float, webpos: float, webthickness: float):
        super().__init__(parent, material)
        self.midpos = midpos
        self.flangewidth = flangewidth
        self.flangethickness = flangethickness
        self.webpos = webpos
        self.webthickness = webthickness
        self.interior = None

        self._update_geometry(parent.interior)

    def _update_geometry(self, exterior):
        self.interior = exterior
        start = self.midpos - self.flangewidth/2
        end = self.midpos + self.flangewidth/2
        cutbox = shpl_geom.box(start, exterior.bounds[1]-3,
                               end, exterior.bounds[3]+3)

        cutgeom = cutbox.intersection(exterior)

        offsetside = self._get_inside_direction(exterior)

        flangegeoms = []
        offsetlines = []

        for cutline in cutgeom.geoms:
            flangegeoms.append(
                self._create_offset_box(cutline, self.flangethickness,
                                        offsetside)
            )
            offsetline = cutline.parallel_offset(self.flangethickness,
                                                 side=offsetside)
            offsetlines.append(offsetline)

        webpos = start + self.webpos * self.flangewidth

        line1 = offsetlines[0]
        line2 = offsetlines[1]
        xdist1 = abs(line1.coords[0][0]-line2.coords[-1][0])
        xdist2 = abs(line1.coords[0][0]-line2.coords[0][0])

        if xdist1 < xdist2:
            webbox = shpl_geom.LinearRing([*line2.coords, *line1.coords])
        else:
            webbox = shpl_geom.LinearRing([*line1.coords[::-1], *line2.coords])

        webstart = webpos - self.webthickness/2
        webend = webpos + self.webthickness/2

        webcutbox = shpl_geom.box(webstart, exterior.bounds[1]-3,
                             webend,  exterior.bounds[1]+3)

        web = webcutbox.intersection(shpl_geom.Polygon(webbox))

        self.geometry = shpl_geom.GeometryCollection([*flangegeoms, web])

        self._cut_elements = [cutgeom]

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
                if ii < 2:
                    material = flangematerial
                else:
                    material = webmaterial
                
                tmpmass = geom.area*material.ρ 
                mass += tmpmass
                cg += np.array(geom.centroid) * tmpmass

            cg /= mass

            return shpl_geom.Point(cg), mass
        else:
            return super().massproperties



class BoxSpar(_AbstractBaseStructure):
    def __init__(self, parent, material, midpos: float, width: float,
                 flangeheigt, webwidth):
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

        abox = shpl_geom.box(start, exterior.bounds[1]-3, end, exterior.bounds[3]+3)

        intersection = abox.intersection(shpl_geom.Polygon(exterior))

        for elem in self.cut_elements:
            if not intersection.intersects(elem):
                continue
            raise Exception('Cannot use cut element in cut element')

        anotherbox = shpl_geom.box(start+self.webwidth,
                              exterior.bounds[1]-3,
                              end-self.webwidth,
                              exterior.bounds[3]+3)

        web = intersection-anotherbox
        flangebox = intersection.intersection(anotherbox)

        flangcutout = shpl_geom.box(start, flangebox.bounds[1]+self.flangeheigt,
                               end, flangebox.bounds[3]-self.flangeheigt)

        flange = flangebox-flangcutout

        self.geometry = shpl_geom.GeometryCollection([*web.geoms, *flange.geoms])

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
                if ii < 2:
                    material = webmaterial
                else:
                    material = flangematerial
                
                tmpmass = geom.area*material.ρ 
                mass += tmpmass
                cg += np.array(geom.centroid) * tmpmass
            mass = (self.geometry.geoms[0].area+self.geometry.geoms[1].area)*webmaterial.ρ \
                    + (self.geometry.geoms[2].area+self.geometry.geoms[3].area)*flangematerial.ρ 
            cg /= mass

            return shpl_geom.Point(cg), mass
        else:
            return super().massproperties


class MassAnalysis:
    def __init__(self, parent):
        self.parent = parent

    def massproperties(self):
        mass = 0.0
        cg = np.zeros(2)
            
        current = self.parent
        while not isinstance(current, SectionBase):

            cur_cg, cur_mass = current.massproperties

            mass += cur_mass
            cg += cur_mass * np.array(cur_cg)

            current = current.parent
        
        return shpl_geom.Point(cg/mass), mass
        

def otherside(side):
    if side == 'left':
        return 'right'
    else:
        return 'left'
