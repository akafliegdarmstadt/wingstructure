"""Module for representation of wing section's structure

This module defines classes which allow definition of a wing section as chain of multiple structure elements.
Definition of a structure starts with SectionBase, which stores airfoil coordinates. Beginning with this exterior
geometry structure elements can be added from out to inside. Possible Elements are Layers, covering the whole
surface, Reinforcements for higher stiffness locally and Spar elements (ISpar and BoxSpar).

Example
-------
.. code-block:: python

   import numpy as np
   from wingstructure import structure

   # Load cooardinates
   coords = np.loadtxt('ah93157.dat', skiprows=1) * 1.2
   sectionbase = structure.SectionBase(coords)

   # Define Material
   import collections
   Material = collections.namedtuple('Material', ['ρ'])

   carbonfabric = Material(1.225)
   foam = Material(1.225)
   sandwich = Material(1.225)

   # Add layers
   outerlayer = structure.Layer(sectionbase, carbonfabric, 5e-4)
   core = structure.Layer(outerlayer, foam, 1e-2)
   innerlayer = structure.Layer(core, carbonfabric, 5e-4)

   # Add Spar
   spar = structure.ISpar(parent=innerlayer, 
                          material={'flange': carbonfabric, 'web': sandwich},
                          midpos=0.45,
                          flangewidth=0.2,
                          flangethickness=0.03,
                          webpos=0.5,
                          webthickness=0.02)

   # Analyse Mass
   massana = structure.MassAnalysis(spar)
   cg, mass = massana.massproperties()
"""

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

    def area(self):
        return self.geometry.area

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
            offsetline2 = line.parallel_offset(thickness, side=_oderside(side))
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
    """Foundation for section's wing structure description
    
    Parameters
    ----------
    airfoil_coordinates : np.ndarray
        airfoils coordinates from dat file
    
    Attributes
    ----------
    interior: shapely.geometry.LinearRing
        representation of airfoil as shapely geometry
    
    """

    def __init__(self, airfoil_coordinates):
        super().__init__()
        self._geometry = shpl_geom.LinearRing(airfoil_coordinates)

    @property
    def interior(self):
        return self._geometry
    
    def _repr_svg_(self):
        return self._geometry._repr_svg_()


class Layer(_AbstractBaseStructure):
    """Layer of constant thickness representation
    
    Parameters
    ----------
    parent : 
        Structure Element
    material :
        Material definition
    thickness : float
        The Layers geometric thickness

    Attributes
    ----------
    interior :
        shapely geometry representation of interior
    
    """

    def __init__(self, parent, material, thickness=0.0):
        super().__init__(parent, material)
        self._thickness = thickness

        self.interior = None

        self._update_geometry(parent.interior)

    def _update_geometry(self, exterior):
        
        inside_direction = self._get_inside_direction(exterior)

        self.interior = exterior.parallel_offset(self._thickness,
                                                 side=inside_direction)

        if self.interior.type == 'MultiLineString':
                    self.interior = shpl_geom.LinearRing(self.interior.geoms[0])
        else:
            self.interior = shpl_geom.LinearRing(self.interior)

        self.geometry = shpl_geom.Polygon(exterior)-shpl_geom.Polygon(self.interior)

    def _centerline(self):
        exterior = self.parent.interior

        inside_direction = self._get_inside_direction(exterior)

        centerline = exterior.parallel_offset(self._thickness/2.0,
                                              side=inside_direction)
        
        return centerline

    def middle_circumference(self):
        """Calculate circumference of centerline
        
        Returns
        -------
        float
            circumference value
        """

        return self._centerline().length

    def enclosed_area(self):
        """Calculate area enclosed by centerline

        Returns
        -------
        float
            ecnlosed area
        """
        
        cline = self._centerline()

        return shpl_geom.Polygon(cline).area

    @property
    def thickness(self):
        return self._thickness


class Reinforcement(Layer):
    """Local reinforcement structure
    
    Parameters
    ----------
    parent :
        Parent structure Element
    material :
        Material Defintion
    thickness :
        reinforcement thickness
    limits :
        bounds for reinforcement in chordwise direction

    Attributes
    ----------
    interior :
        shapely geometry representation of interior

    """

    def __init__(self, parent, material, thickness=0.0, limits=None):
        _AbstractBaseStructure.__init__(self, parent, material)  #TODO fix
        self._thickness = thickness
        self._limits = limits

        self.interior = None

        self._update_geometry(parent.interior)
    
    def _update_geometry(self, exterior):
        limited_box = shpl_geom.box(self._limits[0], exterior.bounds[1]*1.1,
                                    self._limits[1], exterior.bounds[3]*1.1)

        intersection = limited_box.intersection(exterior)

        side = self._get_inside_direction(exterior)

        tmp_geometries = []

        tmp_interior = shpl_geom.Polygon(exterior)

        for ageo in intersection.geoms:
            tmp_geometry = self._create_offset_box(ageo, self._thickness, side, 
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
    """Representation for I or double-T Spar
    
    Parameters
    ----------
    parent :
        Parent structure element.
    material :
        Material defintion.
    midpos : float
        Chord position of spar (absolute)
    flangewidth : float
        width of flanges
    flangethickness : float
        thickness of flanges (symmetric)
    webpos : float
        position of web at flange (relative, 0 - left, 1 - right)
    webthickness : float
        thickness of web


    Attributes
    ----------
    interior : 
        shapely geometry representation of interior
    
    Raises
    ------
    Exception
        when geometry cannot be created with choosen parameters

    """

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

    def webpos_abs(self):
        return self.midpos + (self.webpos - 0.5) * self.flangewidth

    def _update_geometry(self, exterior):
        self.interior = exterior
        start = self.midpos - self.flangewidth/2
        end = self.midpos + self.flangewidth/2

        if start<exterior.bounds[0] or end > exterior.bounds[2]:
            raise Exception('flange width is too large')

        cutbox = shpl_geom.box(start, exterior.bounds[1]-3,
                               end, exterior.bounds[3]+3)

        cutgeom = cutbox.intersection(exterior)

        self.tmp = cutgeom

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

        webpos_abs = start + self.webpos * self.flangewidth

        line1 = offsetlines[0]
        line2 = offsetlines[1]
        xdist1 = abs(line1.coords[0][0]-line2.coords[-1][0])
        xdist2 = abs(line1.coords[0][0]-line2.coords[0][0])

        if xdist1 < xdist2:
            webbox = shpl_geom.LinearRing([*line2.coords, *line1.coords])
        else:
            webbox = shpl_geom.LinearRing([*line1.coords[::-1], *line2.coords])

        webstart = webpos_abs - self.webthickness/2
        webend = webpos_abs + self.webthickness/2

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
    """Representation for box spar
    
    Parameters
    ----------
    parent :
        Parent structure element.
    material :
        Material defintion.
    midpos : float
        Chord position of spar (absolute)
    width : float
        width of spar
    flangeheight : float
        thickness of flanges (symmetric)
    webwidth : float
        thickness of web


    Attributes
    ----------
    interior : 
        shapely geometry representation of interior

    """

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
    """Analyse Mass of Structure
    
    Parameters
    ----------
    parent
        last element of structure elment chain to be analysed
    
    """

    def __init__(self, parent):
        self.parent = parent

    @property
    def massproperties(self):
        mass = 0.0
        cg = np.zeros(2)
            
        current = self.parent
        while not isinstance(current, SectionBase):

            cur_cg, cur_mass = current.massproperties

            mass += cur_mass
            cg += cur_mass * np.array(cur_cg)

            current = current.parent
        
        return cg/mass, mass


class LineIdealisation:
    """Idealisation of section for structural analysis

    Currently only section consisting of a I-Spar and one Layer can
    be idealized.

    Parameters
    ----------
    parent
        last element of structure feature chain to be analysed
    """

    def __init__(self, parent):
        
        self.parent = parent

        self.geometries = [None, None, None]

        self._chekc()

        self._update_geometry()

    def _chekc(self):

        self.spar = self.parent
        self.shell = self.spar.parent

        if not isinstance(self.spar, ISpar) or not isinstance(self.shell, Layer):
            raise Exception('Only limited section definition allowed for Idealization.')


    def _update_geometry(self):
        from shapely import geometry as shpl_geom
        bb = self.shell.parent.interior.bounds
        cutbox = shpl_geom.box(bb[0], bb[1], self.spar.webpos_abs(), bb[3])

        shell_center = self.shell._centerline()

        def opt_linemerge(geom):
            from shapely import ops
            if geom.type == 'MultiLineString':
                return ops.linemerge(geom)
            return geom

        geom_left = opt_linemerge(shell_center.intersection(cutbox))
        geom_right = opt_linemerge(shell_center.difference(cutbox))

        mid_start = np.array(geom_left.coords[0])
        mid_end = np.array(geom_left.coords[-1])

        coords_mid = mid_start + np.outer((mid_end-mid_start), np.linspace(0,1,40)).T

        self.geometries = np.array(geom_right.coords), coords_mid, np.array(geom_left.coords)
    
    @property
    def geometry(self):

        return self.geometries

    def _generate_datatuple(self, valtuple):

        datagenerator = (valtuple[i]*np.ones_like(self.geometries[i]) for i in range(3))

        return tuple(datagenerator)

    @property
    def thickness(self):
        
        shell_t = self.shell.thickness
        web_t = self.spar.webthickness

        return self._generate_datatuple((shell_t, web_t, shell_t))
    
    @property
    def youngsmoduli(self):
        shell_E = self.shell.material.E
        web_E = self.spar.material['web'].E

        return self._generate_datatuple((shell_E, web_E, shell_E))

    @property
    def shearmoduli(self):
        shell_G = self.shell.material.G
        web_G = self.shell.material.G

        return self._generate_datatuple((shell_G, web_G, shell_G))


def _oderside(side):
    if side == 'left':
        return 'right'
    else:
        return 'left'
