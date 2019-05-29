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

   # create layers
   outerlayer = structure.Layer(carbonfabric, 5e-4)
   core = structure.Layer(foam, 1e-2)
   innerlayer = structure.Layer(carbonfabric, 5e-4)

   # create Spar
   spar = structure.ISpar(material={'flange': carbonfabric, 'web': sandwich},
                          midpos=0.45,
                          flangewidth=0.2,
                          flangethickness=0.03,
                          webpos=0.5,
                          webthickness=0.02)

    # add to sectionbase
    sectionbase.extend([outerlayer, core, innerlayer, spar])

   # Analyse Mass
   massana = structure.MassAnalysis(sectionbase)
   cg, mass = massana.massproperties()
"""

from abc import ABC, abstractmethod
from functools import wraps

import numpy as np
from collections import namedtuple
from numpy.linalg import norm
from shapely import geometry as shpl_geom
from shapely import ops
from shapely.algorithms import cga

from .section_helper import _get_inside_direction, _create_offset_box, _updatedecorator, _refine_interior


class SectionBase:
    """Foundation for section's wing structure description

    Implements a list like interface for the assembling of wing section's interior.
    
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
        self._geometry = shpl_geom.LinearRing(airfoil_coordinates)
        self._features = []
        self._featuregeometries = []

    @_updatedecorator
    def append(self, featuredef):
        self._features.append(featuredef)

    @_updatedecorator
    def insert(self, index, newfeature):
        self._features.insert(index, newfeature)

    @_updatedecorator
    def extend(self, featuredefs):
        for featuredef in featuredefs:
            self.append(featuredef)

    @_updatedecorator
    def remove(self, featuredef):
        if featuredef not in self._features:
            raise Exception('No feature {} found.'.format(featuredef))

        self._features.remove(featuredef)

    @_updatedecorator
    def pop(self):
        self._features.pop()

    def _update_callback(self, updated_feature):
        try:
            first_idx = self._features.index(updated_feature) #TODO: Error handling
        except ValueError:
            raise ValueError('Feature {} not in {}s features.'.format(updated_feature, self))

        self._update()
    
    def _update(self):

        # reset geometries
        self._featuregeometries = []

        self.tmp = []

        exterior = self._geometry

        for feature in self._features:
            tmp, tmp_geometry = feature._calculate_geometry(exterior)
            
            self.tmp.append(tmp)

            exterior=tmp
            self._featuregeometries.append(tmp_geometry)
    
    def _repr_svg_(self):

        shply_collection = shpl_geom.GeometryCollection([self._geometry, *self._featuregeometries])

        svg = shply_collection._repr_svg_()

        return rework_svg(svg, 1000, 250)

    def exportgeometry(self, refpoint=np.zeros(2)):

        geoms = []

        for feature, geometry in zip(self._features, self._featuregeometries):
            geoms.extend(exportstructure(geometry, feature.material, refpoint))
        
        return geoms

    def __getitem__(self, idx):
        return self._features[idx] # TODO rework


class Layer:
    """Layer of constant thickness representation
    
    Parameters
    ----------
    material :
        Material definition
    thickness : float
        The Layers geometric thickness

    Attributes
    ----------
    interior :
        shapely geometry representation of interior
    
    """

    def __init__(self, material, thickness=0.0):
        self.material = material
        self.thickness = thickness

    def _calculate_geometry(self, exterior):

        inside_direction = _get_inside_direction(exterior)

        interior = exterior.parallel_offset(self.thickness, side=inside_direction)

        if interior.type == 'MultiLineString':
            interior = shpl_geom.LinearRing(interior.geoms[0])
        else:
            interior = shpl_geom.LinearRing(interior)

        geometry = shpl_geom.Polygon(exterior)-shpl_geom.Polygon(interior)

        return interior, geometry

    def _centerline(self):

        inside_direction = _get_inside_direction(self.exterior)

        centerline = self.exterior.parallel_offset(self._thickness/2.0,
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


class CompositeLayer(Layer):
    pass #TODO: implement


class Reinforcement(Layer):
    """Local reinforcement structure
    
    Parameters
    ----------
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

    def __init__(self, material, thickness=0.0, limits=None):
        super().__init__(material, thickness)
        self._limits = np.array(limits)
    
    def _calculate_geometry(self, exterior):
        
        limited_box = shpl_geom.box(self._limits[0], exterior.bounds[1]*1.1,
                                    self._limits[1], exterior.bounds[3]*1.1)

        intersection = limited_box.intersection(exterior)

        side = _get_inside_direction(exterior)

        tmp_geometries = []

        tmp_interior = shpl_geom.Polygon(exterior)

        for ageo in intersection.geoms:
            tmp_geometry = _create_offset_box(ageo, self.thickness, side, 
                                              symmetric=False)
            tmp_interior -= tmp_geometry
            tmp_geometries.append(tmp_geometry)

        geometry = shpl_geom.GeometryCollection(tmp_geometries)

        interior = shpl_geom.LinearRing(tmp_interior.exterior)
        
        return _refine_interior(interior), geometry

    def exportgeometry(self, refpoint=np.zeros(2)):

        return [geom2array(geo, refpoint)._replace(material=self.material) for geo in self.geometry.geoms]

class MassAnalysis:
    """Analyse Mass of Structure
    
    Parameters
    ----------
    sectionbase
        section definiton to be analyzed
    
    """

    def __init__(self, sectionbase):
        self.secbase = sectionbase

    @property
    def massproperties(self):
        mass = 0.0
        cg = np.zeros(2)
            
        for feature in self.secbase._features:
            cur_cg, cur_mass = feature.massproperties

            mass += cur_mass
            cg += cur_mass * np.array(cur_cg)
        
        return cg/mass, mass


def _oderside(side):
    if side == 'left':
        return 'right'
    else:
        return 'left'


def rework_svg(svg:str, width:float, height:float=100.0, stroke_width:float=None)->str:
    """Prettify svgs generated by shapely
    
    Parameters
    ----------
    svg : str
        svg definitiom
    width : float
        with of returned svg
    height : int, optional
        height of returned svg (the default is 100)
    stroke_width : float, optional
        width of lines in svg (the default is None, which [default_description])
    
    Raises
    ------
    Exception
        raised when input svg string is invalid
    
    Returns
    -------
    str
        prettified svg string
    """

    import re
    
    if stroke_width is None:
        
        search_res = re.search(r'viewBox="((?:-?\d+.\d+\s*){4})"', svg)
        
        if search_res is None:
            raise Exception('Invalid SVG!')
            
        bounds = [float(val) for val in search_res.groups()[0].split()]
        Δy = bounds[3]-bounds[1]
        Δx = bounds[2]-bounds[0]
        
        stroke_width = min(Δx,Δy)/100
    
    svg = re.sub(r'width="\d+\.\d+"',
                 'width="{:.1f}"'.format(width),
                 svg,
                 count=1)
    
    svg = re.sub(r'height="\d+\.\d+"',
                 'heigth="{:.1f}"'.format(height),
                 svg,
                 count=1)
    
    svg = re.sub(r'stroke-width="\d+\.\d+"',
                 'stroke-width="{:f}"'.format(stroke_width),
                 svg)

    return svg


Structure = namedtuple('Structure', ['exterior','interiors','material'])


def exportstructure(geometry, material, refpoint=np.zeros(2)):
    """Helper function to export polygon geometry from shapely to numpy arrays based format
    
    Parameters
    ----------
    geometry : shapely.Polygon
        Geometrie to be exported (must be a Polygon)
    refpoint : np.array, optional
        reference point (2D), will be substracted from coordinates 
        (the default is np.zeros(2)
    
    Raises
    ------
    ValueError
        Wrong geometry type given..
    
    Returns
    -------
    ArrayStruc
        Custom geometry and material container
    """

    def coordsclockwise(linearring):
        if cga.signed_area(linearring) < 0.0:
            return np.array(linearring.coords)
        else:
            return np.array(linearring.coords)[::-1]

    if geometry.type in ('GeometryCollection', 'MultiPolygon'):
        res = []

        for geo in geometry.geoms:
            res.extend(exportstructure(geo, material, refpoint))
        
        return res

    elif geometry.type != 'Polygon':
        raise ValueError('Geometry must be of type \'Polygon\' or \'GeometryCollection\', not {}'.format(geometry.type))

    exterior = coordsclockwise(geometry.exterior) - refpoint.flat

    interiors = [coordsclockwise(interior) - refpoint.flat for interior in geometry.interiors]

    return [Structure(exterior, interiors, material)]
    