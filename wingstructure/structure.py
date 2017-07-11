from shapely.geometry import LinearRing, Polygon, box, LineString
from shapely.ops import cascaded_union, linemerge
from shapely.algorithms import cga
import numpy as np


class WingSection(object):

    def __init__(self, wing, airfoil_coords_db, skins, spar):
        
        self.wing_definition = wing
        self.coords_db = airfoil_coords_db
        self.geometries = dict()
        self.skins = skins
        self.spar = spar
        self.span_pos = 0.0
        self.mass = 0.0
        self.area = 0.0
        self.cg = np.array((0.0, 0.0))
        
        self.inter = None
        
    def update_geometry(self, span_pos):
         
        chord_length = self.wing_definition.chord_at(span_pos)
        coords = self.calculate_airfoil_coords(span_pos) * chord_length
        
        # TODO: turn profil according section alpha
        cg_list = []
        
        # Reset
        self.geometries = dict()
        self.mass = 0.0
        self.area = 0.0
        self.cg = np.array((0.0, 0.0))
        
        # Recalculate
        outline = LinearRing(coords)
        self.interior = LinearRing(outline)
        
        #  Iterate over skins
        for ii, skin in enumerate(self.skins):            
            print('Iteration {}: type={}'.format(ii, type(self.interior)))
                        
            if 'chord-bounds' in skin.keys():
                chord_bounds = chord_length*np.array(skin['chord-bounds'])
                
                self._create_reinforcements(chord_bounds, skin['thickness'], skin['rho'], cg_list, name = 'part_skin_{}'.format(ii))
            else:
                print('create skin_{}'.format(ii))
                tmp_interior = self.interior.parallel_offset(skin['thickness'], side=self._get_inside_direction(self.interior))
                
                if tmp_interior.type == 'MultiLineString':
                    tmp_interior = LinearRing(tmp_interior.geoms[0])
                else:
                    tmp_interior = LinearRing(tmp_interior)
                
                tmp_poly = Polygon(self.interior).difference(Polygon(tmp_interior))            
                self.geometries['skin_{}'.format(ii)] = tmp_poly                
                self.interior = tmp_interior                
                self.area += tmp_poly.area                
                tmp_mass = tmp_poly.area * skin['rho']
                self.mass += tmp_mass
                
                cg_list.append({'mass': tmp_mass, 'point': np.array(tmp_poly.centroid.coords[0])})
                    
        # create spar
        bounds = chord_length*np.array((self.spar['position']-self.spar['width']/2, self.spar['position']+self.spar['width']))
        self._create_reinforcements(bounds, self.spar['thickness'], self.spar['rho'], cg_list, name='flange')
        
        # Calculate Center of Gravity
        cg_sum = np.zeros(2)
        for cg_point in cg_list:
            cg_sum += cg_point['point'] * cg_point['mass']
        self.cg = cg_sum/self.mass
    
    def calculate_airfoil_coords(self, span_pos):
        
        airfoils = self.wing_definition.airfoils_at(span_pos)
        
        if len(airfoils) == 1:
            return self.coords_db[airfoils[0]]
        else:
            from .geometry import Airfoil
            coords1 = self.coords_db[airfoils[0]]
            coords2 = self.coords_db[airfoils[1]]
            
            af1 = Airfoil(coords1)
            af2 = Airfoil(coords2)
            
            beta = airfoils['beta']
            
            af3 = af1.interpolate(af2, beta)
            
            return af3.coords
    
    def _create_reinforcements(self, bounds, thickness, rho, cg_list, name='noname'):
        print(' create reinforcement {}'.format(name))
        print('  bounds: {}\n   thickness:{}'.format(bounds, thickness))
        # Create Box for intersecting with interior
        abox = box(bounds[0], self.interior.bounds[1]*1.1, bounds[1], self.interior.bounds[3]*1.1)
        intersection = abox.intersection(self.interior) 
        # create copy of interior
        tmp_interior = LinearRing(self.interior)
        # determin inside direction for Offsetting
        side = self._get_inside_direction(self.interior)
        
        self.inter=intersection
        for part, geom in zip(('upper','lower'),intersection.geoms):
            print('create '+name+'_{}'.format(part))
            tmp_poly = self._create_offset_box(geom, thickness, side)                                        
            self.geometries[name+'_{}'.format(part)] = tmp_poly
            tmp_interior = Polygon(tmp_interior).difference(tmp_poly).exterior                                        
            self.area += tmp_poly.area                    
            tmp_mass = tmp_poly.area * rho
            self.mass += tmp_mass                
            cg_list.append({'mass': tmp_mass, 'point': np.array(tmp_poly.centroid.coords[0])})
            
        self.interior = tmp_interior
        
    def _create_offset_box(self, line, thickness, side):
        
        offsetline = line.parallel_offset(thickness, side=side)
        
        if side == 'left':
            connect1 = LineString((line.coords[-1],offsetline.coords[-1]))
            connect2 = LineString((line.coords[0], offsetline.coords[0]))
            return Polygon(linemerge((line,connect1,offsetline,connect2)))
        else:
            connect1 = LineString((line.coords[-1],offsetline.coords[0]))
            connect2 = LineString((line.coords[0], offsetline.coords[-1]))
            return Polygon(linemerge((offsetline,connect1,line,connect2)))
        
    def _get_inside_direction(self, linearring):
        """Gets the inside direction for parallel offset (left or right) from signed area of geometry"""
        tmparray = np.array(linearring.coords)        
        top_id = np.argmax(tmparray[:,1])        
        bottom_id = np.argmin(tmparray[:,1])
        
        if cga.signed_area(linearring) > 0:
            return 'left'
        else:
            return 'right'
        
    def plot(self):
        
        from matplotlib import pyplot as plt
        
        for key, geometry in self.geometries.items():
            
            if key.startswith('skin_'):
                coords_ext = np.array(geometry.exterior.coords)
                plt.plot(coords_ext[:,0], coords_ext[:,1], 'black')
                coords_int = np.array(geometry.interiors[0].coords)
                plt.plot(coords_int[:,0], coords_int[:,1], 'black')
                
            elif key.startswith('part_skin_') or key.startswith('flange'):
                coords_ext = np.array(geometry.exterior.coords)
                plt.plot(coords_ext[:,0], coords_ext[:,1], 'black')
                
            plt.plot(self.cg[0], self.cg[1], 'r+')
            
            plt.axis('equal')

