from shapely.geometry import LinearRing, Polygon, box, LineString
from shapely.ops import cascaded_union, linemerge
from shapely.algorithms import cga
import numpy as np

def getExtrapoledLine(p1,p2, length):
    'Creates a line extrapoled in p1->p2 direction'
    line = LineString(p1,p2)
    EXTRAPOL_RATIO = length/line.strength
    a = p1
    b = (p1[0]+EXTRAPOL_RATIO*(p2[0]-p1[0]), p1[1]+EXTRAPOL_RATIO*(p2[1]-p1[1]) )
    return LineString([a,b])

class WingSection(object):

    def __init__(self, wing, airfoil_coords_db, skins, spar, trailing_web):
        
        self.wing_definition = wing
        self.coords_db = airfoil_coords_db
        self.geometries = dict()
        self.skins = skins
        self.spar = spar
        self.trailing_web = trailing_web
        self.span_pos = 0.0
        self.mass = 0.0
        self.area = 0.0
        self.cg = np.array((0.0, 0.0))
        
        self.inter = None
        self.intera = []
        
    def update_geometry(self, span_pos):
         
        chord_length = self.wing_definition.chord_at(span_pos)
        coords = self.calculate_airfoil_coords(span_pos) * chord_length
        flap_depth = self.wing_definition.get_flap_depth(span_pos)
        
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
        
        #  ceate skins
        for ii, skin in enumerate(self.skins):            
            print('Iteration {}: type={}'.format(ii, type(self.interior)))
            
            if 'span-bounds' in skin.keys():
                if span_pos < skin['span-bounds'][0] or span_pos >= skin['span-bounds'][1]:
                    continue
            
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
                

                if flap_depth > 0.0:
                    print('cut flap')
                    tmp_poly = tmp_poly.difference(box(chord_length*(1-flap_depth), tmp_poly.bounds[1], tmp_poly.bounds[2], tmp_poly.bounds[3]))
                     
                self.geometries['skin_{}'.format(ii)] = tmp_poly                
                self.interior = tmp_interior                
                self.area += tmp_poly.area                
                tmp_mass = tmp_poly.area * skin['rho']
                self.mass += tmp_mass
                
                cg_list.append({'mass': tmp_mass, 'point': np.array(tmp_poly.centroid.coords[0])})
                    
        # create spar
        bounds = chord_length*np.array((self.spar['position']-self.spar['flange_width']/2, self.spar['position']+self.spar['flange_width']))
        self._create_reinforcements(bounds, self.spar['flange_thickness'], self.spar['flange_rho'], cg_list, name='flange')
        # web
        bounds = chord_length*(self.spar['web_position'] + np.array((-1,1))*(self.spar['web_thickness']))
        abox = box(bounds[0], self.interior.bounds[1], bounds[1], self.interior.bounds[3])
        web = Polygon(self.interior).intersection(abox)
        tmp_mass = web.area * self.spar['web_rho']
        self.geometries['web'] = web
        self.area += web.area
        self.mass += tmp_mass
        cg_list.append({'mass': tmp_mass, 'point': np.array(web.centroid.coords[0])})
        
        # create trailing web
        if flap_depth > 0.0:
            wing_end = chord_length*(1-flap_depth)
            web_end = wing_end - self.trailing_web['offset']
            web_start = web_end - self.trailing_web['thickness']
            abox = box(web_start, self.interior.bounds[1], web_end, self.interior.bounds[3])
            tmp_poly = abox.intersection(Polygon(self.interior))
            tmp_mass = tmp_poly.area * self.trailing_web['rho']    
            self.geometries['web_trailing'] = tmp_poly
            self.intera.append(tmp_poly)
            print('ajsdföjsad:{}'.format(tmp_poly.type))
            self.area += tmp_poly.area
            self.mass += tmp_mass
            abox = box(web_start, *self.interior.bounds[1:])
            self.interior = self.interior.difference(abox)
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
    
    def _create_reinforcements(self, bounds, thickness, rho, cg_list, bevel=0.0, name='noname'):
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
            if bevel > 0.0: 
                tmp_poly = self._create_offset_box(geom, thickness, side, bevel=bevel)
            elif bevel < 0.0:
                raise Exception('Negativ bevel angle not alowed')
            else:
                tmp_poly = self._create_offset_box(geom, thickness, side)                                        
            self.geometries[name+'_{}'.format(part)] = tmp_poly
            tmp_interior = Polygon(tmp_interior).difference(tmp_poly).exterior                                        
            self.area += tmp_poly.area                    
            tmp_mass = tmp_poly.area * rho
            self.mass += tmp_mass                
            cg_list.append({'mass': tmp_mass, 'point': np.array(tmp_poly.centroid.coords[0])})
            
        self.interior = tmp_interior
        
    def _create_offset_box(self, line, thickness, side, bevel=0.0):
        
        offsetline = line.parallel_offset(thickness, side=side)
        
        if bevel > 0.0:
            
            raise Exception('Beveling not yet implemented')
        
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
                
                if len(geometry.interiors) > 0:
                    coords_int = np.array(geometry.interiors[0].coords)
                    plt.plot(coords_int[:,0], coords_int[:,1], 'black')
                
            elif key.startswith('part_skin_') or key.startswith('flange') or key.startswith('web'):
                coords_ext = np.array(geometry.exterior.coords)
                plt.plot(coords_ext[:,0], coords_ext[:,1], 'black')
                
            plt.plot(self.cg[0], self.cg[1], 'r+')
            
            plt.axis('equal')

