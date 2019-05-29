import numpy as np

from shapely.algorithms import cga
from shapely import geometry as shpl_geom
from shapely import ops

def _get_inside_direction(linearring):
    """Gets the inside direction for parallel offset (left or right)
    from signed area of geometry"""
    
    if cga.signed_area(linearring) > 0:
        return 'left'
    else:
        return 'right'


def _create_offset_box(line, thickness, side, bevel=0.0,
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


def _updatedecorator(method):
    def wrappedcall(self, *args, **kwargs):
        method(self, *args, **kwargs)
        self._update()
    return wrappedcall


def _refine_interior(interior):
    from numpy.linalg import norm

    if interior.type == 'MultiLineString':
        interior = shpl_geom.LinearRing(interior.geoms[0])
    else:
        interior = shpl_geom.LinearRing(interior)

    # find bugs in interior and remove them
    pts = np.array(interior)[:-1, :]
    del_pts = []

    for ii in range(len(pts)):
        vec1 = pts[ii]-pts[ii-1]
        vec2 = pts[(ii+1) % len(pts)]-pts[ii]

        if(norm(vec1/norm(vec1)+vec2/norm(vec2)) < 1e-3):
            del_pts.append(ii)

    pts2 = np.delete(pts, del_pts, axis=0)

    return shpl_geom.LinearRing(pts2)