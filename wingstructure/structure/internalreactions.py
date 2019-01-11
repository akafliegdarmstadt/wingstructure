import numpy as np


def calc_lineloadresultants(ys, q):
    """Calculate resultants of line loads acting on individual segments
    """
    # calculate element lengths
    Δys = np.diff(ys)
    
    # initialize arrays for
    
    # force resultants
    Q = np.zeros_like(Δys)
    
    # resultants attack point
    y_res = np.zeros_like(Δys)
    
    # iterate over parts of wing
    for i in range(1,len(y_res)):
        
        # trapez rule to get resultant
        Q[i-1] = Δys[i-1] * (q[i]+q[i-1])/2
        # center of trapez as attack point of resultant
        y_res[i-1] = ys[i] + np.abs(Δys[i-1])/3 * np.abs((q[i]+2*q[i-1]) / (q[i]+q[i-1]))
        
        if ys[i] > y_res[i-1]  or ys[i-1] < y_res[i-1]  :
            print('Achtung Achtung')
        
    return y_res, Q


def combine_loads( loadslist ):
    loads = np.vstack(loadslist)
    
    # sort by y coord
    loads = loads[np.flip(loads[:,1].argsort())]
    
    return loads


def calc_moments(coords, discrete_loads):
    """Calculates cummulated forces and moments
    """
    
    # number of grid points
    n = coords.shape[0]
    
    # split up loads
    coords_Qs = discrete_loads[:,:3]
    Qs = discrete_loads[:,3:]
    
    # initialize arrays with zeros for cumulated forces and moments
    Q_k = np.zeros((n, 3))
    M = np.zeros_like(Q_k)
    
    # iterate over grid points
    j = 0
    for i in range(1, n):

        # initialize cumulated force with previous one
        Q_k[i,:] += Q_k[i-1,:]

        # initialize section moment with
        # previous sections moment + moment of previous cumulated force
        a = (coords[i-1] - coords[i]).T
        M[i,:] += M[i-1,:] + np.cross(Q_k[i-1, :], a)

        while np.abs(coords_Qs[j, 1]) > np.abs(coords[i, 1]):
            # sum up forces
            Q_k[i,:] += Qs[j,:]

            # sum up resulting moments
            ## lever
            a = (coords_Qs[j, :] - coords[i, :]).T

            ## cross product for moment
            M[i,:] += np.cross(Qs[j,:],a)

            j += 1

            # break if j exceeds resultants size
            if j >= Qs.shape[0]: break
    
    return Q_k, M