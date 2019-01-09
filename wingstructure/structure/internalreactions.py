import numpy as np


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