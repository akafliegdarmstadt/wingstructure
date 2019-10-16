import numpy as np


def calc_lineloadresultants(ys, q):
    """Calculate resultants of loads for piecewise linear load distributin
    
    Parameters
    ----------
    ys : numpy array, list
        grid points
    q : numpy array, list
        field values
    
    Returns
    -------
    array
        discrete resultant forces and coordinates [[x, y, z, Q_x, Q_y, Q_z], ...]
    """

    # calculate element lengths
    Δys = np.diff(ys)
    
    # initialize arrays for
    
    # force resultants
    Q = []
    
    # resultants attack point
    y_res = []

    # segments force resultant acts on
    segs = []
    
    # iterate over parts of wing
    for i in range(1,len(ys)):
        
        if q[i]==0 and q[i-1]==0:
            # Nothing to do..
            continue
        elif (q[i]>=0 and q[i-1]>=0) or (q[i]<=0 and q[i-1]<=0):
            # trapez rule to get resultant
            Q.append(Δys[i-1] * (q[i]+q[i-1])/2)
            # center of trapez as attack point of resultant
            y_res.append(ys[i-1] + np.abs(Δys[i-1])/3 * np.abs((q[i-1]+2*q[i]) / (q[i]+q[i-1])))
            segs.append(i-1)
        else:
            # sign changes from q[i-1] to q[i]
            # cannot be captured by single resultant within this section
            # -> resultants of the two triangles are used
            y_0 = ys[i-1] + Δys[i-1] * np.abs(q[i]/q[i-1]) / (1 + np.abs(q[i]/q[i-1]))

            #print(f'y_0 = {y_0, Δys, ys[i-1]}')

            # left triangle
            Q.append(0.5*(y_0-ys[i-1])*q[i-1])
            y_res.append(ys[i-1] + (y_0 - ys[i-1])/3.0)
            segs.append(i-1)

            # right triangle
            Q.append(0.5*(ys[i] - y_0)*q[i])
            y_res.append(ys[i] - (ys[i]-y_0)/3.0)
            segs.append(i-1)
    
    loads = np.zeros((len(y_res), 7))
    loads[:, 1] = y_res
    loads[:, -2] = Q
    loads[:, -1] = segs
    
    return loads


def calc_discretemoments(ys, m, axis=0):
    """Determine discrete moments from moment distribution
    
    Parameters
    ----------
    ys : array
        grid points
    m : array
        moment distribution values
    axis : int, optional
        axis the moments act, by default 0
    
    Returns
    -------
    array
        discrete moemnts
    """

    m = np.array(m)

    # calculate element lengths
    Δys = np.diff(ys)

    M = np.zeros((len(Δys), 4))
    M[:, axis] = 0.5 * (m[1:]+m[:-1]) * Δys

    M[:, 3] = np.arange(0,len(Δys))

    return M


def transform_forces(flatwing, forces, rotate=False, inline=False):
    """transform forces from flat to three dimensional wing
    
    Parameters
    ----------
    flatwing: FlatWing
        instance of flattend wing
    forces : array
        resultant loads array [[x, y, z, Q_x, N, Q_y]..]
    rotate : bool, optional
        switch for rotation of loads, by default False
    inline : bool, optional
        modify forces if True or return new transformed_forces array

    Returns
    -------
    array
        transformed_forces if inline=False

      ..caution::
        when using inline option make sure forces have floating point dtype 
    """

    ys = flatwing.ys
    dy = np.diff(ys)

    if not inline:
        forces = np.array(forces, copy=True, dtype=np.float)

    for j, load in enumerate(forces):
        y = load[1]

        # find last value smaller than y in ys
        i = np.searchsorted(ys, np.abs(y))

        # position in 3D
        pos1 = np.array(flatwing.basewing.sections[i-1].pos)
        pos2 = np.array(flatwing.basewing.sections[i].pos)

        # calculate relativ position between sections
        f = (np.abs(y)-ys[i-1])/dy[i-1]

        # interpolate position
        forces[j, :3] = pos1 + (pos2-pos1) * f

        if rotate:
            rotmat = np.eye(3)

            n0 = np.array((0,1,0))
            n1 = pos2-pos1
            n1 /= np.linalg.norm(n1)

            cosφ = np.dot(n0, n1)
            sinφ = np.sqrt(1-cosφ**2)

            rotmat[1:,1:] = [[cosφ, sinφ], [-sinφ, cosφ]]

            if y<0.0:
                rotmat = rotmat.T

            forces[j, 3:-1] = forces[j, 3:-1] @ rotmat

    if not inline:
        return forces


def transform_moments(flatwing, moments, ys, inline=False):
    """Transform moments into wing coordinate system
    
    Parameters
    ----------
    flatwing : FlatWing
        discription of flattend wing
    moments : array
        discrete moments [[Mx, My, Mz, segment],...]
    ys: array
        grid points
    inline: bool
        modify moments or return new array

    Returns
    -------
    array
        transformed moments, only if inline=False
    """

    from scipy.interpolate import interp1d

    if inline:
        transformed_moments = moments
    else:
        transformed_moments = np.copy(moments)

    positions = np.array([(sec.pos.y, sec.pos.z) for sec in flatwing.basewing.sections])

    normals = np.diff(positions, axis=0, append=0.0)
    normals[:-1,:] /= np.linalg.norm(normals[:-1,:])

    cosφ = interp1d(flatwing.ys, np.dot((1,0), normals.T), kind='previous')

    midy = ys[:-1] + np.diff(ys)/2

    for i, (_,_,_,seg_id) in enumerate(moments):
        y = midy[int(seg_id)]
        ccosφ = cosφ(abs(y))
        csinφ = np.sqrt(1-ccosφ**2)

        rotmat = np.eye(3)
        rotmat[1:,1:] = [[ccosφ, csinφ], [-csinφ, ccosφ]]

        if y<0.0:
            rotmat = rotmat.T
        transformed_moments[i,:-1] = transformed_moments[i,:-1] @ rotmat
    
    if not inline:
        return transformed_moments
   

def get_nodes(wing, ys, chordpos=0.25):
    """calculate grid points from wing
    
    Parameters
    ----------
    wing : Wing
        a object describing a wing
    ys : array
        grid points on flattend span
    chordpos: float
        relativ position in chordwise direction
    
    Returns
    -------
    array
        nodes [[x,y,z],...]
    """
    from numpy import diff, linalg
    from scipy.interpolate import interp1d

    pos = [(sec.pos.x + sec.chord*chordpos, sec.pos.y, sec.pos.z) \
                for sec in wing.sections]

    pos = np.array(pos)

    secy = np.cumsum(linalg.norm(np.diff(pos[:,1:], axis=0, prepend=0), axis=1))

    return interp1d(secy, pos, axis=0)(ys)


def solve_equilibrium(nodes, forces=np.zeros((1,7)), moments=np.zeros((1,4)), prescribed={0: np.zeros(6)}):
    """Solve static equilibrium in unbranched stick model
    
    Parameters
    ----------
    nodes : array
        coordinates of nodes: [[x, y, z], ...]
    forces : array
        forces with point of attack and segment: [[x, y, z, fx, fy, fz, seg], ...]
    moments: array
        discrete moments [[mx, my, mz, seg], ...]
    free_node : int
        designate node without loads
    
    Returns
    -------
    array
        internal loads at nodes [[Qx, Qy Qz, Mx, My, Mz], ...]
    """

    free_node = 2

    n = len(nodes)
    A = np.zeros([6*n,6*n])
    b = np.zeros(6*n)

    # equilibrium conditions (-x_0 + x_1 = 0)
    for i in range(6*(n-1)):
        A[i][i] = -1
        A[i][i+6] = 1

    # moments resulting from internal forces
    for i in range(n-1):
        lx, ly, lz = nodes[i+1] - nodes[i]
        idx = 6*i

        ## cross product lever x internal force
        A[idx+3][idx+7] = -lz
        A[idx+3][idx+8] = ly

        A[idx+4][idx+6] = lz
        A[idx+4][idx+8] = -lx

        A[idx+5][idx+7] = lx
        A[idx+5][idx+6] = -ly

    # force or moment free boundary conditions
    #for i in range(6):
    #    #idx = 6*free_node
    #    A[6*(n-1) + i][6*free_node + i] = 1

    for node, prescribed_values in prescribed.items():
        for i, prescribed_value in enumerate(prescribed_values):
            if prescribed_value is None:
                continue
            A[6*(n-1) + i][6*node + i] = 1
            b[6*(n-1) + i] = prescribed_value

    # right hand side / external loads
    for i in range(n-1):
        # forces
        for l in forces[forces[:,-1] == i]:
            # force
            b[6*i:6*i+3] += l[3:-1]
    
            # resulting moment (cross product)
            M = np.cross(l[:3] - nodes[i], l[3:-1]) 
            b[6*i+3:6*i+6] += M
        # moments
        for m in moments[moments[:,-1] == i]:
            b[6*i+3:6*i+6] += m[:-1]

    x = np.linalg.solve(A, -b)

    return np.reshape(x, (len(x)//6, 6))


def internalloads2spar(internalloads, sparnodes):
    """Transform internal loads to spar coordinate system
    
    Parameters
    ----------
    internalloads : array
        internal loads (global cartesian coordinate system)
    sparnodes : array
        coordinates of spar nodes
    
    Returns
    -------
    array
        transformed internal loads [[Qn, Q1 Q2, Mt, Mb1, Mb2], ...]
    """

    internalloads = np.copy(internalloads)
    
    for i, load in enumerate(internalloads):

        # get current spar normal

        ## if last node - use same direction vector as before
        if i == sparnodes.shape[0]-1:
            i -= 1
        
        ns = _get_normal(*sparnodes[i:i+2,:])
        
        # remove x component
        ns[0] = 0.0
        ns /= np.linalg.norm(ns)

        # collect all direction vectors
        n1 = np.array([1.0, 0.0, 0.0])
        n2 = ns
        n3 = np.cross(n2, n1)

        # build rotation matrix
        rotmat = np.vstack((n1,n2,n3)).T

        # do transformation
        internalloads[i,:3] = load[:3] @ rotmat
        internalloads[i, 3:] = load[3:] @ rotmat

    return internalloads

def _get_normal(p1, p2):
    n = p2-p1
    n0 = n/np.linalg.norm(n)

    return n0