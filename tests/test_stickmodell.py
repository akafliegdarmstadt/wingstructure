import pytest
import numpy as np


def test_lineloadresultants():
    from wingstructure.structure.stickmodel import calc_lineloadresultants

    loads = calc_lineloadresultants(
        (-2.0, -1.0, 1.0, 2.0, 3.0, 4.0, 5.0), 
        (1.0, 1.0, -1.0, 1.0, 2.0, 0.0, 0.0)
    )

    # check that all not involved values are zero
    assert (loads[:,0]==0.0).all()
    assert (loads[:,2]==0.0).all()
    assert (loads[:,3]==0.0).all()
    assert (loads[:,4]==0.0).all()

    # check positions of resulting forces
    assert np.isclose(loads[:,1], (-1.5, -2/3, 2/3, 1+1/6, 2-1/6, 2+5/9, 3+1/3)).all()

    # check values of resulting forces
    assert np.isclose(loads[:,5], (1, 0.5, -0.5, -0.25, 0.25, 1.5, 1)).all()

    # check segment assigment
    assert np.all(loads[:,-1] == [0, 1, 1, 2, 2 ,3, 4])


def test_discretemoments():
    from wingstructure.structure.stickmodel import calc_discretemoments

    ys = np.array((-2.0, -1.0, 1.0, 2.0, 3.0, 4.0, 5.0))
    m = (1.0, 1.0, -1.0, 1.0, 2.0, 0.0, 0.0)

    moments = calc_discretemoments(ys, m, axis=0)

    assert np.isclose(
        moments[:, 0],
        [1, 0, 0, 1.5, 1, 0]
    ).all()



def test_transformload():
    from wingstructure.data.wing import Wing, FlatWing
    from wingstructure.structure.stickmodel import transform_loads
    
    wing = Wing()
    wing.append()
    wing.append(pos=(0.0, 1.0, 1.0))

    flatwing = FlatWing(wing)

    loads = np.array([
        [0,1,0,0,0,1,0]
    ])

    tloads = transform_loads(flatwing, loads)

    assert np.isclose(
            tloads[0, :3], 
            [0, 1/np.sqrt(2), 1/np.sqrt(2)]
        ).all()

def test_transformload_with_rotation():
    from wingstructure.data.wing import Wing, FlatWing
    from wingstructure.structure.stickmodel import transform_loads 

    wing = Wing()
    wing.append()
    wing.append(pos=(0.0, 1.0, 1.0))

    flatwing = FlatWing(wing)

    loads = np.array([
        [0,1,0,0,0,1,0]
    ])

    tloads = transform_loads(flatwing, loads, rotate=True)

    # check correct transformation of attack point
    assert np.isclose(
            tloads[0, :3], 
            [0, 1/np.sqrt(2), 1/np.sqrt(2)]
        ).all()

    # check transformation/rotation of force vector
    assert np.isclose(
        tloads[0 ,3:-1],
        [0, -1/np.sqrt(2), 1/np.sqrt(2)]
    ).all()

def test_getnodes():
    from wingstructure.data.wing import Wing
    from wingstructure.structure.stickmodel import get_nodes

    ys = np.array([0.0, 0.5, 1.0, 1.2, 2.0])

    awing = Wing()

    awing.append()
    awing.append(pos=(0.0, 1.0, 0.0))
    awing.append(pos=(0.0, 1.0, 1.0))

    res = get_nodes(awing, ys)

    assert np.isclose(
        res,
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.2],
            [0.0, 1.0, 1.0]
        ]
    ).all()


def test_solve2Ds():
    from wingstructure.structure.stickmodel import solve_equilibrium

    nodes = np.array([[0,0,0], [0,1,0], [0,2,0]])
    forces = np.array([
        [0, 0.5, 0, 0, 1, 1, 0],
        [0, 1.5, 0, 0, 0, 2, 1]
    ])

    moments = np.array([
        [0, 0, 1, 0],
        [0, 0, 2, 1]
    ])

    sol = solve_equilibrium(nodes, forces, moments, free_node=2)

    assert np.isclose(sol, [
        [0, 1, 3, 3.5, 0, 3],
        [0, 0, 2, 1, 0, 2],
        [0, 0, 0, 0, 0, 0]
    ]).all()
