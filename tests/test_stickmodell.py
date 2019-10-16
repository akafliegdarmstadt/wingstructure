import pytest
import numpy as np


@pytest.fixture
def flatwing():
    from wingstructure.data.wing import Wing
    from wingstructure.wingloads import FlatWing

    wing = Wing()
    wing.append()
    wing.append(pos=(0.0, 1.0, 1.0))

    flatwing = FlatWing(wing)

    return flatwing


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


def test_transform_forces(flatwing):
    from wingstructure.structure.stickmodel import transform_forces

    loads = np.array([
        [0,1,0,0,0,1,5]
    ])

    tloads = transform_forces(flatwing, loads)

    assert np.isclose(
            tloads[0, :3], 
            [0, 1/np.sqrt(2), 1/np.sqrt(2)]
        ).all()


def test_transform_forces_with_rotation(flatwing):
    from wingstructure.structure.stickmodel import transform_forces 

    loads = np.array([
        [0.0,  1, 0, 0,  1, 0, 5],
        [0.0,  1, 0, 0, -1, 0, 4000]
    ])

    tloads = transform_forces(flatwing, loads, rotate=True)

    # check that inline method works as expected
    transform_forces(flatwing, loads, rotate=True, inline=True)
    assert np.isclose(tloads, loads).all()

    # check correct transformation of attack point
    assert np.isclose(
            tloads[0, :3], 
            [0, 1/np.sqrt(2), 1/np.sqrt(2)]
    ).all()

    # check transformation/rotation of force vector
    assert np.isclose(
        tloads[0 ,3:-1],
        [0, 1/np.sqrt(2), 1/np.sqrt(2)]
    ).all()

    assert np.isclose(
        tloads[1 ,3:-1],
        [0, -1/np.sqrt(2), -1/np.sqrt(2)]
    ).all()


def test_transformmoments(flatwing):
    from wingstructure.structure.stickmodel import transform_moments

    # grid points
    ys = np.array([-np.sqrt(2), -1.0, 0, 1.0, np.sqrt(2)])

    # moments on left and right side
    moments = np.array([
        [1.0, 2.0, 2.0, 0],
        [1.0, -2.0, 2.0, 3]
    ])

    # do transformation
    tmoments = transform_moments(flatwing, moments, ys)

    # check inline option
    transform_moments(flatwing, moments, ys, inline=True)

    assert np.isclose(tmoments, moments).all()

    # check for equalitiy of absolute values for left and right side
    assert np.isclose(*np.abs(tmoments[:,:-1])).all()

    # check left side
    assert np.isclose(tmoments[0,:-1], [1.0, np.sqrt(8), 0.0]).all()


def test_getnodes():
    from wingstructure.data.wing import Wing
    from wingstructure.structure.stickmodel import get_nodes

    ys = np.array([0.0, 0.5, 1.0, 1.2, 2.0])

    awing = Wing()

    awing.append()
    awing.append(pos=(0.0, 1.0, 0.0))
    awing.append(pos=(0.0, 1.0, 1.0))

    res = get_nodes(awing, ys, chordpos=0.0)

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


def test_solve():
    """
    test equilibrium solver
    """
    from wingstructure.structure.stickmodel import solve_equilibrium

    # straight bar with two segments, along y-direction
    nodes = np.array([[0,0,0], [0,1,0], [0,2,0]])

    # two forces acting on bar
    # first one: middle of first segment, magnitue sqrt(2), angle 45°
    # second one: middle of second segment, magnitue 2 in z direction
    forces = np.array([
        [0, 0.5, 0, 0, 1, 1, 0],
        [0, 1.5, 0, 0, 0, 2, 1]
    ])

    # two discrete moments acting on bar
    # first one: magnitude 1 at first segment (index 0)
    # second one: magnitude 2 at second segment (index 1)
    moments = np.array([
        [0, 0, 1, 0],
        [0, 0, 2, 1]
    ])

    # actual solving, setting forces and moments at node 2 to zero (prescribed)
    sol = solve_equilibrium(nodes, forces, moments, prescribed={2:np.zeros(6)})

    # check with manual calculation
    assert np.isclose(sol, [
        [0, 1, 3, 3.5, 0, 3],
        [0, 0, 2, 1, 0, 2],
        [0, 0, 0, 0, 0, 0]
    ]).all()


def test_2spar():
    from wingstructure.structure.stickmodel import internalloads2spar

    loads = np.array([
        [0, 0, 1, 0, 1, 0],
        [0, 0, 1, 0, 1, 0]
    ], dtype=float)

    sparnodes = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 1.0]
    ])

    tloads = internalloads2spar(loads, sparnodes)

    val = np.sqrt(2)/2

    assert np.isclose(tloads[0,:],
        [0, val, -val, 0, val, val]).all()