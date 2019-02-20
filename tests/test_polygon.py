import pytest
import numpy as np

from wingstructure.structure import polygon


@pytest.fixture
def square():
    return np.array([(0.0,0.0),(0.0,1.0),(1.0,1.0),(1.0,0.0)])


def test_calcarea(square):
    # simple square area
    A = polygon.calcarea(square)

    assert np.isclose(A, 1.0)


def test_calcstaticmoments(square):
    # simple square
    Ss = polygon.calcstaticmoments(square)

    assert np.isclose(Ss, 1/2).all()


def test_intertiamoments(square):
    # simple square
    Is = polygon.calcinertiamoments(square-np.array([0.5,0.5]))

    assert np.isclose(Is, (1/12,1/12,0.0)).all()
    
    
@pytest.fixture
def rectangle():
    return np.array([(0.0,0.0),(0.0,2.0),(1.0,2.0),(1.0,0.0)])


def test_calcarea(rectangle):
    A = polygon.calcarea(rectangle)

    assert np.isclose(A, 2.0)


def test_intertiamoments(rectangle):
    Is = polygon.calcinertiamoments(rectangle-np.array([0.5,1.0]))

    assert np.isclose(Is, (2/3,1/6,0.0)).all()
    
    
@pytest.fixture
def triangle():
    return np.array([(0.0,0.0),(1.0,1.0),(1.0,0.0)])


def test_calcarea(triangle):
    A = polygon.calcarea(triangle)

    assert np.isclose(A, 0.5)


def test_intertiamoments(triangle):
    Is = polygon.calcinertiamoments(triangle-np.array([2/3,1/3]))

    assert np.isclose(Is, (1/36,1/36,-1/72)).all()
    
    
@pytest.fixture
def triangle():
    return np.array([(0.0,0.0),(1.5,1.0),(2.0,0.0)])


def test_calcarea(triangle):
    A = polygon.calcarea(triangle)

    assert np.isclose(A, 1.0)


def test_intertiamoments(triangle):
    Is = polygon.calcinertiamoments(triangle-np.array([7/6,1/3]))

    assert np.isclose(Is, (1/18,13/72,-1/36)).all()