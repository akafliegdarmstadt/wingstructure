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