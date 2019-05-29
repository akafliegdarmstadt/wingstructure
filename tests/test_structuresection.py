import pytest
import numpy as np
from wingstructure.structure import section, material
from wingstructure.structure.section import MassAnalysis


@pytest.fixture
def airfoilcoords():
    import numpy as np

    # load airfoil coordinates
    return np.loadtxt('docs/usage/FX 61-184.dat', skiprows=1, delimiter=',')


def test_structurecreation(airfoilcoords):
    # create material
    amat = material.IsotropicMaterial(1.225, 210e3, 50e3)

    # create sectionbase instance
    secbase = section.SectionBase(airfoilcoords)
    
    # create shell layers
    outerlayer = section.Layer(amat, thickness=0.001)
    innerlayer = section.Layer(amat, thickness=0.001)

    # add shell layers to secbase
    secbase.extend([outerlayer, innerlayer])


    # remove first shell layer
    secbase.remove(outerlayer)

    # remove first shell layer second time -> Exception
    with pytest.raises(Exception):
        secbase.remove(outerlayer)

    # add again
    secbase.insert(0, outerlayer)

    # change thickness of outerlayer
    outerlayer.thickness = outerlayer.thickness*2
