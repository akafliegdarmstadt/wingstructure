import pytest
from wingstructure.structure import section, material, MassAnalysis


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

    # create I-spar and add to secbase
    ispar = section.ISpar(amat, 0.5, 0.2, 0.01, 0.5, 0.01)
    secbase.append(ispar)

    # remove I-spar and add boxspar
    secbase.pop()
    boxspar = section.BoxSpar(amat, 0.5, 0.2, 0.01, 0.01)
    secbase.append(boxspar)

    # remove first shell layer
    secbase.remove(outerlayer)

    # create a massanalysis
    massana = MassAnalysis(secbase)
    res = massana.massproperties

    assert res is not None
