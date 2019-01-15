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
    
    outerlayer = section.Layer(secbase, amat, thickness=0.001)
    innerlayer = section.Layer(outerlayer, amat, thickness=0.001)

    ispar = section.ISpar(innerlayer, amat, 0.5, 0.2, 0.01,
                                0.5, 0.01)

    boxspar = section.BoxSpar(innerlayer, amat, 0.5, 0.2, 0.01,
                                    0.01)
def test_MassAnalysis(airfoilcoords):
    # create material
    amat = material.IsotropicMaterial(1.225, 210e3, 50e3)

    # create sectionbase instance
    secbase = section.SectionBase(airfoilcoords)
    
    outerlayer = section.Layer(secbase, amat, thickness=0.001)
    innerlayer = section.Layer(outerlayer, amat, thickness=0.001)

    ispar = section.ISpar(innerlayer, amat, 0.5, 0.2, 0.01,
                                0.5, 0.01)

    boxspar = section.BoxSpar(innerlayer, amat, 0.5, 0.2, 0.01,
                                    0.01)

    massana = MassAnalysis(secbase)
    cg, mass = massana.massproperties

    #TODO: write a meaningful test
    assert mass != None
