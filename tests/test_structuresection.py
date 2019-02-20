import pytest
import numpy as np
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

    # create a massanalysis
    massana = MassAnalysis(secbase)
    res = massana.massproperties

    assert np.isclose(res[1], 0.011681585464705432)
    assert all(np.isclose(res[0], [0.49406408, 0.02683857]))
    # remove I-spar and add boxspar
    secbase.pop()
    boxspar = section.BoxSpar(amat, 0.5, 0.2, 0.01, 0.01)
    secbase.append(boxspar)

    # remove first shell layer
    secbase.remove(outerlayer)

    # remove first shell layer second time -> Exception
    with pytest.raises(Exception):
        secbase.remove(outerlayer)

    # add again
    secbase.insert(0, outerlayer)

    # change thickness of outerlayer
    outerlayer.thickness = outerlayer.thickness*2

    # analyse mass again
    res = massana.massproperties

    assert np.isclose(res[1], 0.010542072442709407)
    assert all(np.isclose(res[0], [0.48091987, 0.02773624]))
