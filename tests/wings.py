import pytest

@pytest.fixture
def d38wing():
    from wingstructure.data.wing import Wing

    awing = Wing()

    awing.append(chord=0.943, airfoil='FX 61-184')
    awing.append(chord=0.754, airfoil='FX 61-184', pos=(0.0, 4.51, 0.0))
    awing.append(chord=0.377, airfoil='FX 60-12', pos=(0.134, 7.5, 0.0))
    
    return awing

@pytest.fixture
def d43wing():
    from wingstructure.data.wing import Wing

    awing = Wing()

    awing.append((0,0,0), 1.12)
    awing.append((-68e-3,4.0,0.21996), 1.028)
    awing.append((51e-3,7.223,0.42874), 0.673)
    awing.append((0.1753225, 8.5, 0.6328), 0.454537)
    awing.append((0.245, 9.0, 0.7424), 0.36)
    
    return awing