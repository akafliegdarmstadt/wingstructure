from wingstructure.mass import Masspoint


def test_addition():
    mp1 = Masspoint(1.0, (1.0, 0.0, 0.0))
    mp2 = Masspoint(1.0, (-1.0, 0.0, 0.0))

    mpsum = mp1+mp2

    assert mpsum.mass == 2.0
    assert mpsum.point[0] == 0.0


def test_summation():
    mp1 = Masspoint(1.0, (1.0, 0.0, 0.0))
    mp2 = Masspoint(1.0, (-1.0, 0.0, 0.0))
    mp3 = Masspoint(1.0, (0.0, 0.0, 0.0))

    mpsum = sum((mp1, mp2, mp3))

    assert mpsum.mass == 3.0
    assert mpsum.point[0] == 0.0
