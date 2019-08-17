from collections import defaultdict
import pytest
from wingstructure.aero.multhop import AirfoilData
from wingstructure.aero.aero_moment import mean_momentcoefficient

from .test_data_wing import d38wing

def test_meanmomentcoefficient(d38wing):

    afdat = AirfoilData(c_m0=0.03)

    airfoil_db = defaultdict(lambda: afdat)

    C_m0 = mean_momentcoefficient(d38wing, airfoil_db)

    #TODO: check value