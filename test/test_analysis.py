from wingstructure import geometry, LiftAndMomentAnalysis as LiftAnalysis, Airfoil

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize, integrate
import pandas as pd
import yaml

x = np.array([0, 0, 0.8926, 235.0122, 460.9768,
              585.6099])*1e-3 # in m
y = np.array([0.0, 325, 749.4176, 5646.6968,
              7141.6452, 7490.1669])*1e-3 # in m
c = np.array([671, 671, 671, 447, 224, 100])*1e-3 # in m
twists = np.zeros(len(x))
airfoils = ['dummy']*len(twists)

# Flügel Objekt erzeugen
wing = geometry.WingExt.create_from_planform(y, c, x, twists, airfoils)

# Flügelklappen hinzufügen
wing.set_flap('QR_innen', 0.749, 5.646, [0.15, 0.15])
wing.set_flap('QR_aussen', 5.646, 7.1416452, [0.25,0.35])

# Bremsklappe definieren
wing.set_airbrake(2,2+1.3)

# Höhenleitwerkshebelarm festlegen
x_hlw = 6.0 #m
x_sp = 0.2 #m

klappenstellungen = {
    'Kl+': {
        'QR_innen': [np.radians(20)]*2
    },
    'Kl-': {
        'QR_innen': [-np.radians(10)]*2
    }
}

massen=(450.0, 260.0, 450.0, 300.0)

# Profildatenbank anlegen
dummyairfoil = Airfoil(cm0 = 0.025) # nach JAR22 betragsmäßig der minimal zulässige Wert

airfoildb = {'dummy': dummyairfoil}

liftana = LiftAnalysis(wing, airfoildb)

alpha, c_l, c_ = liftana.calculate(1.5, klappenstellungen['Kl-'])

print(np.rad2deg(alpha))