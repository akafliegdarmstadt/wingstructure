====================
Structure Definition
====================

This document descibes methods for mass estimation of the wing structure. As first step numpy and the wingstructure.structure module are imported:

.. code-block:: python

   import numpy as np
   from wingstructure import structure

With the modules loaded airfoil coordinates 
are loaded and scaled. Afterwards a *SectionBase* instance is created:

.. code-block:: python

   coords = np.loadtxt('ah93157.dat', skiprows=1) * 1.2
   sectionbase = structure.SectionBase(coords)

Furthermore we need some materials defined:

.. code-block:: python

   import collections
   Material = collections.namedtuple('Material', ['ρ'])

   carbonfabric = Material(1.225)
   foam = Material(1.225)
   sandwich = Material(1.225)

Starting from the *SectionBase* object the structure will be generated beginning with the outside surface. Elements that can be
added are *Layer*, *Reinforcement*, *BoxSpar* and *ISpar*. A *Layer* covers the full surface and the only parameters
are material and thickness. In the following a sandwich shell is defined:

.. code-block:: python

   outerlayer = structure.Layer(sectionbase, carbonfabric, 5e-4)
   core = structure.Layer(outerlayer, foam, 1e-2)
   innerlayer = structure.Layer(core, carbonfabric, 5e-4)

The previous structure element is always passed as a first argument. Now a spar can be defined. Here we choose a
I-spar instead of the box-spar.

.. code-block:: python

   spar = structure.ISpar(parent=innerlayer, 
                          material={'flange': carbonfabric, 'web': sandwich},
                          midpos=0.45,
                          flangewidth=0.2,
                          flangethickness=0.03,
                          webpos=0.5,
                          webthickness=0.02)
   
After the definiton of the structure, it can be analysed. One currently available analysis is the mass anlysis:

.. code-block:: python

   massana = structure.MassAnalysis(spar)
   cg, mass = massana.massproperties()

Here cg is the center of gravity and mass the section mass per unit wingspan.

