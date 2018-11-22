===============
Getting Started
===============

This document gives a basic overview over the wingstructure package. Therefore a basic analysis is done.
As a first step the input data will be defined.

Input Data Definition
=====================

For the definition of data the *yaml* format is used. A data format designed to be human readable with a syntax 
similiar to *python*. For simplicity and clarity reasons we restrict ourselves 
to a subset called `strictyaml <https://hitchdev.com/strictyaml/>`_. 

The input yaml file for wingstructure contains a map as root structure. This root map
can have the keys *geometry*, *aerodynamic*, *configurations* and *mass*. Values for those
keys are described in the following. The data used is taken from the `D-38 <https://www.akaflieg.tu-darmstadt.de/d-38>`_ sailplane.
We will start with the geometry defintion:

Geometry Definition
"""""""""""""""""""

An important information is the wing geometry of our plane. To store it inside the *yaml*-file
a *wing* entry is added within the *geometry*-content:

.. code-block:: yaml

   geometry:
       wing:
           ....

The planform defined in the following table 

.. list-table:: D-38 Wing Planform
   :widths: 10 10 10 10 10 20
   :header-rows: 1

   * - x (offset)
     - y (span-position)
     - z (->dihedral)
     - chord length
     - twist
     - airfoil
   * - 0.0
     - 0.0
     - 0.0
     - 0.943
     - 0.0
     - FX 61-184
   * - 0.0
     - 4.5
     - 0.0
     - 0.754
     - -1.13
     - FX 61-184
   * - 0.134
     - 7.5
     - 0.0
     - 0.377
     - -3.86
     - FX 60-126

is transcribed to the yaml-format:

.. code-block:: yaml

   geometry:
       wing:
           sections:
              - pos:
                    x: 0.0
                    y: 0.0
                    z: 0.0
                chord: 0.943
                airfoil: FX 61-184
              - pos:
                    y: 4.5
                chord: 0.754
                twist: -1.13
                airfoil: FX 61-184
              - pos:
                    x: 0.134
                    y: 7.5
                chord: 0.377
                twist: -3.86
                airfoil: FX 60-126
            
The position (*pos*) can contain *x*, *y* und *z*-coordinates. Not listet values are set to zero.

Aerodynamic Data
""""""""""""""""
The next section treated is *aerodynamic*. Its content descibes the overall aerodynamic characteristics
of the sailplane and aerodynamic characteristics of the airfoils can be defined.

.. code-block:: yaml

   aerodynamic:
      c_dmin: 0.002
      airfoils:
        - name: 'FX 60-126'
          α_0: -0.15
          c_m0: 0.025
        - name: 'FX 61-184'
          α_0: -0.004
          c_m0: 0.25

Here a minimal drag for the whole airplane and characteristics for the two airfoils used are defined.
The airfoil characteristics include a zero lift angle of attack (α_0) and the constant moment coefficient
(cm_0) with regarding to the c/4-line.

.. tip::

   Additional information on input data definition can be found in :doc:`inputdata`.

Load Input Data
===============

You can use the standard *yaml* module to load the data. *wingstructure* has a specific function 
for that, which also validates the imported input data.

.. code-block:: python

   import wingstructure.data as wsdata

   inputdata = wsdata.loaddata('inputdata.yaml')



Aerodynamic analysis
====================
