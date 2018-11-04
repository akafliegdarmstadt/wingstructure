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

Basic Yaml
----------

In *yaml* there are two concepts to structure data

 * maps
 * sequences

A *map* is *yaml*'s equivalent to a *python* *dictionary*:

.. code-block:: yaml

   a: 1.0
   b: text

The *sequences* are comparable to *python* *lists* or *tuples* and are marked with leading hyphens:

.. code-block:: yaml

   - item1
   - item2
   - 10.0

This structures can be abitrary nested. To indicate nesting indentation is used (spaces, no tabs):

.. code-block:: yaml

    a:
      - item1a
      - item2a
    b:
      - 5
      - 2.0
      - text
      - c: 2.0
        d: 3.0
        e: text2

Example Definition
^^^^^^^^^^^^^^^^^^

The input yaml file for wingstructure contains a map as root structure. This root map
can have the keys *geometry*, *aerodynamic*, *configurations* and *mass*. Values for those
keys are described in the following. The data used is taken from the `D-38 <https://www.akaflieg.tu-darmstadt.de>`_ sailplane.
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
                airfoil: FX 61-184
              - pos:
                    x: 0.134
                    y: 7.5
                chord: 0.377
                airfoil: FX 60-126
            
The position (*pos*) can contain *x*, *y* und *z*-coordinates. Not listet values are set to zero.

Load data in Python
-------------------

Object Orientet Interfaces
--------------------------


Aerodynamic analysis
====================
