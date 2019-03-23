Basic examples
==============

Creating a Crosssection
-----------------------

At first you have to load the geometry and material properties of the crosssection.
Geometry, youngs moduli, shear moduli and thickness have to be tuples with a length of 3.
Each tuple contains the geometries or material properties of left side, center and right side of the crosssection.
Create the crosssection with
generate_testsection = lambda : beam.Crosssection(geometry_section, youngsmoduli_section, shearmoduli_section, thickness_section)

note: in the future just one file for the geometry and material properties of all the crosssections of the wing

Importing forces and moments
----------------------------

wie geschieht das jetzt?
optimal wäre, wenn man einfach die Kräfte und Momente an jeder crosssection importieren würde

Calculate shearflow
-------------------

At first the shearforce has to be transformed into the principal axis:
Q_new = beam.transform(np.array((Q)), testsections[0].Θ)

Then you can calculate the shearflow:
shearflow = testsections[0].shearflow_total(Q_new, T)

Calculate twist of crosssection
-------------------------------

testsections.twist(Q,T)



