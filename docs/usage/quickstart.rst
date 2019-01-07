===============
Getting Started
===============

This document is intended as starting point for wingstructure usage.
The wingstructure package has three main features:

 * **representation of wing geometry**

   Hereby the representation of a wing geometry in yaml files and a
   object oriented representation in python is ment.

 * **calculation of lift distribution**

   The distribution of lift and aerodynamic moments can
   be calculated using multhopps quadrature method.

 * **estimation of wing structure mass and structural analysis (sectionwise)**

   Starting from an airfoil *dat* file various
   structure elements (shell, spar) can be added. The mass and
   some structural quantities can be calculated.

A brief introduction for those functions is given in this document. For more detailed information have a
look at the Jupyter notebook `examples <https://github.com/akafliegdarmstadt/wingstructure/tree/master/examples>`_.
Some topics are already covered by specific documenatation pages
in deeper detail. Those are linked within the following sections.

Installation
============

Before wingstructure can be used, it has to be installed.
This can be done with *pip*:

.. code-block:: shell

   pip install https://github.com/akafliegdarmstadt/wingstructure/archive/master.zip


Basic examples
==============

Creating a wing
---------------

wingstructure can read geometry definition from *yaml* files, which simplifies
the reuse of those definitions in multiple scripts. The wing of `D-38 <https://www.akaflieg.tu-darmstadt.de/d-38/>`_
can be described with the following yaml file:

.. literalinclude:: D-38.yaml
   :caption: D-38.yaml
   :language: yaml
   :linenos:

To make the document extensible the wing geometry is stored inside
of a *wing* object, which itself is placed within *geometry*.
Starting with line 3 the wing planform is defined through sections.
The content should be rather self explanatory.

After the *sections* control surfaces are defined.

This definition can be loaded with any yaml parser, but it is recommended
to use wingstructure's :py:meth:`wingstructure.data.loaddata` function:

.. ipython:: ipython

   In [1]: from wingstructure import data

   In [2]: definition = data.loaddata('usage/D-38.yaml')

For a more convenient access to the data an
:py:class:`wingstructure.data.Wing` object is created.

.. ipython::

   In [3]: wing = data.Wing.create_from_dict(definition['geometry']['wing'])

   @savefig wingplot.png width=4in
   In [4]: wing.plot()

Besides plotting the wing object has some useful properties. For example
calculating the mean aerodynamic chord:

.. ipython::

   In [5]: wing.mac

Further information can be found in :py:class:`wingstructure.data.Wing`.

Calculating wing's lift distribution
------------------------------------

The wing is now defined and can be analysed. Therefore the analysis
submodule of wingstructure is imported. Within this submodule you find
an LiftAnalyis class, which will be used now.

.. ipython::
   :okwarning:
   
   In [6]: from wingstructure import analysis

   In [7]: liftana = analysis.LiftAnalysis(wing)

   In [8]: alpha, c_ls = liftana.calculate(1.0)

This calculates the angle of attack and local lift coefficients. Those
are plotted in the following. The associated span positions are stored
inside liftana (:py:attr:`wingstructure.analysis.LiftAnalysis.calc_ys`).

.. ipython::
   
   In [9]: from matplotlib import pyplot as plt
   
   In [10]: plt.figure();

   @savefig lift.png width=4in
   In [11]: plt.plot(liftana.calc_ys, c_ls)
      ....: plt.xlabel('span position')
      ....: plt.ylabel('local lift coefficient $c_l$')
      ....: plt.grid();

With creation of liftana automatically distributions for all control surfaces
are created.
The method :py:meth:`wingstructure.analysis.LiftAnalysis.calculate` combines
those to calculate local lift coefficient and angle of attack for a
requested lift. Its parameters allow also the enabling of spoilers and
deflection of other control surfaces. Have a look at
`this analysis <https://github.com/akafliegdarmstadt/wingstructure/blob/master/examples/Analysis_Example.ipynb>`_
and the documentation :py:meth:`wingstructure.analysis.calculate`.


Estimating wing section's mass
------------------------------

Calculating the mass of a wing section is not that closely connected to the
previous freatures. You start with the import of the
:py:mod:`wingstructure.structure`
module and loading of an airfoil *dat* file. 
For loading :py:func:`numpy.loadtxt` is used. 

.. ipython::
   
   In [12]: from wingstructure import structure

   In [13]: import numpy as np

   In [14]: coords = np.loadtxt('usage/FX 61-184.dat', skiprows=1, delimiter=',')

The loaded airfoil file can for examle be found 
`here <http://airfoiltools.com/airfoil/seligdatfile?airfoil=fx61184-il>`_

To analyse a section with this airfoil as outline an instance of
:py:obj:`wingstructure.structure.SectionBase` is created. 

.. ipython::

   In [14]: secbase = structure.SectionBase(coords)

Before filling this outline with structure we define a material:

.. ipython::

   In [15]: import collections as col

   In [16]: Material = col.namedtuple('Material', ['ρ'])

   In [17]: basematerial = Material(ρ=0.3e3)

Now the structure can be created. Definition begins from the outside
inwards. First feature is a constant layer (for example a fabric). We continue
with a spar definition and analyse the section's mass.

.. ipython::

   In [18]: layer = structure.Layer(secbase, basematerial, 3e-3)

   In [19]: spar = structure.ISpar(layer,
      ....:                        {'flange': basematerial,
      ....:                         'web': basematerial},
      ....:                          0.5,
      ....:                        0.300,
      ....:                        3e-2,
      ....:                        0.8,
      ....:                        1e-3)

   In [20]: massana = structure.MassAnalysis(spar)
      ....: massana.massproperties

For a more complex example have a look at the corresponding
`notebook <https://github.com/akafliegdarmstadt/wingstructure/blob/master/examples/Experimental_Mass_and_Structure.ipynb>`_
.