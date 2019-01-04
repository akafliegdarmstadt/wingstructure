===============
Getting Started
===============

This document is intended as starting point for wingstructure usage.
The wingstructure package has three principal functions:

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

   @savefig plot_simple.png width=4in
   In [4]: wing.plot()

Besides plotting the wing object has some useful properties. For example
calculating the mean aerodynamic chord:

.. ipython::

   In [5]: wing.mac

Further information can be found in :py:class:`wingstructure.data.Wing`.

Calculating wing's lift distribution
------------------------------------



Estimating wing section's mass
------------------------------

