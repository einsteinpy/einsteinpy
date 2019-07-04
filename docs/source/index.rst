EinsteinPy - Making Einstein possible in Python
===============================================

.. image:: https://einsteinpy.org/img/wordmark.png
   :target: https://einsteinpy.org/
   :alt: EinsteinPy logo
   :width: 675px
   :align: center


**EinsteinPy** is an open source pure Python package dedicated to problems arising
in General Relativity and gravitational physics, such as goedesics plotting for
Schwarzschild, Kerr and Kerr Newman space-time model, calculation of Schwarzschild
radius, calculation of Event Horizon and Ergosphere for Kerr space-ime. Symbolic
Manipulations of various tensors like Metric, Riemann, Ricci and Christoffel Symbols
is also possible using the library. EinsteinPy also features Hypersurface Embedding of
Schwarzschild space-time, which will soon lead to modelling of Gravitational Lensing!
It is released under the MIT license.


View `source code`_ of EinsteinPy!

.. _`source code`: https://github.com/einsteinpy/einsteinpy


Key features of EinsteinPy are:

* Geometry analysis and trajectory calculation in vaccum solutions of Einstein's field equations
 
 * Schwarzschild space-time
 * Kerr space-time
 * Kerr-Newman space-time

* Various utilities related to above geometry models

 * Schwarzschild Radius
 * Event horizon and ergosphere for Kerr black hole
 * Maxwell Tensor and electromagnetic potential in Kerr-Newman space-time
 * And much more!

* Symbolic Calculation of various quantities in GR

 * Christoffel Symbols
 * Riemann Curvature Tensor
 * Simplification of symbolic expressions

* Static Geodesic Plotting
* Coordinate conversion with unit handling

 * Spherical/Cartesian Coordinates
 * Boyer-Lindquist/Cartesian Coordinates


And more to come!


Einsteinpy is developed by an open community. Release
announcements and general discussion take place on our `mailing list`_
and `chat`_.

.. _`mailing list`: https://groups.io/g/einsteinpy-dev
.. _`chat`: https://riot.im/app/#/room/#einsteinpy:matrix.org

.. include:: form.rst


The `source code`_, `issue tracker`_ and `wiki`_ are hosted on GitHub, and all
contributions and feedback are more than welcome. You can test EinsteinPy in your
browser using binder, a cloud Jupyter notebook server:

.. image:: https://img.shields.io/badge/launch-binder-e66581.svg?style=flat-square
   :target: https://beta.mybinder.org/v2/gh/einsteinpy/einsteinpy/master?filepath=index.ipynb

.. _`source code`: https://github.com/einsteinpy/einsteinpy
.. _`issue tracker`: https://github.com/einsteinpy/einsteinpy/issues
.. _`wiki`: https://github.com/einsteinpy/einsteinpy/wiki/



EinsteinPy works on recent versions of Python and is released under
the MIT license, hence allowing commercial use of the library.

.. code-block:: python

    from einsteinpy.plotting import StaticGeodesicPlotter
    a = StaticGeodesicPlotter(mass)
    a.plot(r,v)

Contents
--------

.. toctree::
    :maxdepth: 2

    getting_started
    user_guide
    metric
    jupyter
    changelog
    dev_guide
    api/index
