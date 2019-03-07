EinsteinPy - Making Einstein possible in Python
===============================================

.. image:: http://einsteinpy.github.io/img/logo.png
   :target: http://einsteinpy.github.io/
   :alt: EinsteinPy logo
   :width: 675px
   :align: center


**EinsteinPy** is an open source pure Python package dedicated to problems arising
in General Relativity and relativistic physics, such as goedesics plotting for
schwartzschild space-time model, calculation of schwartzschild radius for any mass given.
Features like visualisation of geodesics of curved black holes and 3D visualisations
are some of the features which are planned. It is released under the MIT license.


View `source code`_ of EinsteinPy!

.. _`source code`: https://github.com/einsteinpy/einsteinpy


Key features of EinsteinPy are:

* Schwarzschild Geometry Analysis and trajectory calculation
* Symbolic Calculation of various tensors in GR

 * Christoffel Symbols
 * Riemann Curvature Tensor

* Static Geodesic Plotting
* Velocity of Coordinate time w.r.t proper time
* Easy Calculation of Schwarzschild Radius
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
    jupyter
    changelog
    api/index
