User Guide
==========

Plotting and Analysing Schwarzschild Geodesics: :py:class:`~einsteinpy.metric.schwarzschild.Schwarzchild`
---------------------------------------------------------------------------------------------------------

The core of einsteinpy as of now are the :py:class:`~einsteinpy.metric.schwarzschild.Schwarzchild` objects
inside the :py:mod:`einsteinpy.metric` module. They store all the required information to define the geometry
according to the matrix:

* The mass of the body acting as the central body of the geometry, for example the Earth.
* The position and velocity vectors for any particle.
* The time at which the geomtery is defined.

First of all, we have to import the relevant modules and classes:

.. code-block:: python

    from astropy import units as u

    from einsteinpy.plotting import StaticGeodesicPlotter
    from einsteinpy.metric import Schwarzchild

From Spherical Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    pos_vec = [147.09e6 * u.km, np.pi / 2 * u.rad, np.pi * u.rad]
    vel_vec = [0 * u.km / u.s, 0 * u.rad / u.s, omega]
    m = 1.989e30 * u.kg

    ss = StaticGeodesicPlotter(m)
    ss.plot(pos_vec, vel_vec)

One can even plot the perihelion of mercury using this.
