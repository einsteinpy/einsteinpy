User guide
##########

Defining the geometry: :py:class:`~einsteinpy.metric` objects
*************************************************************

EinsteinPy provides a way to define the background geometry on which the code would deal with the dynamics. These geometry has a central operating quantity known as metric tensor  and encapsulate all the geometrical and topological information about the 4d spacetime in them.

* The central quantity required to simulate trajectory of a particle in a gravitational field is christoffel symbols.
* EinsteinPy provides an easy to use interface to calculate these symbols.

Schwarzschild metric
====================

EinsteinPy provides an easy interface for calculating time-like geodesics in Schwarzschild Geometry.

First of all, we import all the relevant modules and classes :

    .. code-block:: python

        import numpy as np
        from astropy import units as u
        from einsteinpy.coordinates import SphericalDifferential, CartesianDifferential
        from einsteinpy.metric import Schwarzschild

From position and velocity in Spherical Coordinates
---------------------------------------------------

There are several methods available to create :py:class:`~einsteinpy.metric.schwarzschild.Schwarzschild` objects. For example, if we have the position and velocity vectors we can use :py:meth:`~einsteinpy.metric.schwarzschild.Schwarzschild.from_spherical`:

    .. code-block:: python

        M = 5.972e24 * u.kg
        sph_coord = SphericalDifferential(306.0 * u.m, np.pi/2 * u.rad, -np.pi/6*u.rad,
                                  0*u.m/u.s, 0*u.rad/u.s, 1900*u.rad/u.s)
        obj = Schwarzschild.from_spherical(sph_coord, M , 0* u.s)

From position and velocity in Cartesian Coordinates
---------------------------------------------------
For initializing with Cartesian Coordinates, we can use :py:class:`~einsteinpy.metric.schwarzschild.Schwarzschild.from_cartesian`:

    .. code-block:: python

        cartsn_coord = CartesianDifferential(.265003774 * u.km, -153.000000e-03 * u.km,  0 * u.km,
                          145.45557 * u.km/u.s, 251.93643748389 * u.km/u.s, 0 * u.km/u.s)
        obj = Schwarzschild.from_cartesian(cartsn_coord, M , 0* u.s)

Calculating Trajectory/Time-like Geodesics
------------------------------------------
After creating the object we can call :py:class:`~einsteinpy.metric.schwarzschild.Schwarzschild.calculate_trajectory`

    .. code-block:: python

        end_tau = 0.01 # approximately equal to coordinate time
        stepsize = 0.3e-6
        ans = obj.calculate_trajectory(end_lambda=end_tau, OdeMethodKwargs={"stepsize":stepsize})
        print(ans)

    .. code-block:: python

        (array([0.00000000e+00, 2.40000000e-07, 2.64000000e-06, ...,
            9.99367909e-03, 9.99607909e-03, 9.99847909e-03]), array([[ 0.00000000e+00,  3.06000000e+02,  1.57079633e+00, ...,
                0.00000000e+00,  0.00000000e+00,  9.50690000e+02],
            [ 2.39996635e-07,  3.05999885e+02,  1.57079633e+00, ...,
                -9.55164950e+02,  1.32822112e-17,  9.50690712e+02],
            [ 2.63996298e-06,  3.05986131e+02,  1.57079633e+00, ...,
                -1.05071184e+04,  1.46121838e-16,  9.50776184e+02],
            ...,
            [ 9.99381048e-03,  3.05156192e+02,  1.57079633e+00, ...,
                8.30642520e+04, -1.99760372e-12,  9.55955926e+02],
            [ 9.99621044e-03,  3.05344028e+02,  1.57079633e+00, ...,
                7.34673728e+04, -2.01494258e-12,  9.54780155e+02],
            [ 9.99861041e-03,  3.05508844e+02,  1.57079633e+00, ...,
                6.38811856e+04, -2.03252073e-12,  9.53750261e+02]]))

Return value can be obtained in Cartesian Coordinates by :

    .. code-block:: python

        ans = obj.calculate_trajectory(end_lambda=end_tau, OdeMethodKwargs={"stepsize":stepsize}, return_cartesian=True)


Utilities: :py:class:`~einsteinpy.utils`
****************************************

EinsteinPy provides a great set of utility functions which are frequently used in general and numerical relativity.

* Conversion of Coordinates (both position & velocity)

 * Cartesian/Spherical
 * Cartesian/Boyer-Lindquist

* Symbolic calculation of relevant terms in GR

 * Christoffel Symbols
 * Riemann Curvature Tensor

* Calculation of Schwarzschild Geometry related quantities

 * Schwarzschild Radius
 * Rate of change of coordinate time w.r.t. proper time

Coordinate Conversion
=====================

In a short example, we would see coordinate conversion between Cartesian and Boyer-Lindquist Coordinates.

Using the functions:

* :py:class:`~einsteinpy.utils.bl_coord_transforms.CartesianToBL_pos`
* :py:class:`~einsteinpy.utils.bl_coord_transforms.CartesianToBL_vel`
* :py:class:`~einsteinpy.utils.bl_coord_transforms.BLToCartesian_pos`
* :py:class:`~einsteinpy.utils.bl_coord_transforms.BLToCartesian_vel`

    .. code-block:: python

        import numpy as np
        from astropy import units as u
f       from einsteinpy.coordinates import BoyerLindquistDifferential, CartesianDifferential, Cartesian, BoyerLindquist

        a = 0.5 * u.km

        pos_vec = Cartesian(.265003774 * u.km, -153.000000e-03 * u.km,  0 * u.km)

        bl_pos = pos_vec.to_bl(a)
        print(bl_pos)

        cartsn_pos = bl_pos.to_cartesian(a)
        print(cartsn_pos)

        pos_vel_coord = CartesianDifferential(.265003774 * u.km, -153.000000e-03 * u.km,  0 * u.km,
                                  145.45557 * u.km/u.s, 251.93643748389 * u.km/u.s, 0 * u.km/u.s)

        bl_coord = pos_vel_coord.bl_differential(a)
        bl_coord = np.array(bl_coord)
        bl_vel = bl_coord[3:]
        print(bl_vel)

        cartsn_coord = bl_coord.cartesian_differential(a)
        cartsn_coord = np.array(cartsn_coord)
        cartsn_vel = cartsn_coord[3:]
        print(cartsn_vel)


    .. code-block:: python

        [ 200.  -100.    20.5]
        [224.54398697   1.47937288  -0.46364761]

Symbolic Calculations
=====================
EinsteinPy also supports smbolic calculations in :py:class:`~einsteinpy.utils.christoffel`

    .. code-block:: python

        import sympy
        from einsteinpy.utils import christoffel

        syms = sympy.symbols('t r theta phi')
        kch = christoffel.christkerr_christoffels()
        print(sympy.simplify(kch[0][0][1]))

    .. code-block:: python

        R*(-a**4*cos(theta)**2 - a**2*r**2*cos(theta)**2 + a**2*r**2 + r**4)/(2*(a**2*cos(theta)**2 + r**2)**2*(-R*r + a**2 + r**2))


Future Plans
============
* Support for null-geodesics in different geometries
* Ultimate goal is providing numerical solutions for Einstein's equations for arbitarily complex matter distribution.
* Relativistic hydrodynamics
