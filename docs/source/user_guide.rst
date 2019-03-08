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

        from einsteinpy.metric import Schwarzschild

From position and velocity in Spherical Coordinates
---------------------------------------------------

There are several methods available to create :py:class:`~einsteinpy.metric.schwarzschild.Schwarzschild` objects. For example, if we have the position and velocity vectors we can use :py:meth:`~einsteinpy.metric.schwarzschild.Schwarzschild.from_spherical`:

    .. code-block:: python

        M = 5.972e24 * u.kg
        pos_vec = [306.0 * u.m, np.pi/2 * u.rad, -30 * u.deg]
        vel_vec = [0 * u.km/u.s, 0 * u.rad/u.s, 950.69 * u.rad/u.s]
        obj = Schwarzschild.from_spherical(pos_vec, vel_vec, 0 * u.s, M)

From position and velocity in Cartesian Coordinates
---------------------------------------------------
For initializing with Cartesian Coordinates, we can use :py:class:`~einsteinpy.metric.schwarzschild.Schwarzschild.from_cartesian`:

    .. code-block:: python

        pos_vec = [265.003774 * u.m, -153.000000e * u.m,  0 * u.m]
        vel_vec = [145455.57 * u.m/u.s, 251936.43748389 * u.m/u.s, 0 * u.km/u.s]
        obj = Schwarzschild.from_cartesian(pos_vec, vel_vec, 0*u.s, M)

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

        from einsteinpy import utils

        pos_vec = np.array([200, -100, 20.5])
        vel_vec = np.array([-12, 14, 0.5])
        a = 0.5
        bl_pos_vec = utils.CartesianToBL_pos(pos_vec, a)
        bl_vel_vec = utils.CartesianToBL_vel(pos_vec, vel_vec, a)
        cs_pos_vec = utils.BLToCartesian_pos(bl_pos_vec, a)
        cs_vel_vec = utils.BLToCartesian_vel(bl_pos_vec, bl_vel_vec, a)
        print(pos_vec)
        print(bl_pos_vec)

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

* Support for more geometries, like Kerr, Kerr-Newman
* Support for null-geodesics in different geometries
* Ultimate goal is providing numerical solutions for Einstein's equations for arbitarily complex matter distribution.
* Relativistic hydrodynamics