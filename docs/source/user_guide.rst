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

