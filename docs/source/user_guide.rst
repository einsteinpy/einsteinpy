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
        obj = Schwarzschild.from_coords(sph_coord, M , 0* u.s)

From position and velocity in Cartesian Coordinates
---------------------------------------------------
For initializing with Cartesian Coordinates, we can use :py:class:`~einsteinpy.metric.schwarzschild.Schwarzschild.from_cartesian`:

    .. code-block:: python

        cartsn_coord = CartesianDifferential(.265003774 * u.km, -153.000000e-03 * u.km,  0 * u.km,
                          145.45557 * u.km/u.s, 251.93643748389 * u.km/u.s, 0 * u.km/u.s)
        obj = Schwarzschild.from_coords(cartsn_coord, M , 0* u.s)

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


Bodies Module: :py:class:`~einsteinpy.bodies`
*********************************************

EinsteinPy has a module to define the attractor and revolving bodies, using which plotting and geodesic calculation 
becomes much easier.

Importing all the relevant modules and classes :

    .. code-block:: python

        import numpy as np
        from astropy import units as u
        from einsteinpy.coordinates import BoyerLindquistDifferential
        from einsteinpy.metric import Kerr
        from einsteinpy.bodies import Body
        from einsteinpy.geodesic import Geodesic


Defining various astronomical bodies :

    .. code-block:: python

        spin_factor = 0.3 * u.m
        Attractor = Body(name="BH", mass = 1.989e30 * u.kg, a = spin_factor)
        BL_obj = BoyerLindquistDifferential(50e5 * u.km, np.pi / 2 * u.rad, np.pi * u.rad,
                                            0 * u.km / u.s, 0 * u.rad / u.s, 0 * u.rad / u.s,
                                            spin_factor)
        Particle = Body(differential = BL_obj, parent = Attractor)
        geodesic = Geodesic(body = Particle, end_lambda = ((1 * u.year).to(u.s)).value / 930,
                            step_size = ((0.02 * u.min).to(u.s)).value,
                            metric=Kerr)
        geodesic.trajectory  # get the values of the trajectory


Plotting the trajectory :

    .. code-block:: python

        from einsteinpy.plotting import ScatterGeodesicPlotter
        obj = ScatterGeodesicPlotter()
        obj.plot(geodesic)
        obj.show()


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

* :py:class:`~einsteinpy.coordinates.BoyerLindquistDifferential.to_cartesian`
* :py:class:`~einsteinpy.coordinates.CartesianDifferential.to_bl`

    .. code-block:: python

        import numpy as np
        from astropy import units as u
        from einsteinpy.coordinates import BoyerLindquistDifferential, CartesianDifferential, Cartesian, BoyerLindquist

        a = 0.5 * u.km

        pos_vec = Cartesian(.265003774 * u.km, -153.000000e-03 * u.km,  0 * u.km)

        bl_pos = pos_vec.to_bl(a)
        print(bl_pos)

        cartsn_pos = bl_pos.to_cartesian(a)
        print(cartsn_pos)

        pos_vel_coord = CartesianDifferential(.265003774 * u.km, -153.000000e-03 * u.km,  0 * u.km,
                                  145.45557 * u.km/u.s, 251.93643748389 * u.km/u.s, 0 * u.km/u.s)

        bl_coord = pos_vel_coord.bl_differential(a)
        bl_coord = bl_coord.si_values()
        bl_vel = bl_coord[3:]
        print(bl_vel)

        cartsn_coord = bl_coord.cartesian_differential(a)
        cartsn_coord = cartsn_coord.si_values()
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
        from einsteinpy.symbolic import SchwarzschildMetric, ChristoffelSymbols

        m = SchwarzschildMetric()
        ch = ChristoffelSymbols.from_metric(m)
        print(ch[1,2,:])

    .. code-block:: python

        [0, 0, -r*(-a/r + 1), 0]


Future Plans
============

* Support for null-geodesics in different geometries
* Ultimate goal is providing numerical solutions for Einstein's equations for arbitarily complex matter distribution.
* Relativistic hydrodynamics
