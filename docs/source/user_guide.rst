User guide
##########

Defining the geometry: :py:class:`~einsteinpy.metric` objects
*************************************************************

EinsteinPy provides a way to define the background geometry, on which the code would deal with the relativistic dynamics. This geometry has a central operating quantity, known as the Metric Tensor, that encapsulates all the geometrical and topological information about the 4D spacetime.

* EinsteinPy provides a :py:class:`~einsteinpy.metric.BaseMetric` class, that has various utility functions and a proper template, that can be used to define custom Metric classes. All pre-defined classes in :py:class:`~einsteinpy.metric` derive from this class.
* The central quantity required to simulate trajectory of a particle in a gravitational field are the metric derivatives, that can be succinctly written using Christoffel Symbols.
* EinsteinPy provides an easy to use interface to calculate these symbols.
* BaseMetric also provides support for ``f_vec`` and ``perturbation``, where ``f_vec`` corresponds to the RHS of the geodesic equation and ``perturbation`` is a linear Kerr-Schild Perturbation, that can be defined on the underlying metric.
* Note that, EinsteinPy does not perform physical checks on ``perturbation`` currently, and so, users should exercise caution while using it.

We provide an example below, showing how to calculate Time-like Geodesics in Schwarzschild spacetime.

Schwarzschild Metric
====================

EinsteinPy provides an intuitive interface for calculating time-like geodesics in Schwarzschild spacetime.

First of all, we import all the relevant modules and classes:

    .. code-block:: python

        import numpy as np

        from einsteinpy.coordinates.utils import four_position, stacked_vec
        from einsteinpy.geodesic import Geodesic
        from einsteinpy.metric import Schwarzschild


Defining initial parameters and our Metric Object
-------------------------------------------------

Now, we define the initial parameters, that specify the Schwarzschild metric and our test particle.

    .. code-block:: python

        M = 6e24  # Mass
        t = 0.  # Coordinate Time (has no effect in this case, as Schwarzschild metric is static)
        x_vec = np.array([130.0, np.pi / 2, -np.pi / 8])  # 3-Position of test particle
        v_vec = np.array([0.0, 0.0, 1900.0])  # 3-Velocity of test particle

        ms_cov = Schwarzschild(M=M) # Schwarzschild Metric Object
        x_4vec = four_position(t, x_vec) # Getting Position 4-Vector
        ms_cov_mat = ms_cov.metric_covariant(x_4vec) # Calculating Schwarzschild Metric at x_4vec
        init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True) # Contains 4-Pos and 4-Vel


Calculating Trajectory/Time-like Geodesic
-----------------------------------------
After creating the metric object and the initial vector, we can use :py:class:`~einsteinpy.geodesic.Geodesic` to create a Geodesic object, that automatically calculates the trajectory. 

    .. code-block:: python

        # Calculating Geodesic
        geod = Geodesic(metric=ms_cov, init_vec=init_vec, end_lambda=0.002, step_size=5e-8)
        # Getting a descriptive summary on geod
        print(geod)

    .. code-block:: python

        Geodesic Object:

        Metric = ((
        Name: (Schwarzschild Metric),            
        Coordinates: (S),            
        Mass: (6e+24),            
        Spin parameter: (0),            
        Charge: (0),            
        Schwarzschild Radius: (0.008911392322942397)
        )),            

        Initial Vector = ([ 0.00000000e+00  1.30000000e+02  1.57079633e+00 -3.92699082e-01
        1.00003462e+00  0.00000000e+00  0.00000000e+00  1.90000000e+03]),            

        Trajectory = ([[ 0.00000000e+00  1.20104339e+02 -4.97488462e+01 ...  9.45228078e+04
        2.28198245e+05  0.00000000e+00]
        [ 4.00013846e-08  1.20108103e+02 -4.97397110e+01 ...  9.36471118e+04
        2.28560931e+05 -5.80379473e-14]
        [ 4.40015231e-07  1.20143810e+02 -4.96475618e+01 ...  8.48885265e+04
        2.32184177e+05 -6.38424865e-13]
        ...
        [ 1.99928576e-03  1.29695466e+02 -6.52793459e-01 ...  1.20900076e+05
        2.46971585e+05 -1.86135457e-10]
        [ 1.99968577e-03  1.29741922e+02 -5.53995726e-01 ...  1.11380963e+05
        2.47015864e+05 -1.74024168e-10]
        [ 2.00008578e-03  1.29784572e+02 -4.55181739e-01 ...  1.01868292e+05
        2.47052855e+05 -1.61922169e-10]])


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

        from einsteinpy.plotting import GeodesicPlotter
        obj = GeodesicPlotter()
        obj.plot(geodesic)
        obj.show()


Utilities: :py:class:`~einsteinpy.utils`
****************************************

EinsteinPy provides a great set of utility functions which are frequently used in general and numerical relativity.

* Conversion of Coordinates (both position & velocity)

 * Cartesian/Spherical
 * Cartesian/Boyer-Lindquist

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
EinsteinPy also supports symbolic calculations in
:py:class:`~einsteinpy.symbolic`

    .. code-block:: python

        import sympy
        from einsteinpy.symbolic import SchwarzschildMetric, ChristoffelSymbols

        m = SchwarzschildMetric()
        ch = ChristoffelSymbols.from_metric(m)
        print(ch[1,2,:])

    .. code-block:: python

        [0, 0, -r*(-a/r + 1), 0]


    .. code-block:: python

        import sympy
        from einsteinpy.symbolic import SchwarzschildMetric, EinsteinTensor

        m = SchwarzschildMetric()
        G1 = EinsteinTensor.from_metric(m)
        print(G1.arr)

    .. code-block:: python

        [[a*c**2*(-a + r)/r**4 + a*c**2*(a - r)/r**4, 0, 0, 0], [0, a/(r**2*(a - r)) + a/(r**2*(-a + r)), 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]


Future Plans
============

* Support for null-geodesics in different geometries
* Ultimate goal is providing numerical solutions for Einstein's equations for arbitrarily complex matter distribution.
* Relativistic hydrodynamics
