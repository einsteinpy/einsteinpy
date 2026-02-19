Integrators
===========

This module provides numerical integration methods for solving Ordinary Differential
Equations (ODEs) in the context of numerical relativity and relativistic mechanics.
These integrators are used to solve the geodesic equation and other dynamical equations
in general relativity.

Available Integrators
--------------------

The module includes several integration methods:

1. **RK4naive**: A simple 4th-order Runge-Kutta method for educational purposes
2. **RK45**: A 5th-order Runge-Kutta method with adaptive step size control
3. **Cash-Karp RK45**: An alternative adaptive method with error estimation
4. **FANTASY**: Uses Forward Mode Automatic Differentiation for metric derivatives

The FANTASY integrator supports orders 2, 4, 6, and 8 and is particularly useful
when working with custom metrics where analytical derivatives are difficult to compute.

Quick Start
-----------

.. code-block:: python

    from einsteinpy.integrators import RK45
    import astropy.units as u

    # Define your ODE function
    def rhs(t, state):
        # Your differential equations here
        return state

    # Create and use integrator
    integrator = RK45(
        func=rhs,
        t0=0 * u.s,
        y0=[1, 0, 0, 0],
        t_bound=100 * u.s
    )

.. toctree::
    :maxdepth: 2

    runge_kutta
    fantasy
