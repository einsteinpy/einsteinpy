Integrators
===========

This module provides numerical integration methods for solving Ordinary Differential 
Equations (ODEs) in the context of numerical relativity and relativistic mechanics. 
These integrators are used to solve the geodesic equation and other dynamical equations 
in general relativity.

Available integrators include:

- **RK4naive**: A simple 4th-order Runge-Kutta method for educational purposes
- **RK45**: A more sophisticated 5th-order Runge-Kutta method with adaptive step size control
- **GeodesicIntegrator**: Specialized integrator for computing geodesic trajectories in curved spacetime

These integrators work with functions that define the right-hand side of ODEs and can 
handle both scalar and vector-valued systems.

.. toctree::
    :maxdepth: 2

    runge_kutta
    fantasy
