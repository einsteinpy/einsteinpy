Geodesic
========

This module provides classes for computing geodesic trajectories in curved spacetime.
Geodesics are the paths that particles follow when influenced only by gravity, representing
the generalization of straight lines to curved spacetime.

The module defines:

- **Geodesic**: Base class for all geodesic calculations
- **TimeLikeGeodesic**: Geodesics for massive particles (follows proper time)
- **NullGeodesic**: Geodesics for massless particles like photons (follows affine parameter)

These classes are essential for simulating particle trajectories around black holes
and studying gravitational lensing, orbital dynamics, and other relativistic phenomena.

Quick Start
-----------

.. code-block:: python

    from einsteinpy.metric import Schwarzschild
    from einsteinpy.geodesic import TimeLikeGeodesic
    import astropy.units as u

    # Create a metric
    metric = Schwarzschild(M=1 * u.kg)

    # Define initial position and four-velocity
    # Using Boyer-Lindquist coordinates
    position = [0, 10, 0, 0]  # r, theta, phi in SI units
    four_velocity = [0, 0, 1, 0]  # dr/dt, dtheta/dt, dphi/dt, etc.

    # Create and compute geodesic
    geodesic = TimeLikeGeodesic(metric, position, four_velocity)
    trajectory = geodesic.calculate_trajectory(
        steps=1000,
        step_size=0.1 * u.s
    )

.. automodule:: einsteinpy.geodesic.geodesic
    :members:
    :show-inheritance:
