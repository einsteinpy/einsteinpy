Coordinates module
==================

This module contains classes for representing positions and velocities in various
coordinate systems commonly used in general relativity and relativistic mechanics.
It provides support for:

- **Position representations**: Cartesian, Spherical, and Boyer-Lindquist coordinate systems
- **Velocity representations**: Differential classes for each coordinate system
- **Coordinate transformations**: Utilities for converting between different coordinate systems

The module uses SI units internally and integrates with Astropy's units system for
dimensional analysis and unit conversions.

Quick Start
-----------

Here's a simple example of creating and using coordinates:

.. code-block:: python

    from einsteinpy.coordinates import Cartesian
    import astropy.units as u

    # Create a position in Cartesian coordinates
    position = Cartesian(
        t=0 * u.s,
        x=1 * u.m,
        y=2 * u.m,
        z=3 * u.m
    )

    # Convert to Spherical coordinates
    spherical = position.to_spherical()

Coordinate Systems
------------------

The module supports the following coordinate systems:

1. **Cartesian (x, y, z)**: Standard 3D Cartesian coordinates
2. **Spherical (r, θ, φ)**: Radial distance and angles
3. **Boyer-Lindquist (r, θ, φ)**: Coordinates used in Kerr spacetime

Each coordinate system has a corresponding differential class for representing
velocities in that system.

.. toctree::
    :maxdepth: 2

    core
    differential
