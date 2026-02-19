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

.. toctree::
    :maxdepth: 2

    core
    differential
