Metric module
=============

This module contains classes for defining and working with metric tensors that describe
the geometry of spacetime in general relativity. A metric tensor defines how distances
and times are measured in curved spacetime, and is essential for computing geodesics
and understanding gravitational effects.

The module provides:

- **BaseMetric**: A base class for defining custom metric tensors
- **Schwarzschild**: The Schwarzschild metric for non-rotating, uncharged black holes
- **Kerr**: The Kerr metric for rotating black holes
- **KerrNewman**: The Kerr-Newman metric for rotating, charged black holes

Each metric class allows you to compute covariant and contravariant forms of the metric
tensor, Christoffel symbols, and other geometric quantities needed for relativistic
calculations.

Quick Start
-----------

.. code-block:: python

    from einsteinpy.metric import Schwarzschild
    import astropy.units as u

    # Create a Schwarzschild metric at a given position
    metric = Schwarzschild(M=1 * u.kg)

    # Get metric tensor components
    g = metric()
    g_lower = metric.lower()
    g_upper = metric.upper()

.. toctree::
    :maxdepth: 2

    base_metric
    schwarzschild
    kerr
    kerrnewman
