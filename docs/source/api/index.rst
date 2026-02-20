EinsteinPy API
==============

Welcome to the API documentation for EinsteinPy. This reference guide provides detailed 
information about the classes, methods, and modules available in EinsteinPy. Navigate 
through the module sections below to learn about the API.

If you find anything missing or have suggestions for improvement, please open an 
`issue in the repo`_.

.. _`issue in the repo` : https://github.com/einsteinpy/einsteinpy/issues

Module Overview
---------------

The following graph shows the main modules available in EinsteinPy and their relationships:

.. graphviz::

    digraph {
        "einsteinpy" -> "integrators", "metric", "utils", "plotting", "constant", "units", "coordinates", "geodesic", "bodies", "hypersurface", "symbolic", "rays"

    }

.. toctree::
    :maxdepth: 2

    integrators/integrators_index
    metric/metric_index
    symbolic/symbolic_index
    hypersurface/hypersurface_index
    rays/rays_index
    utils/utils_index
    plotting/plotting_index
    coordinates/coordinates_index
    constant
    units
    bodies
    geodesic/geodesic_index
    examples
