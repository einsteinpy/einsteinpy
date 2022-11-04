EinsteinPy API
==============

Welcome to the API documentation of EinsteinPy. Please navigate through the given modules to get to
know the API of the classes and methods. If you find anything missing, please open an `issue in the repo`_ .

.. _`issue in the repo` : https://github.com/einsteinpy/einsteinpy/issues

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
