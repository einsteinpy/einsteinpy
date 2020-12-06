Jupyter notebooks
=================
   
.. toctree::
   :maxdepth: 1
   :caption: hypersurface
   
   This module provides basic computational functions which are essential for modelling the space-like hypersurface for any space-time geometry. Currently, the    module has an implementation for Schwarzschild geometry with conversion functions from Schwarzschild coordinates to 3-D spherical coordinates and some plotting utilities.
   
   /examples/Plotting spacial hypersurface embedding for schwarzschild spacetime.ipynb

.. toctree::
   :maxdepth: 1
   :caption: rays
   
   The rays module of EinsteinPy tries to focus on how light behaves around heavy mass object. The major functionality this moduleprovides is the utility functions for carrying out shadow related computations for different black hole spacetimes.
   
   /examples/Shadow cast by an thin emission disk around a black hole.ipynb

.. toctree::
   :maxdepth: 1
   :caption: metric, bodies, geodesic, coordinates
   
   * Metric - This module captures all the geometric and causal structure of specific space-times, and can calculate trajectories in a given space-time. The module derives its name from Metric Tensor, a tensorial quantity used to describe differential lengths, especially in curvedspace-time.
   * Body - The Body data type helps define an attractor i.e. the central black-hole or the heavy mass object and later the objects under itsinfluence. The bodies can have their characteristics defined then and there thus helping the user create the system easily.
   * Geodesic - A geodesic is the path that an-accelerating particle would follow. In the plane, the geodesics are straight lines. On the sphere, thegeodesics are great circles. The Geodesic data type helps define the geodesic of a body in the presence of an attractor or more aptly,according to the system user created. TheseGeodesicobjects can be defined for any configuration and metric.
   * Co-Ordinates - The coordinates subpackage provides support for representing and transforming coordinates. EinsteinPy provides support in Carte-sian, Spherical and Boyer-Lindquist Coordinate systems. All of these are inter-convertible and can be used anytime.
   
   /examples/Animations in EinsteinPy.ipynb
   /examples/Using Geodesics (Back-ends & Plotting).ipynb
   /examples/Visualizing Event Horizon and Ergosphere (Singularities) of Kerr Metric or Black Hole.ipynb
   /examples/Visualizing Frame Dragging in Kerr Spacetime using EinsteinPy!.ipynb
   /examples/Visualizing Precession in Schwarzschild Spacetime.ipynb
   
.. toctree::
   :maxdepth: 1
   :caption: symbolic
   
   /examples/Einstein Tensor symbolic calculation.ipynb
   /examples/Lambdify symbolic calculation.ipynb
   /examples/Playing with Contravariant and Covariant Indices in Tensors(Symbolic).ipynb
   /examples/Predefined Metrics in Symbolic Module.ipynb
   /examples/Ricci Tensor and Scalar Curvature symbolic calculation.ipynb
   /examples/Symbolically Understanding Christoffel Symbol and Riemann Curvature Tensor using EinsteinPy.ipynb
   /examples/Weyl Tensor symbolic calculation.ipynb
   
