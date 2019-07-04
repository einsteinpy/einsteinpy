What's new
==========

einsteinpy 0.2.0 - 2019-07-15
-----------------------------

This release will bring a lots and lots of features with a more intuitive API, a better
plotting support, more metrics introduced and what not. This is a major release which will
short term support version and will be supported only till June 2020. This major release
comes before Python in Astronomy 2019 workshop.

Part of this release is sponsored by European Space Agency, through Summer of Code in Space
(SOCIS) 2019 program.

Features
........

* Kerr Metric
* Kerr-Newman Metric
* Coordinates Module with Boyer Lindquist Coordinates and transformation
* Bodies Module
* Defining Geodesics with ease!
* Animated plots
* Intuitive API for plotting
* Schwarzschild Hypersurface Embedding
* Interactive Plotting
* Environment-aware plotting and exceptional support for iPython Notebooks!
* Support for Tensor Algebra in General Relativity
* Symbolic Manipulation of Metric Tensor, Riemann Tensor and Ricci Tensor
* Support for Index Raising and Lowering in Tensors
* Numerical Calculation and Symbolic Manupulation of Christoffel Symbols
* Calculations of Event Horizon and Ergosphere of Kerr Black holes!


Contributors
............

This is the complete list of the people that contributed to this release, with a + sign indicating first contribution.

* Shreyas Bapat
* Ritwik Saha
* Bhavya Bhatt
* Sofía Ortín Vela+
* Raphael Reyna+
* Prithvi Manoj Krishna+
* Manvi Gupta+
* Divya Gupta+
* Yash Sharma+
* Shilpi Jain+
* Rishi Sharma+
* Varun Singh+
* Alpesh Jamgade+
* Saurabh Bansal+
* Tanmay Rustagi+
* Abhijeet Manhas+
* Ankit Khandelwal+
* Tushar Tyagi+
* Hrishikesh Sarode
* Naman Tayal+
* Ratin Kumar+
* Govind Dixit+
* Jialin Ma+

Bugs Fixed
----------

* `Issue #115`_: Coordinate Conversion had naming issues that made them confusing!
* `Issue #185`_: isort had conflicts with Black
* `Issue #210`_: Same notebook had two different listings in Example Gallery
* `Issue #264`_: Removing all relative imports
* `Issue #265`_: New modules were lacking API Docs
* `Issue #266`_: The logo on documentation was not rendering
* `Issue #267`_: Docs were not present for Ricci Tensor and Vacuum Metrics
* `Issue #277`_: Coordinate Conversion in plotting module was handled incorrectly


.. _`Issue #115`: https://github.com/einsteinpy/einsteinpy/issues/115
.. _`Issue #185`: https://github.com/einsteinpy/einsteinpy/issues/185
.. _`Issue #210`: https://github.com/einsteinpy/einsteinpy/issues/210
.. _`Issue #264`: https://github.com/einsteinpy/einsteinpy/issues/264
.. _`Issue #265`: https://github.com/einsteinpy/einsteinpy/issues/265
.. _`Issue #266`: https://github.com/einsteinpy/einsteinpy/issues/266
.. _`Issue #267`: https://github.com/einsteinpy/einsteinpy/issues/267
.. _`Issue #277`: https://github.com/einsteinpy/einsteinpy/issues/277

Backwards incompatible changes
..............................

* The old :code:`StaticGeodesicPlotter` has been renamed to
  :py:class:`einsteinpy.plotting.senile.StaticGeodesicPlotter`, please adjust
  your imports accordingly
* The old :code:`ScatterGeodesicPlotter` has been renamed to
  :py:class:`einsteinpy.plotting.senile.ScatterGeodesicPlotter`, please adjust
  your imports accordingly.
* :py:class:`einsteinpy.metric.Schwarzschild`,
  :py:class:`einsteinpy.metric.Kerr`, and
  :py:class:`einsteinpy.metric.KerrNewman` now have different signatures for
  class methods, and they now explicitly support :py:mod:`einsteinpy.coordinates`
  coordinate objects. Check out the notebooks and their respective documentation.
* The old `coordinates` conversion in :py:mod:`einsteinpy.utils` has been deprecated.
* The old `symbolic` module in :py:mod:`einsteinpy.utils` has been moved to
  :py:mod:`einsteinpy.symbolic`.

einsteinpy 0.1.0 - 2019-03-08
-----------------------------

This is a major first release for world's first actively maintained python library
for General Relativity and Numerical methods. This major release just comes before
the Annual AstroMeet of IIT Mandi, AstraX. This will be a short term support version
and will be supported only until late 2019.

Features
........

* Schwarzschild Geometry Analysis and trajectory calculation
* Symbolic Calculation of various tensors in GR

 * Christoffel Symbols
 * Riemann Curvature Tensor

* Static Geodesic Plotting
* Velocity of Coordinate time w.r.t proper time
* Easy Calculation of Schwarzschild Radius
* Coordinate conversion with unit handling

 * Spherical/Cartesian Coordinates
 * Boyer-Lindquist/Cartesian Coordinates


Contributors
............

This is the complete list of the people that contributed to this release, with a + sign indicating first contribution.

* Shreyas Bapat+
* Ritwik Saha+
* Bhavya Bhatt+
* Priyanshu Khandelwal+
* Gaurav Kumar+
* Hrishikesh Sarode+
* Sashank Mishra+
