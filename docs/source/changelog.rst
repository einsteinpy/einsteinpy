What's New
==========

einsteinpy 0.4.0 - 2021-05-06
-----------------------------

This major release brings a lot of improvements to the numerical module of the project. COVID19 is proving
to be very difficult for India and we are trying to cope up. Please forgive us for any issues you faced with
the v0.3.1 and the documentation.

Changes
.......

* [#506]: Tests moved outside of the package.
* [#510]: Added a utility function for Outer Product and Contraction of Tensors in the symbolic module.
* [#512]: Refactored metric, utils and geodesic modules, added metric.BaseMetric class
* [#512]: Fixed #113, Fixed cyclomatic complexity issue in metric.Schwarzschild and metric.Kerr classes
* [#512]: Fixed #141, Refactored utils and merged most utilities into metric itself
* [#512]: Fixed #410, Improved __str__ and __repr__ messages for Geodesic objects
* [#512]: Fixed #507, Fixed a mathematical inaccuracy in metric.Schwarzschild class
* [#512]: Fixed #508, Removed a stray scaling factor in metric.KerrNewman class
* [#512]: Fixed #514, Replaced Spin Parameter with Length Parameter throughout metric module
* [#512]: Fixed #515, Renamed "maxwell" to "em"
* [#521]: Refactored coordinates, geodesic and metric modules, added support for 4-Vectors
* [#521]: Fixed #517, Removed Spin Parameter from bodies
* [#521]: Fixed #524, Fixed breakage, caused due to isort changes
* [#521]: Partially fixed #523, Fixed Schwarzschild and Kerr
* [#527]: Added support for Null Geodesics in Kerr & Schwarzschild Spacetimes
* [#527]: Added new features to plotting.geodesic
* [#527]: Dropped support for Python 3.6
* [#547]: Fixed #516, Added __all__ across modules
* [#547]: Fixed #552, Renamed missing attributes in einsteinpy.plotting.geodesic.static
* [#547]: Increased test coverage for einsteinpy.ijit
* [#551]: Fixed #526, Exceptions module added, CoordinateError class implemented
* [#565]: Fixed #549, Updated einsteinpy.symbolic Jupyter Notebooks
* [#548]: Fixed #36, Added support for animated Geodesic plots
* [#548]: Fixed #40, Added support for Order 4, 6 & 8 Geodesic Integrators
* [#548]: Fixed #105, Added support for simulating Null Geodesics in Kerr & Schwarzschild geometries
* [#548]: Fixed #122, Schwarzschild & Kerr algorithms validated
* [#548]: Fixed #367, Scaling issues fixed for Frame Dragging
* [#548]: Fixed #535, Moved to a pure python geodesic solver, Julia dependency removed
* [#548]: Minor edits to documentation for geodesic and plotting.geodesic modules
* [#571]: Fixed #570, Updated Master to Main
* [#573]: Fixed bug in Riemann Tensor calculation from Christoffels

Contributors
............

This is the complete list of the people that contributed to this release, with a + sign indicating first contribution.

* Shreyas Bapat
* Jyotirmaya Shivottam
* Bibek Gautam
* Qbiwan+ (GitHub Username)
* Aditya Prashant Dalvi+
* Aditya Prakash+
* aweinr4+ (GitHub Username)

einsteinpy 0.3.1 - 2021-01-16
-----------------------------

This release is a minor patch release for fixing a minor Debian issue.

Contributors
............

This is the complete list of the people that contributed to this release, with a + sign indicating first contribution.

* Shreyas Bapat


einsteinpy 0.3.0 - 2020-05-05
-----------------------------

This major release would bring some very important improvements. This release fixes a very crucial
bug with sympy. Fixes coordinate conversions so they don't fail on edge cases anymore.

EinsteinPy now uses GitHub Actions for macOS builds. Big changes to the plotting module.

The release comes for the paper of EinsteinPy. The release marks the beginning of Google Summer of Code 2020.
The release also brings a new rays module, which will form the base for null geodesics in future releases.

Features
........

* Loads of Predefined Metrics
* Sympy version incompatibilities handled
* Numba as a default installation
* Lorentz Transform for Einstein Tensor
* Lorentz Transform to Tensor Class
* Hypersurface Plotting API similar to the common plotting API
* Find Function in Predefined Metrics
* Increased Code Coverage
* New rays module
* Plotting Black Hole Shadows
* Coordinate Subscripting
* Supports Python 3.8, dropping support fpr Python 3.5
* numpy moveaxis replaced with sympy permutedims
* name parameter in Metric Tensor
* Tags to Tensor names

Contributors
............

This is the complete list of the people that contributed to this release, with a + sign indicating first contribution.

* Shreyas Bapat
* Ritwik Saha
* Manvi Gupta
* Micky Yun Chan+
* DylanBrdt+ (GitHub Username)
* Vineet Gandham+
* Pratyush Kerhalkar+
* Bhavam Vidyarthi+
* Khalid Shaikh+
* Rohit Sanjay+
* Saurabh+
* Raahul Singh+
* Nimesh Vashishtha+
* Shamanth R Nayak K+
* Arnav Das+
* Gim Seng Ng+
* Nihar Gupte+
* Suyash Salampuria+
* Atul Mangat+
* Ganesh Tarone+
* Shreyas Kalvankar+
* Swastik Singh+
* Jyotirmaya Shivottam+
* Sitara Srinivasan+
* Aayush Gautam+
* Zac Yauney+
* Gagan-Shenoy+
* Bibek Gautam+
* Erin Allard+
* Suyog Garg+


einsteinpy 0.2.1 - 2019-11-02
-----------------------------

This minor release would bring improvements and new feature additions to the already existing symbolic calculations module along
with performance boosts of order of 15x.

This release concludes the SOCIS 2019 projects of Sofía Ortín Vela (ortinvela.sofia@gmail.com) and Varun Singh(varunsinghs2021@gmail.com).

Part of this release is sponsored by European Space Agency, through Summer of Code in Space
(SOCIS) 2019 program.

Features
........

* New tensors in symbolic module

  * Ricci Scalar
  * Weyl Tensor
  * Stress-Energy-Momentum Tensor
  * Einstein Tensor
  * Schouten Tensor

* Improvement in performance of current tensors
* Lambdify option for tensors
* Support for vectors at arbitrary space-time symbolically as 1\ :sup:`st` order tensor.
* Support for scalars at arbitrary space-time symbolically as 0\ :sup:`th` order tensor.
* Addition of constants sub-module to symbolic module
* Improvement in speed of Geodesic plotting
* Move away from Jupyter and Plotly Widgets
* New Plotting Framework

Contributors
............

This is the complete list of the people that contributed to this release, with a + sign indicating first contribution.

* Shreyas Bapat
* Ritwik Saha
* Sofía Ortín Vela
* Varun Singh
* Arnav Das+
* Calvin Jay Ross+  


einsteinpy 0.2.0 - 2019-07-15
-----------------------------

This release brings a lot of new features for the EinsteinPy Users. 

A better API, intuitive structure and easy coordinates handling! This major release
comes before Python in Astronomy 2019 workshop and brings a lots of cool stuff. 

Part of this release is sponsored by ESA/ESTEC - Adv. Concepts & Studies Office
(European Space Agency), through Summer of Code in Space (SOCIS) 2019 program.

This is a short-term supported version and will be supported only until December 2019. 
For any feature request, write a mail to developers@einsteinpy.org describing what you need.

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
* Numerical Calculation and Symbolic Manipulation of Christoffel Symbols
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
..........

* `Issue #115`_: Coordinate Conversion had naming issues that made them confusing!
* `Issue #185`_: Isort had conflicts with Black
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
