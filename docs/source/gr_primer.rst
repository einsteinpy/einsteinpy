General Relativity Primer
=========================

This page provides a high-level overview of some of the concepts in 
General Relativity. For more detailed treatments, see the `references`_ 
at the end of this page.

Einstein's Field Equations
--------------------------
Einstein's Field Equations (EFE) relate local spacetime curvature 
with local energy and momentum. In short, these equations determine the metric tensor 
of a spacetime, given the arrangement of stress-energy. The EFE are given as follows:

.. math::
    \boxed{R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}}.

Here, :math:`R_{\mu\nu}` is the Ricci tensor, :math:`R` is the 
scalar curvature (trace of Ricci tensor), :math:`g_{\mu\nu}` 
is the metric tensor, :math:`\Lambda` is the cosmological constant and 
lastly, :math:`T_{\mu\nu}` denotes the stress-energy tensor. 
All other variables hold their usual meaning. If we introduce the 
Einstein tensor :math:`G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu}`, 
then the EFE can be written as:

.. math::
    \boxed{G_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu} - \Lambda g_{\mu\nu}}.

These equations form a ten-component tensor equation, which collectively 
denotes a system of coupled non-linear partial differential equations. These are 
usually intractable to approach analytically. However, under certain conditions, 
the EFE can be simplified and solved to yield *exact* solutions. 
For example, if the Einstein tensor is assumed to vanish, which implies that the
stress-energy tensor also vanishes, then the EFE reduces to the vacuum Einstein 
equations, whose solutions are called *vacuum solutions*. The Minkowski, 
Schwarzschild and Kerr spacetimes are some examples of vacuum solutions.

Similarly, if an electromagnetic field is assumed to be present and also 
the only source of non-gravitational energy, then the EFE reduces to the 
source-free Maxwell-Einstein equations, whose exact solutions are termed 
*electrovacuum solutions*. The Kerr-Newman and Reissner-Nordström spacetimes are 
examples of electrovacuum solutions.

Metric Tensor
-------------
The metric tensor denotes a solution to the EFE. As such, it is a fundamental 
entity in general relativity, that captures the geometry of spacetime. It is 
a symmetric, indefinite rank-2 tensor, which can be represented by a matrix. 
It can be thought of as characterizing the differential line element for 
a given geometry:

.. math::
  \mathrm{d} s^2 = g_{\mu\nu}\mathrm{d}x^{\mu}\mathrm{d}x^{\nu}

To get the matrix representation, the coefficients, :math:`g_{\mu\nu}`, in the above equation 
are substituted into the corresponding places designated by the indices, 
:math:`\mu` and :math:`\nu`, in a matrix. For example, the metric tensor in the spherical-polar 
coordinate system can be written as:

.. math::

  g_{\mu\nu} = \begin{pmatrix}
    1 & 0 & 0 \\
    0 & r^2 & 0 \\
    0 & 0 & r^2\sin^2\theta
  \end{pmatrix}

Note that the off-diagonal components are 0, since we are using an 
orthogonal basis. However, it is not always the case and in general, 
:math:`g_{\mu\nu} \ne 0`, for :math:`\mu \ne \nu`.

The metric tensor is also used to define the operations of lowering 
and raising indices. For this, we also need the *inverse* metric tensor, :math:`g^{\mu\nu}` 
(or *contravariant* metric tensor to be precise), which is defined using the following identity:

.. math::
  g^{\mu\nu}g_{\nu\lambda} = \delta^{\mu}_{\lambda}

The inverse metric tensor can then be used to raise indices, while the metric tensor 
can be used to lower indices.

Notion of Curved Space
----------------------
Imagine a bug traveling on a sheet of paper folded into a cone. The 
bug can't see up and down. So, it lives in a 2D world. But it still 
experiences the curvature of the surface (space), as, after a long journey, 
it would return to the starting position.

Mathematically, the curvature of space is characterized by the rank-4 
Riemann Curvature tensor, whose *contraction* gives the rank-2 Ricci 
tensor. Taking the trace of the Ricci tensor yields the rank-0 Ricci 
Scalar or Scalar Curvature.

Straight lines in Curved Space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Straight lines are used to describe the shortest path between two points in Euclidean 
space. *Geodesics* extend this notion of "shortest path" to curved spaces. In GR, geodesics 
represent the shortest path between two points in a possibly curved spacetime. They are 
completely characterized by the metric tensor and its derivatives, resulting in a set 
of coupled ordinary differential equations:

.. math::
  \boxed{\frac{\mathrm{d}^2x^{\mu}}{\mathrm{d}s^2} + \Gamma^{\mu}_{\alpha\beta}\frac{\mathrm{d}x^{\alpha}}{\mathrm{d}s}\frac{\mathrm{d}x^{\beta}}{\mathrm{d}s} = 0}

Here, :math:`x` denotes the geodesic and the derivatives are taken with respect to :math:`s`, an 
*affine* parameter that uniquely parameterizes :math:`x`. :math:`\Gamma^{\mu}_{\alpha\beta}` denotes 
*Christoffel symbols of the second kind*, which is essentially an array of partial derivatives of the 
metric tensor, computed with respect to the coordinate basis. It can be thought of as a "connection" 
between the derivates of nearby points along the geodesic and is given as:

.. math::
  \Gamma^{\mu}_{\alpha\beta} = \frac{1}{2}g^{\mu\gamma}\left(\frac{\partial g_{\alpha\beta}}{\partial x^{\gamma}} + \frac{\partial g_{\alpha\gamma}}{\partial x^{\beta}} - \frac{\partial g_{\beta\gamma}}{\partial x^{\alpha}}\right)

The Riemann Curvature tensor encapsulates the idea of *curvature* in GR. It can be written in 
a condensed notation using the Christoffel symbols and their derivaties:

.. math::
  R^{\mu}_{\alpha\beta\gamma} = \partial_\gamma\Gamma^{\mu}_{\alpha\beta} - \partial_\beta\Gamma^{\mu}_{\alpha\gamma} + \Gamma^{\mu}_{\alpha\delta}\Gamma^{\delta}_{\beta\gamma} - \Gamma^{\mu}_{\beta\delta}\Gamma^{\delta}_{\alpha\gamma}

A space with zero curvature implies that the Riemann tensor is zero and vice-versa. Since, 
we are dealing with tensors, i.e. objects that operate independently of basis choice, this statement 
holds for any coordinate system, i.e. :math:`R = 0` in one coordinate system implies :math:`R = 0` 
and by extension, zero curvature in all coordinate systems.

The Ricci tensor is another geometrical object in GR that is related to curvature. It can be obtained 
by contracting the first and third indices of the Riemann tensor:

.. math::
  R_{\mu\nu} = R^{\rho}_{\mu\rho\nu}

:math:`R_{\mu\nu}` can be thought of as quantifying the deformation of a shape as it is 
translated along a given geodesic. The trace of the Ricci tensor gives the Scalar Curvature, :math:`R`:

.. math::
  R = g^{\mu\nu}R_{\mu\nu}

:math:`R` relates the volume of infinitesimal geodesic balls in curved space to that in Euclidean space.

This was a short and superficial look into some of the basic quantities that are used to characterize 
the structure of spacetime in General Relativity. Readers, who are interested in gaining a deeper 
understanding, are strongly recommended to peruse the resources listed in `References`_.

----

References
----------

* Wikipedia

  * `Einstein's Field Equations <https://en.wikipedia.org/wiki/Einstein%27s_field_equations>`_
  * `Metric Tensor <https://en.wikipedia.org/wiki/Metric_tensor>`_
  * `Raising and lowering indices <https://en.wikipedia.org/wiki/Raising_and_lowering_indices>`_
  * `Riemann Curvature Tensor <https://en.wikipedia.org/wiki/Riemann_curvature_tensor>`_
  * `Ricci Tensor <https://en.wikipedia.org/wiki/Ricci_tensor>`_
  * `Scalar Curvature <https://en.wikipedia.org/wiki/Scalar_curvature>`_
  * `Geodesic equation <https://en.wikipedia.org/wiki/Geodesics_in_general_relativity>`_
  * `Christoffel Symbols <https://en.wikipedia.org/wiki/Christoffel_symbols>`_
  * `Levi-Civita Connection <https://en.wikipedia.org/wiki/Levi-Civita_connection>`_

* General Relativity Textbooks (with links to public copies)

  * `Gravitation <https://archive.org/details/gravitation0000misn>`_ by Charles W. Misner, John Archibald Wheeler, and Kip Thorne 
  * `General Relativity <https://archive.org/details/generalrelativit0000wald>`_ by Robert Wald
  * `Spacetime and Geometry: An Introduction to General Relativity <https://worldcat.org/title/1112495919>`_ by Sean Carroll
  * `Black Hole Physics <https://archive.org/details/blackholephysics0000frol>`_ by Valeri P. Frolov and Igor D. Novikov
  * `General Relativity <http://lightandmatter.com/genrel/>`_ by Benjamin Crowell
  