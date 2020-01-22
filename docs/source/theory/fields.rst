==================================================
Gravitation and Electromagnetism as a Field Theory
==================================================

The action for macroscopic gravity, electromagnetism and (possibly) charged
relativistic dust is:

.. math::

    S = S_H + S_M  + S_{EM} + S_q

where:

.. math::

    S_H[g^{\mu\nu}] = {c^4\over 16\pi G} \int R \sqrt{-g}  d^4 x

    S_M[g^{\mu\nu}, x^{\mu}] = -c\int \rho \sqrt{v_\mu v^\mu} \sqrt{-g}  d^4 x

    S_{EM}[g^{\mu\nu}, A^{\mu}] =
        -{1\over4\mu_0} \int F_{\alpha\beta} F^{\alpha\beta} \sqrt{-g}  d^4 x

    S_q[x^\mu, A^{\mu}] = -\int \rho_{EM} v^\mu A_\mu \sqrt{-g}  d^4 x

where :math:`x^\mu` is the field of the matter, :math:`A^{\mu}` is the electromagnetic
field and :math:`g^{\mu\nu}` is the gravitational field. We vary with respect to each
of them to obtain (interacting) equations of motion. :math:`c` is the speed of light,
:math:`G` is the gravitational constant, :math:`\mu_0` the permeability of vacuum. :math:`\rho`
is the mass density of the dust, :math:`\rho_{EM}` is the charge density of the dust,
:math:`v^\mu={ d x^\mu\over \d\tau}` is 4-velocity of the dust,
:math:`F_{\alpha\beta}=\nabla_\alpha A_\beta-\nabla_\beta A_\alpha`
is the electromagnetic field tensor, :math:`R` is the Ricci scalar.

Gravitation
===========

We vary with respect to :math:`g^{\mu\nu}`. By changing the metric, we also change
the invariant volume element (thus also :math:`\rho`), so we need to be careful to
vary properly. We start with :math:`S_H`:

.. math::

    \delta S_H = \delta {c^4\over 16\pi G} \int R \sqrt{-g}  d^4 x =

        = {c^4\over 16\pi G} \int
            (\delta g^{\mu\nu}) R_{\mu\nu} \sqrt{-g}
            +g^{\mu\nu} (\delta R_{\mu\nu}) \sqrt{-g}
            +R (\delta \sqrt{-g})
             d^4 x=

        = {c^4\over 16\pi G} \int
            (\delta g^{\mu\nu}) R_{\mu\nu} \sqrt{-g}
            +g^{\mu\nu} \left(
                \nabla_\lambda(\delta \Gamma^\lambda_{\nu\mu})
                -\nabla_\nu(\delta \Gamma^\lambda_{\lambda\mu})
                \right)\sqrt{-g}
            +R (
            -\frac{1}{2} \sqrt{-g}\, g_{\mu\nu} (\delta g^{\mu\nu}))
             d^4 x=

        = {c^4\over 16\pi G} \int
            (\delta g^{\mu\nu}) R_{\mu\nu} \sqrt{-g}
            + \left(
                \nabla_\lambda g^{\mu\nu}(\delta \Gamma^\lambda_{\nu\mu})
                -\nabla_\nu g^{\mu\nu}(\delta \Gamma^\lambda_{\lambda\mu})
                \right)\sqrt{-g}
            -\frac{1}{2} R g_{\mu\nu} \sqrt{-g}\,
                (\delta g^{\mu\nu})
             d^4 x=

        = {c^4\over 16\pi G} \int
            (\delta g^{\mu\nu}) R_{\mu\nu} \sqrt{-g}
            -\frac{1}{2} R g_{\mu\nu} \sqrt{-g}\,
                (\delta g^{\mu\nu})
             d^4 x=

        = {c^4\over 16\pi G} \int \left( R_{\mu\nu} -\frac{1}{2} R g_{\mu\nu} \right)
                (\delta g^{\mu\nu}) \sqrt{-g}
             d^4 x

Variation of :math:`S_M` is:

.. math::

    \delta S_M = -c \delta \int \rho \sqrt{v_\mu v^\mu} \sqrt{-g}  d^4 x =

        = - \delta \int c \sqrt{ {p}_\mu  {p}^\mu}  d^4 x =

        = - \int c {\delta(g^{\mu\nu}  {p}_\mu  {p}_\nu)
            \over 2\sqrt{ {p}_\alpha  {p}^\alpha}}  d^4 x =

        = - \int c {  {p}_\mu  {p}_\nu
            \over 2\sqrt{ {p}_\alpha  {p}^\alpha}}
            \delta(g^{\mu\nu}) d^4 x =

        = - \int c { \rho v_\mu \rho v_\nu
            \sqrt{-g}^2
            \over 2 \rho c \sqrt{-g} }
             \delta(g^{\mu\nu}) d^4 x =

        = - \int \frac{1}{2} \rho v_\mu v_\nu
             \delta(g^{\mu\nu}) \sqrt{-g}  d^4 x

The variation of :math:`S_{EM}` is:

.. math::

    \delta S_{EM} = -\delta \int {1\over 4\mu_0} F_{\alpha\beta} F^{\alpha\beta}
            \sqrt{-g} d^4 x =

        = -\delta \int {1\over 4\mu_0} g^{\alpha\lambda} g^{\beta\rho}
            F_{\alpha\beta} F_{\lambda\rho} \sqrt{-g} d^4 x =

        = -{1\over 4\mu_0} \int  \left(\delta (g^{\alpha\lambda} g^{\beta\rho})
            F_{\alpha\beta} F_{\lambda\rho} \sqrt{-g}
            + g^{\alpha\mu} g^{\beta\rho}
            F_{\alpha\beta} F_{\lambda\rho} \left(\delta \sqrt{-g}
            \right)
            \right) d^4 x =

        = -{1\over 4\mu_0} \int  \left(2(\delta g^{\alpha\lambda}) g^{\beta\rho}
            F_{\alpha\beta} F_{\lambda\rho} \sqrt{-g}
            + g^{\alpha\lambda} g^{\beta\rho}
            F_{\alpha\beta} F_{\lambda\rho} \left(-\frac{1}{2} \sqrt{-g}
            g_{\mu\nu} (\delta g^{\mu\nu})
            \right)
            \right) d^4 x =

        = -{1\over 4\mu_0} \int  \left(2(\delta g^{\alpha\lambda})
            F_{\alpha\beta} F_\lambda{}^\beta
            -\frac{1}{2} F_{\alpha\beta} F^{\alpha\beta}
            g_{\mu\nu} (\delta g^{\mu\nu})
            \right) \sqrt{-g}  d^4 x =

        = -{1\over 2\mu_0} \int  \left(
            F_{\mu\beta} F_\nu{}^\beta
            -{1\over 4} F_{\alpha\beta} F^{\alpha\beta}
            g_{\mu\nu}
            \right) (\delta g^{\mu\nu}) \sqrt{-g}  d^4 x

The variation of :math:`\delta S_q=0`.

The equations of motion are:

.. math::

    {c^4\over 16\pi G} \left( R_{\mu\nu} -\frac{1}{2} R g_{\mu\nu} \right)
        - \frac{1}{2} \rho v_\mu v_\nu
        -{1\over 2\mu_0} \left(
            F_{\mu\beta} F_\nu{}^\beta
            -{1\over 4} F_{\alpha\beta} F^{\alpha\beta}
            g_{\mu\nu}
            \right) = 0

We rearrange:

.. math::

    R_{\mu\nu} -\frac{1}{2} R g_{\mu\nu}
        =
        {8\pi G\over c^4} \rho v_\mu v_\nu
        +{8\pi G\over c^4}{1\over\mu_0} \left(
            F_{\mu\beta} F_\nu{}^\beta
            -{1\over 4} F_{\alpha\beta} F^{\alpha\beta}
            g_{\mu\nu}
            \right)

We define the stress energy tensor as:

.. math::
    :label: einstein-eq

    R_{\mu\nu} -\frac{1}{2} R g_{\mu\nu} = {8\pi G\over c^4} T_{\mu\nu}

.. math::
    :label: stress-energy-tensor-formula

    T_{\mu\nu} = - {2\over\sqrt{-g}}{\delta (S_M + S_{EM} + S_q)
        \over \delta g^{\mu\nu}}

And we get:

.. math::
    :label: stress-energy-tensors

    T_{\mu\nu} = T^M_{\mu\nu} + T^{EM}_{\mu\nu}

    T^M_{\mu\nu} = \rho v_\mu v_\nu

    T^{EM}_{\mu\nu} = {1\over \mu_0} \left(
            F_{\mu\beta} F_\nu{}^\beta
            -{1\over 4} F_{\alpha\beta} F^{\alpha\beta}
            g_{\mu\nu} \right)

The equations :eq:`einstein-eq` are called Einstein's equations and
the equations :eq:`stress-energy-tensors` are stress energy tensors for the
relativistic dust and electromagnetism.
The equation :eq:`stress-energy-tensor-formula` is the stress
energy tensor corresponding to the given action.
Sometimes it is not possible to write
an action for more complex matter (perfect fluid,
Navier-Stokes equations for fluid, ...) in which case we cannot
use :eq:`stress-energy-tensor-formula`, but we can still specify the
stress energy tensor directly and :eq:`einstein-eq` are the equations of motion.

Electromagnetism
================

We vary with respect to :math:`A^\mu`.
The variation of :math:`\delta S_H=0`. The variation of :math:`\delta S_M = 0`.
The variation of :math:`S_{EM}` is:

.. math::

    \delta S_{EM} =

        = -{1\over4\mu_0} \delta \int F_{\mu\nu} F^{\mu\nu}
            \sqrt{-g}  d^4 x =

        = -{1\over2\mu_0} \int F^{\mu\nu} (\delta F_{\mu\nu})
            \sqrt{-g}  d^4 x =

        = -{1\over\mu_0} \int F^{\mu\nu} (\delta \partial_\nu A_\mu)
            \sqrt{-g}  d^4 x =

        = -{1\over\mu_0} \int F^{\mu\nu} \partial_\nu (\delta A_\mu)
            \sqrt{-g}  d^4 x =

        = {1\over\mu_0} \int \partial_\nu (F^{\mu\nu}\sqrt{-g})
            (\delta A_\mu)  d^4 x =

        = {1\over\mu_0} \int \left(
            {1\over\sqrt{-g}}\partial_\nu (F^{\mu\nu}\sqrt{-g}) \right)
            (\delta A_\mu) \sqrt{-g} d^4 x =

        = {1\over\mu_0} \int \nabla_\mu F^{\mu\nu} (\delta A_\nu)
            \sqrt{-g}  d^4 x

The variation of :math:`S_q` is:

.. math::

    \delta S_q =

        =-\delta\int \rho_{EM} v^\nu A_\nu \sqrt{-g}  d^4 x =

        =-\int \rho_{EM} v^\nu (\delta A_\nu) \sqrt{-g}  d^4 x =

The equation of motion is:

.. math::

    {1\over\mu_0} \nabla_\mu F^{\mu\nu} - \rho_{EM} v^\nu = 0

Rearranging:

.. math::

    \nabla_\mu F^{\mu\nu} = \mu_0 \rho_{EM} v^\nu

Relativistic Dust
=================

We vary the whole action with respect to :math:`x^\mu`.
The variation of :math:`\delta S_H=0`.
The variation of :math:`S_M` is:

.. math::

    \delta S_M
        = -c\delta \int \rho \sqrt{v_\mu v^\mu} \sqrt{-g}  d^4 x =

        = - \delta \int c \sqrt{ {p}_\mu  {p}^\mu}  d^4 x =

        = - \int c {\delta(g^{\mu\nu}  {p}_\mu  {p}_\nu)
            \over 2\sqrt{ {p}_\alpha  {p}^\alpha}}  d^4 x =

        = - \int c { 2 g^{\mu\nu}  {p}_\mu (\delta  {p}_\nu)
            \over 2\sqrt{ {p}_\alpha  {p}^\alpha}}  d^4 x =

        = - \int c {  {p}_\mu \over \sqrt{ {p}_\alpha  {p}^\alpha}}
            (\delta  {p}^\mu)  d^4 x =

        = - \int c {  {p}_\mu \over \sqrt{ {p}_\alpha  {p}^\alpha}}
            \partial_\nu \left( {p}^\nu(\delta x^\mu) -  {p}^\mu (\delta x^\nu)\right)
                 d^4 x =

        = \int c \partial_\nu \left({  {p}_\mu \over
        \sqrt{ {p}_\alpha  {p}^\alpha}}
                \right)
            \left( {p}^\nu(\delta x^\mu) -  {p}^\mu (\delta x^\nu)\right)
                 d^4 x =

        = \int c \left(
            \partial_\nu \left({  {p}_\mu \over \sqrt{ {p}_\alpha  {p}^\alpha}} \right)
            -\partial_\mu \left({  {p}_\nu \over \sqrt{ {p}_\alpha  {p}^\alpha}} \right)
            \right)
             {p}^\nu(\delta x^\mu)
                 d^4 x =

        = \int c \left(
            \nabla_\nu \left({  {p}_\mu \over \sqrt{ {p}_\alpha  {p}^\alpha}} \right)
            -\nabla_\mu \left({  {p}_\nu \over \sqrt{ {p}_\alpha  {p}^\alpha}} \right)
            \right)
             {p}^\nu(\delta x^\mu)
                 d^4 x =

        = \int \left( \nabla_\nu v_\mu -\nabla_\mu v_\nu \right)
            \rho v^\nu (\delta x^\mu) \sqrt{-g}
                 d^4 x =

        = \int \rho (\nabla_\nu v_\mu) v^\nu (\delta x^\mu) \sqrt{-g}
                 d^4 x

The variation of :math:`\delta S_{EM}=0`. The variation of :math:`S_q` is:

.. math::

    \delta S_q
        = - \delta \int \rho_{EM} v^\mu A_\mu \sqrt{-g}  d^4 x =

        = - \delta \int  {j}^\mu A_\mu  d^4 x =

        = - \int (\delta  {j}^\mu) A_\mu  d^4 x =

        = - \int \partial_\nu \left( {j}^\nu (\delta x^\mu)
            -  {j}^\mu (\delta x^\nu)\right) A_\mu  d^4 x =

        = \int \left( {j}^\nu (\delta x^\mu)
            -  {j}^\mu (\delta x^\nu)\right) \partial_\nu A_\mu  d^4 x =

        = \int  {j}^\nu (\delta x^\mu) (\partial_\nu A_\mu -\partial_\mu A_\nu)
             d^4 x =

        = \int \rho_{EM} v^\nu  (\nabla_\nu A_\mu -\nabla_\mu A_\nu)
            (\delta x^\mu) \sqrt{-g}
             d^4 x =

        = -\int \rho_{EM} v^\nu  F_{\mu\nu} (\delta x^\mu) \sqrt{-g}
             d^4 x

The equation of motion is:

.. math::

    \rho (\nabla_\nu v_\mu) v^\nu
        -\rho_{EM} v^\nu  F_{\mu\nu} = 0

Rearranging:

.. math::

    \rho (\nabla_\nu v_\mu) v^\nu = \rho_{EM} v^\nu  F_{\mu\nu}

This is the geodesic equation with Lorentz force.

Equations of Motion
===================

All together, the equations of motion are:

.. math::

    R_{\mu\nu} -\frac{1}{2} R g_{\mu\nu}
        =
        {8\pi G\over c^4} \rho v_\mu v_\nu
        +{8\pi G\over c^4}{1\over\mu_0} \left(
            F_{\mu\beta} F_\nu{}^\beta
            -{1\over 4} F_{\alpha\beta} F^{\alpha\beta}
            g_{\mu\nu}
            \right)

    \nabla_\mu F^{\mu\nu} = \mu_0 \rho_{EM} v^\nu

    \rho (\nabla_\nu v_\mu) v^\nu = \rho_{EM} v^\nu  F_{\mu\nu}

The first equation determines :math:`g_{\mu\nu}` from the given sources (the stress
energy tensors) on the right hand side, that depend on :math:`\rho`, :math:`v^\mu`,
:math:`A^\mu` and :math:`g_{\mu\nu}`. The second equation
determines :math:`A^\mu` from the sources (`\rho_{EM}` and :math:`v^\mu`) and from
:math:`g_{\mu\nu}` (through the covariant derivative).
Finally, the last equation determines :math:`x^\mu` and :math:`v^\mu` from the given fields
:math:`A^\mu` (through the electromagnetic field tensor) and :math:`g_{\mu\nu}` (through
the covariant derivative).

Conservation
------------

We apply covariant 4-divergence and use Bianci identities on the first
equation:

.. math::

    0 = \nabla_\mu T^{\mu\nu} = \nabla_\mu (T^{\mu\nu}_M + T^{\mu\nu}_{EM})

So the total stress energy tensor is conserved. This fact makes the equations
of motion (that follow from the action principle) not all independent. The
third equation can be derived from the fist two as follows.

We calculate:

.. math::

    \nabla_\mu T^{\mu\nu}_{M} = \nabla_\mu (\rho v^\mu v^\nu)

    \nabla_\mu T^{\mu\nu}_{EM} = F^{\alpha\nu} \rho_{EM} v_\alpha

and we get:

.. math::

    \nabla_\mu (\rho v^\mu v^\nu) + F^{\alpha\nu} \rho_{EM} v_\alpha = 0

    \nabla_\mu (\rho v^\mu) v^\nu
    + \rho v^\mu \nabla_\mu v^\nu
        + F^{\alpha\nu} \rho_{EM} v_\alpha = 0

The first term vanishes, because:

.. math::

    v_\nu \nabla_\mu (\rho v^\mu) v^\nu
    + v_\nu \rho v^\mu \nabla_\mu v^\nu
        + v_\nu F^{\alpha\nu} \rho_{EM} v_\alpha = 0

    v_\nu \nabla_\mu (\rho v^\mu) v^\nu
        + v_\nu F^{\alpha\nu} \rho_{EM} v_\alpha = 0

    c^2 \nabla_\mu (\rho v^\mu) + v_\nu F^{\alpha\nu} \rho_{EM} v_\alpha = 0

    c^2 \nabla_\mu (\rho v^\mu) = 0

where we used :math:`v_\nu \nabla_\mu v^\nu=0` (follows from differentiating
:math:`c^2 = v_\nu v^\nu`)
and :math:`v_\nu F^{\alpha\nu}
v_\alpha=0` (contracting symmetric and antisymmetric tensors). We are left
with:

.. math::

    \rho v^\mu \nabla_\mu v^\nu + F^{\alpha\nu} \rho_{EM} v_\alpha = 0

    \rho v^\mu \nabla_\mu v^\nu = -F^{\alpha\nu} \rho_{EM} v_\alpha

    \rho v^\mu \nabla_\mu v^\nu = F^{\nu\alpha} \rho_{EM} v_\alpha

Which is the third equation.
