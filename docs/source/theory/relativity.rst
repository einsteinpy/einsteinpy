.. index:: Relativity

Relativity
==========

Introduction: Why Tensors
-------------------------

This section gives a brief introduction, and in the next sections we derive
everything in detail. The Newton law is:

.. math::

    m {\d^2 {\bf x}\over\d t^2} = {\bf F}

and using a potential for ${\bf F}$, we get:

.. math::

    {\d^2 {\bf x}\over\d t^2} = -\nabla \phi

    {\d^2 x^i\over\d t^2} = -\partial^i \phi

the last two equations are two different equivalent ways to write a tensor
equation in 3D, which means that this equation has the exact same form (is
valid) in any (spatial) coordinate system (rotated, translated, in cartesian
coordinates, spherical coordinates, ...). Each coordinate system has a
different metric, but we can always locally transform into
$g_{ij}=\diag(1, 1, 1)$.

However, if our coordinate transformation depends on time (e.g. a rotating
disk), then the above tensor equation changes (e.g. for the rotating disk, we
get the Coriolis acceleration term), that's because time is treated as a
parameter, not as a coordinate.

To fix this, we need to work in 4D and treat time as a coordinate, so we
introduce $x^0 = ct$ where $c$ is any constant speed (it can be any speed,
doesn't have to be the speed of light). Then in 4D, the above equations are not
tensor equations anymore, because the operator ${\d\over \d t} = c \partial_0$
is not a tensor. The 4D tensor formulation happens to be the geodesic equation:

.. math::

    {\d x^\beta\over\d\lambda}\nabla_\beta {\d x^\alpha\over\d\lambda} = 0

    R_{00} = 4\pi G\rho

    R_{ij} = 0

Which (given that we know how to calculate the Ricci tensor in our coordinates)
is valid in any coordinates, not only rotated, translated, cartesian,
spherical, ..., but also with arbitrary time dependence, e.g. a rotating disk,
accelerating disk, ...

After suitable local coordinate transformation, we can only get two possible
metrics (that connect the time and spatial coordinates): $\diag(-1, 1, 1, 1)$
and $\diag(1, 1, 1, 1)$. Inertial systems have no fictitious forces, so the
metrics is one of the two above (possibly with $c\to\infty$). Transformation
between inertial systems is such a coordinate transformation that leaves the
metric intact, e.g.:

.. math::

     g' = \Lambda^T g \Lambda

There is no coordinate transformation that turns the metric $\diag(-1, 1, 1,
1)$ into $\diag(1, 1, 1, 1)$, so we need to choose either one to describe one
inertial system and then all other inertial systems will automatically have a
metric with the same signature.

The Newton law is valid for small speeds compared to the speed of light, so
when we want to extend the theory for all speeds, we only have 4 options: O(3,
1) with either $c\to\infty$ or $c$ finite and O(4) with either $c\to\infty$ or
$c$ finite. If $c$ is finite, it has to be large enough, so that we still
recover the Newton law for small speeds with the given experimental precision.
All 4 cases give the correct Newton law, but give different predictions for
large speeds. All we need to do to decide which one is correct is to perform
such large speeds (relativistic) experiments. It turns out that all such
relativistic experiments are in agreement with the O(3, 1) case where $c$ is
the (finite) speed of light and with disagreement with the 3 other cases. For
small speeds however (i.e. Newtonean physics), all 4 cases will work, as long
as $c$ is chosen large enough.

Given a tensor equation, we can easily determine, if it transforms correctly
under the Galilean ($c\to\infty$) or Lorentz transformations ($c$ is finite).
All we have to do is to perform the limit $c\to\infty$. For example the Newton
second law is recovered if we do the $c\to\infty$ limit, but Maxwell equations
are only recovered if we choose $c$ to be exactly the speed of light in the
Maxwell equations.

The reason why we write equations as tensor equations in 4D is that we can then
use any coordinates (including any time dependence), i.e. any observer, and the
equations still have the exact same form. So specifying the metrics is enough
to define the coordinates (observer) and since the equations has only one form,
that is all we need. If we write equations only as tensors in 3D, we not only
need to specify the (3D) metrics, but also how the observer accelerates with
respect to some (usually inertial) frame where the equations (let's say Newton
law) is defined and we then need to transform all the time derivatives
correctly. By using tensors in 4D, all those transformations are taken care of
by the standard tensor machinery and all we need to care about is exactly one
observer, defined by its metric tensor.

By choosing the correct metrics and $c$ (i.e. $\diag(-1, 1, 1, 1)$ and $c$ the
speed of light), all equations are then automatically Lorentz invariant. If we
choose $c\to\infty$ (and any metric), we automatically get all equations
Galilean invariant.


High School Formulation
-----------------------


The usual (high school) formulation is the second Newton's law:

.. math::

    m {\d^2 {\bf x}\over\d t^2} = {\bf F}

for some particle of the mass $m$ and position ${\bf x}$. To determine the
force $\bf F$, we have at hand the Newton's law of
gravitation:

.. math::

    |{\bf F}| = G {m_1 m_2\over r^2}

that expresses the magnitude $|\bf F|$ of the force between two particles with
masses $m_1$ and $m_2$ and we also know that the direction of the force is
directly towards the other particle.
We need to take into account all particles in the system, determine the
direction and magnitude of the force due to each of them and sum it up.

College Formulation
-------------------


Unfortunately, it is quite messy to keep track of the direction of the forces
and all the masses involved, it quickly becomes cumbersome for more than 2
particles. For this reason, the better approach is to calculate the force
(field) from the mass density function $\rho$:

.. math::

    \nabla\cdot{\bf F} = -4\pi Gm\,\rho(t, x, y, z)

To see that both formulations are equivalent,
integrate both sides inside some sphere:

.. math::

    \int\nabla\cdot{\bf F}\,\d x\d y\d z = -4\pi Gm_2\int\rho\,\d x\d y\d z

apply the Gauss theorem to the left hand side:

.. math::

    \int\nabla\cdot{\bf F}\,\d x\d y\d z = \int{\bf F}\cdot{\bf n}\,\d S= 4\pi r^2\,{\bf F}\cdot{\bf n}

where ${\bf n}={{\bf r}\over |{\bf r}|}$ and
the right hand side is equal to $-4\pi G m_1m_2$ and we get:

.. math::

    {\bf F}\cdot{\bf n} = -G{m_1m_2\over r^2}

now we multiply both sides with ${\bf n}$, use the fact that
$({\bf F}\cdot{\bf n}){\bf n} ={\bf F}$ (because ${\bf F}$ is spherically
symmetric), and we get the traditional Newton's
law of gravitation:

.. math::

    {\bf F} = -G{m_1m_2\over r^2}{\bf n}


It is useful to deal with a scalar field instead of a vector field (and also
not to have the mass $m$ of the test particle in our equations explicitly), so we
define a gravitational potential by:

.. math::

    {\bf F} = -m\nabla\phi(t, x, y, z)

then the law of gravitation is

.. math::
    :label: grav

    \nabla^2\phi = 4\pi G\rho

and the second law is:

.. math::

    m{\d^2 {\bf x}\over\d t^2} = -m\nabla\phi(t, x, y, z)


Note about units:

.. math::

    [r] = [{\bf x}] = \rm m


.. math::

    [m] = \rm kg


.. math::

    [\rho] = \rm kg\,m^{-3}


.. math::

    [F] = \rm kg\,m\,s^{-2}


.. math::

    [G] = \rm kg^{-1}\,m^3\,s^{-2}


.. math::

    [\phi] = \rm m^2\,s^{-2}

Example
~~~~~~~

Calculate the force acting on a test particle inside an infinitely thin
spherical shell of radius $R$ and surface mass distribution $\sigma(\theta,
\phi)=1$. We need to solve

.. math::
    :label: grav-example1

    \nabla^2\phi = 4\pi G\rho

with

.. math::

    \rho(x, y, z) = \sigma(\theta, \phi) {\delta(R-r)\over r^2}

    r = \sqrt{x^2 + y^2 + z^2}

the Green function of :eq:`grav-example1` is

.. math::

    G({\bf x}, {\bf y}) = {1\over |{\bf x} - {\bf y}|}

so the solution is:

.. math::

    \phi = \int G({\bf x}, {\bf y}) 4\pi G \rho({\bf y}) \d^3 y
        = 4\pi G \int {\rho({\bf y})\over |{\bf x} - {\bf y}|} \d^3 y
        =

    = 4\pi G \int {\sigma(\theta, \phi){\delta(R-r)\over r^2} r^2\sin\theta
        \over \sqrt{
            (x-r\sin\theta\cos\phi)^2 +
            (y-r\sin\theta\sin\phi)^2 +
            (z-r\cos\theta)^2
            }} \d \theta \d \phi \d r =

    = 4\pi G \int {\delta(R-r)\sin\theta
        \over \sqrt{
            (x-r\sin\theta\cos\phi)^2 +
            (y-r\sin\theta\sin\phi)^2 +
            (z-r\cos\theta)^2
            }} \d \theta \d \phi \d r =

    = 4\pi G \int {\sin\theta
        \over \sqrt{
            (x-R\sin\theta\cos\phi)^2 +
            (y-R\sin\theta\sin\phi)^2 +
            (z-R\cos\theta)^2
            }} \d \theta \d \phi =

    = 4\pi G \int {\sin\theta
        \over \sqrt{x^2 + y^2 + z^2 + R^2
            -2R(x\sin\theta\cos\phi + y\sin\theta\sin\phi + z\cos\theta)
            }} \d \theta \d \phi

for symmetry reasons we can set $x=0$, $y=0$ (it can also be done more exactly,
as shown in :ref:`spherical-int-example`):


.. math::

    \phi(0, 0, z)
    = 4\pi G \int_0^{2\pi} \d\phi \int_0^\pi \d\theta {\sin\theta
        \over \sqrt{z^2 + R^2 -2Rz\cos\theta }} =

    = 8\pi^2 G \int_0^\pi \d\theta {\sin\theta
        \over \sqrt{z^2 + R^2 -2Rz\cos\theta }} =

    = 8\pi^2 G \int_{-1}^1 {\d y \over \sqrt{z^2 + R^2 -2Rzy }} =

    = -{4\pi^2 G\over R z} \int_{(R-z)^2}^{(R+z)^2} {\d u \over \sqrt{u}} =

    = -{4\pi^2 G\over R z} \Big[2\sqrt u\Big]_{(R-z)^2}^{(R+z)^2} =

    = -{4\pi^2 G\over R z} \Big[2|R+z| - 2|R-z|\Big] =

    = -{4\pi^2 G\over R z} \Big[4z\Big] =

    = -{16\pi^2 G\over R}

This must hold for all $x$ and $y$ (less than $R$), so:

.. math::

    \phi(x, y, z) = -{16\pi^2 G\over R}

And the force is

.. math::

    {\bf F} = -m\nabla\phi(t, x, y, z) = -m\nabla
        \left(-{16\pi^2 G\over R}\right) = 0

So the force acting on a test particle inside the shell is zero.

Differential Geometry Formulation
---------------------------------


There are still problems with this formulation, because it is not immediatelly
clear how to write those laws in other frames, for example rotating, or
accelerating -- one needs to employ nontrivial assumptions about the systems,
space, relativity principle and it is often a source confusion.
Fortunately there is a way out --- differential geometry. By reformulating the
above laws in the language of the differential geometry, everything will
suddenly be very explicit and clear. As an added bonus, because the special and
general relativity uses the same language, the real differences between all
these three theories will become clear.

We write $x, y, z$ and $t$ as components of one 4-vector

.. math::

     x^\mu = \mat{ct\cr x\cr y\cr z\cr}

In this section, you can imagine $c=1$, but we'll need it later, so we put it
in right now, so that we don't need to rederive all equations again.
Now we need to connect the Newtonian equations to geometry. To do that, we
reformulate the Newton's second law:

.. math::

     {\d^2 x^i\over\d t^2} + \delta^{ij}\partial_j\phi =0

by choosing a parameter $\lambda$ such, that ${\d^2 \lambda\over\d t^2}=0$,
so in general

.. math::

     \lambda = at+b

and

.. math::

    {\d^2\over\d t^2} = a^2{\d^2\over\d \lambda^2}

so

.. math::

     {\d^2 x^i\over\d\lambda^2} + {1\over a^2}\delta^{ij}\partial_j\phi =0

and using the relation ${\d \lambda\over\d a}=a$ we get

.. math::

     {\d^2 x^i\over\d\lambda^2} + \delta^{ij}\partial_j\phi \left({\d t\over\d\lambda}\right)^2 =0

So using $x^0$ instead of $t$, we endup with the following equations:

.. math::

    {\d^2x^0\over\d\lambda^2}=0

    {\d^2 x^i\over\d\lambda^2} + {1\over c^2}\delta^{ij}\partial_j\phi
        \left({\d x^0\over\d\lambda}\right)^2 =0

But this is exactly the geodesic equation for the following Christoffel symbols:

.. math::
    :label: Chris-newton

    \Gamma^i_{00} = {1\over c^2}\delta^{ij}\partial_j\phi

and all other components are zero.

In order to formulate the gravitation law, we now need to express
$\nabla^2\phi$ in terms of geometric quantities like
$\Gamma^\alpha_{\beta\gamma}$ or $R^\alpha{}_{\beta\gamma\delta}$.
We get the only nonzero components of
the Riemann tensor:

.. math::

    R^j{}_{0k0} = -R^j{}_{00k} = {1\over c^2}\delta^{ji}\partial_i\partial_k\phi

we calculate the $R_{\alpha\beta}$ by contracting:

.. math::

    R_{00} = R^\mu{}_{0\mu0} = R^i{}_{0i0} = {1\over c^2}\delta^{ij}\partial_i\partial_j\phi


.. math::

    R_{ij} = 0

comparing with :eq:`grav` we see that the Newton gravitation law is

.. math::

    R_{00} = {4\pi G\over c^2}\rho

    R_{ij} = 0


Thus we have reformulated the Newton's laws in a frame invariant way --- the
matter curves the geometry using the equations:

.. math::

    R_{00} = {4\pi G\over c^2}\rho

    R_{ij} = 0

from which one can (for example) calculate the Christoffel symbols and other
things. The particles then move on the geodesics:

.. math::

    {\d^2 x^\alpha\over\d\lambda^2} + \Gamma^\alpha_{\beta\gamma} {\d x^\beta\over\d\lambda}{\d x^\gamma\over\d\lambda} = 0

Both equations now have the same form in all coordinate systems (inertial or
not) and it is clear how to transform them --- only the Christoffel symbols
(and Ricci tensor) change and we have a formula for their transformation.

Obviously this works for any value of $c$ (as it cancels out in the final
equations of motion) and at this level we don't really need it yet, so we can
set $c=1$ and forget about it. In the next section we will need some constant
in the metric to send to infinity in order to obtain the correct Christoffel
symbols, and we can conveniently just use $c$. Later on we introduce special
relativity and we need to introduce a speed of light and it turns out that we
can again just use $c$ for that without any loss of generality.

Metrics
-------


There is a slight problem with the metrics --- it can be proven that there is
no metrics, that generates the Christoffel symbols above. However, it turns out
that if we introduce an invariant speed $c$ in the metrics, then calculate the
Christoffel symbols (thus they depend on $c$) and then do the limit
$c\to\infty$, we can get the Christoffel symbols above.

In fact, it turns out that there are many such metrics that generate the right
Christoffel symbols. Below we list several similar metrics and the
corresponding Christoffel symbols (in the limit $c\to\infty$), so that we can
get a better feeling what metrics work and what don't and why:

.. math::

    g_{\mu\nu} = \mat{-c^2-2\phi & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & -1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=-\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{-c^2-2\phi & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & -1 & 0\cr 0 & 0 & 0 & -1\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=-\partial_y\phi


.. math::

    \Gamma^3_{00}=-\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{-c^2-2\phi & 0 & 0 & 0\cr 0 & -1 & 0 & 0\cr 0 & 0 & -1 & 0\cr 0 & 0 & 0 & -1\cr}


.. math::

    \Gamma^1_{00}=-\partial_x\phi


.. math::

    \Gamma^2_{00}=-\partial_y\phi


.. math::

    \Gamma^3_{00}=-\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{-c^2+45-2\phi & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{-c^2-2\phi & 0 & 0 & 0\cr 0 & 1-{2\phi\over c^2} & 0 & 0\cr 0 & 0 & 1-{2\phi\over c^2} & 0\cr 0 & 0 & 0 & 1-{2\phi\over c^2}\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{-c^2-2\phi & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{c^2-2\phi & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{c^2-2\phi & 0 & 0 & 0\cr 0 & c^2 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{c^2-2\phi & 0 & 0 & 0\cr 0 & 1 & 0 & {2\phi\over c^2}\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{c^2-2\phi & 0 & 0 & 0\cr 0 & 1 & 0 & c^2\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^1_{00}=-\infty


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{c^2-2\phi & 0 & 0 & 0\cr 0 & 1 & 0 & 5\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi-5\partial_z\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi


.. math::

    g_{\mu\nu} = \mat{c^2-2\phi & 0 & 5 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi

If we do the limit $c\to\infty$ in the metrics itself, all the working metrics
degenerate to:

.. math::

    g_{\mu\nu} = \mat{\pm\infty & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

(possibly with nonzero but finite elements $g_{0i}=g_{i0}\neq0$).
So it seems like any metrics whose limit is
$\diag(\pm\infty, 1, 1, 1)$, generates the correct Christoffel symbols:

.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi

but this would have to be investigated further.

Let's take the metrics
$\diag(-c^2-2\phi, 1-{2\phi\over c^2}, 1-{2\phi\over c^2}, 1-{2\phi\over c^2})$
and calculate the Christoffel symbols (without the limit $c\to\infty$):

.. math::

    \Gamma^0_{\mu\nu}=\begin{pmatrix}- \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{- 2 \phi\left(t,x,y,z\right) - {c}^{2}} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{- 2 \phi\left(t,x,y,z\right) - {c}^{2}} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{- 2 \phi\left(t,x,y,z\right) - {c}^{2}} & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{- 2 \phi\left(t,x,y,z\right) - {c}^{2}}\\- \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{- 2 \phi\left(t,x,y,z\right) - {c}^{2}} & \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(- 2 \phi\left(t,x,y,z\right) - {c}^{2}\right)} & 0 & 0\\- \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{- 2 \phi\left(t,x,y,z\right) - {c}^{2}} & 0 & \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(- 2 \phi\left(t,x,y,z\right) - {c}^{2}\right)} & 0\\- \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{- 2 \phi\left(t,x,y,z\right) - {c}^{2}} & 0 & 0 & \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(- 2 \phi\left(t,x,y,z\right) - {c}^{2}\right)}\end{pmatrix}

    \Gamma^1_{\mu\nu}=\begin{pmatrix}\frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}} & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & 0 & 0\\- \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)}\\0 & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & 0\\0 & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & 0 & \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)}\end{pmatrix}

    \Gamma^2_{\mu\nu}=\begin{pmatrix}\frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}} & 0 & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & 0\\0 & \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & 0\\- \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)}\\0 & 0 & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)}\end{pmatrix}

    \Gamma^3_{\mu\nu}=\begin{pmatrix}\frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}} & 0 & 0 & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)}\\0 & \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & 0 & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)}\\0 & 0 & \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)}\\- \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)} & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{{c}^{2} \left(1 - 2 \frac{\phi\left(t,x,y,z\right)}{{c}^{2}}\right)}\end{pmatrix}

By taking the limit $c\to\infty$, the only nonzero Christoffel symbols are:

.. math::

    \Gamma^1_{00}=\partial_x\phi


.. math::

    \Gamma^2_{00}=\partial_y\phi


.. math::

    \Gamma^3_{00}=\partial_z\phi

or written compactly:

.. math::

    \Gamma^i_{00}=\delta^{ij}\partial_j\phi

So the geodesics equation

.. math::

    {\d^2 x^\alpha\over\d\lambda^2} + \Gamma^\alpha_{\beta\gamma} {\d x^\beta\over\d\lambda}{\d x^\gamma\over\d\lambda} = 0

becomes

.. math::

    {\d^2 x^0\over\d\lambda^2}=0


.. math::

    {\d^2 x^i\over\d\lambda^2} + \delta^{ij}\partial_j\phi \left({\d x^0\over\d\lambda}\right)^2 = 0

From the first equation we get $x^0 = a\lambda+b$, we substitute to the second
equation:

.. math::

    {1\over a^2}{\d^2 x^i\over\d\lambda^2} + \delta^{ij}\partial_j\phi = 0

or

.. math::

    {\d^2 x^i\over\d (x^0)^2} + \delta^{ij}\partial_j\phi = 0


.. math::

    {\d^2 x^i\over\d t^2}=-\delta^{ij}\partial_j\phi

So the Newton's second law *is* the equation of geodesics.

In the above, we
have set $c=1$ in the Christoffel symbols themselves (see the last paragraph
from the last section) and introduced another constant $c$ in the metric
itself. As we can see, the metric will become infinite with this approach in
the limit $c\to\infty$. Another approach is to store this $c$ in the $x^\mu$
vector itself, then the metric stays finite (in fact becomes a diagonal matrix
$\diag(\pm 1, 1, 1, 1)$, thus it gives all the Christoffel symbols equal to
zero, in the limit), but the vector becomes infinite in the limit.

Either way our formalism breaks down, and thus we need to keep $c$ finite and
only do the limit in the final equations (after we don't need differential
geometry anymore). When needed, we can also carefully neglect higher terms in
$c$, that will not appear in the final equations after doing the limit, but one
needs to make sure that no mistake is made.

It is customary to put the constant
$c$ into the vector $x^\mu$ and so we will do so too from this point on.

Conclusion About Metric
-----------------------

We will use the convention to keep $c$ in the 4-vector and the simplest metric
that generates the correct Christoffel symbols is the following:

.. math::

    g_{\mu\nu} = \mat{\pm 1 -{2\phi\over c^2} & 0 & 0 & 0\cr 0 & 1-{2\phi\over c^2} & 0 & 0\cr 0 & 0 & 1-{2\phi\over c^2} & 0\cr 0 & 0 & 0 & 1-{2\phi\over c^2}\cr}

In the limit $c\to\infty$ we get the following nonzero Christoffel symbols (for
both signs in $\pm 1$ above):

.. math::

    \Gamma^i_{00} = {1\over c^2}\delta^{ij}\partial_j\phi

all other symbols contain higher powers of $c$ and thus will not contribute in
the limit $c\to\infty$. The remaining $c^2$ in $\Gamma^i_{00}$ will cancel with
the $c$ in $x^0=ct$ in the final equations.

As seen above, there is some freedom in which metric we can use in order to
obtain the correct Christoffel symbols, but the above metric is the simplest,
so we'll use it from now on.

Einstein's Equations
--------------------

Einstein's equations are derived from the Hilbert action:

.. math::

    S_H = {c^4\over 16\pi G} \int R \sqrt{ |\det g_{\mu\nu}| } \d^4 x
        = {c^4\over 16\pi G} \int g^{\mu\nu} R_{\mu\nu} \sqrt{ |\det g| } \d^4 x

The Lagrangian density $R \sqrt{ |\det g_{\mu\nu}| }$ has to be given, that's
our assumption and everything else is derived from it. In principle it can have
other terms, for example
$\alpha_1 R^2 + \alpha_2 R_{\mu\nu} R^{\mu\nu} + \alpha_3 g^{\mu\nu} \nabla_\mu R \nabla_\nu R + \cdots$
and there are a lot of possibilities and ultimately the exact form of the
Lagrangian has to be decided by experiment.
The Hilbert action is the simplest possible action and it already gives a
theory which agrees with experiment, so that will be our starting point.

Varying it with respect to the metric $g^{\mu\nu}$ we get:

.. math::

    \delta S_H = \delta {c^4\over 16\pi G} \int R \sqrt{ |\det g| } \d^4 x =

        = {c^4\over 16\pi G} \int
            (\delta g^{\mu\nu}) R_{\mu\nu} \sqrt{ |\det g| }
            +g^{\mu\nu} (\delta R_{\mu\nu}) \sqrt{ |\det g| }
            +R (\delta \sqrt{ |\det g| })
            \d^4 x=

        = {c^4\over 16\pi G} \int
            (\delta g^{\mu\nu}) R_{\mu\nu} \sqrt{ |\det g| }
            +g^{\mu\nu} \left(
                \nabla_\lambda(\delta \Gamma^\lambda_{\nu\mu})
                -\nabla_\nu(\delta \Gamma^\lambda_{\lambda\mu})
                \right)\sqrt{ |\det g| }
            +R (
            -\half \sqrt{ |\det g| }\, g_{\mu\nu} (\delta g^{\mu\nu}))
            \d^4 x=

        = {c^4\over 16\pi G} \int
            (\delta g^{\mu\nu}) R_{\mu\nu} \sqrt{ |\det g| }
            + \left(
                \nabla_\lambda g^{\mu\nu}(\delta \Gamma^\lambda_{\nu\mu})
                -\nabla_\nu g^{\mu\nu}(\delta \Gamma^\lambda_{\lambda\mu})
                \right)\sqrt{ |\det g| }
            -\half R g_{\mu\nu} \sqrt{ |\det g| }\,
                (\delta g^{\mu\nu})
            \d^4 x=

        = {c^4\over 16\pi G} \int
            (\delta g^{\mu\nu}) R_{\mu\nu} \sqrt{ |\det g| }
            -\half R g_{\mu\nu} \sqrt{ |\det g| }\,
                (\delta g^{\mu\nu})
            \d^4 x=

        = {c^4\over 16\pi G} \int \left( R_{\mu\nu} -\half R g_{\mu\nu} \right)
                (\delta g^{\mu\nu}) \sqrt{ |\det g| }
            \d^4 x

Where we used the following identities:

.. math::

    \delta \sqrt{ |\det g| } =
            -\half \sqrt{ |\det g| }\, g_{\mu\nu} (\delta g^{\mu\nu})

    \delta R^\rho{}_{\mu\lambda\nu} =
                \nabla_\lambda(\delta \Gamma^\rho_{\nu\mu})
                -\nabla_\nu(\delta \Gamma^\rho_{\lambda\mu})

    \delta R_{\mu\nu} = \delta R^\lambda{}_{\mu\lambda\nu} =
                \nabla_\lambda(\delta \Gamma^\lambda_{\nu\mu})
                -\nabla_\nu(\delta \Gamma^\lambda_{\lambda\mu})

and the fact that the four divergence doesn't contribute to the
integral. By setting $\delta S_H=0$, we get:

.. math::

    {2\over\sqrt{ |\det g| }}{\delta S_H\over\delta g^{\mu\nu}}
        = {c^4\over 8\pi G}(R_{\mu\nu} -\half R g_{\mu\nu}) = 0

Combining the Hilbert action $S_H$ with the action for matter $S_M$ we get:

.. math::

    S = S_H + S_M

Varying this action as above we get:

.. math::

    {2\over\sqrt{ |\det g| }}{\delta S\over\delta g^{\mu\nu}}
        ={c^4\over 8\pi G} \left( R_{\mu\nu} -\half R g_{\mu\nu} \right)
        +
    {2\over\sqrt{ |\det g| }}{\delta S_M\over\delta g^{\mu\nu}} = 0

so:

.. math::

    R_{\mu\nu} -\half R g_{\mu\nu} =
        -{8\pi G \over c^4}
        {2\over\sqrt{ |\det g| }}{\delta S_M\over\delta g^{\mu\nu}}
    ={8\pi G \over c^4} T_{\mu\nu}

Where we set:

.. math::

    T_{\mu\nu} = - {2\over\sqrt{ |\det g| }}{\delta S_M\over\delta g^{\mu\nu}}

This is a definition of the stress energy tensor corresponding to the action
$S_M=\int \L_M \sqrt{ |\det g| }\d^4 x$. We can also write it in terms of the
Lagrangian $\L_M$ directly as:

.. math::

    T_{\mu\nu} = - {2\over\sqrt{ |\det g| }}{\delta S_M\over\delta g^{\mu\nu}}=

        = - {2\over\sqrt{ |\det g| }}{\delta \int \L_M \sqrt{ |\det g| }\d^4 x
                \over\delta g^{\mu\nu}} =

        = - {2\over\sqrt{ |\det g| }}{\int (\delta \L_M) \sqrt{ |\det g| }
            + \L_M \left(\delta \sqrt{ |\det g| }\right)\d^4 x
                \over\delta g^{\mu\nu}} =

        = - {2\over\sqrt{ |\det g| }}{\int \left({\delta \L_M
            \over\delta g^{\mu\nu}}(\delta g^{\mu\nu})\right) \sqrt{ |\det g| }
            + \L_M \left(-\half
                \sqrt{ |\det g| }\, g_{\mu\nu} (\delta g^{\mu\nu})\right)\d^4 x
                \over\delta g^{\mu\nu}} =

        = - {2\over\sqrt{ |\det g| }}\left({\delta \L_M
            \over\delta g^{\mu\nu}} \sqrt{ |\det g| }
            - \half \L_M \sqrt{ |\det g| }\, g_{\mu\nu} \right) =

        = - 2 {\delta \L_M \over\delta g^{\mu\nu}}
            + g_{\mu\nu} \L_M

If this action contains electromagnetic field, we get an electromagnetic
stress energy tensor. For continous matter, we get the stress energy tensor for
continous matter, see the next section. The right hand side of the Einstein's
equations contains the sum of all stress energy tensors (for all fields in the
Lagrangian).


Continuous Distribution of Matter
---------------------------------

The action is:

.. math::

    S_M = -\int \rho c \sqrt{v_\mu v^\mu} \sqrt{ |\det g| } \d^4 x

But it isn't suitable for applying variations because $\rho$ and $v^\mu$ are
not independent quantities. So we write it in terms of a 4-momentum vector
density $\textgoth{p}^\mu$:

.. math::

    p^\mu = \rho v^\mu

    \textgoth{p}^\mu = p^\mu \sqrt{ |\det g| } = \rho v^\mu \sqrt{ |\det g| }

    \sqrt{\textgoth{p}_\mu \textgoth{p}^\mu}
        = \sqrt{\rho v_\mu \sqrt{ |\det g| } \rho v^\mu \sqrt{ |\det g| }}
        = \rho \sqrt{v_\mu v^\mu} \sqrt{ |\det g| }

and the action becomes:

.. math::

    S_M = -\int \rho c \sqrt{v_\mu v^\mu} \sqrt{ |\det g| } \d^4 x
        = -\int c \sqrt{\textgoth{p}_\mu \textgoth{p}^\mu} \d^4 x

We vary $S_M$ with respect to $g^{\mu\nu}$:

.. math::

    \delta S_M
        = - \delta \int c \sqrt{\textgoth{p}_\mu \textgoth{p}^\mu} \d^4 x =

        = - \int c {\delta(g^{\mu\nu} \textgoth{p}_\mu \textgoth{p}_\nu)
            \over 2\sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}} \d^4 x =

        = - \int c { \textgoth{p}_\mu \textgoth{p}_\nu
            \over 2\sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}}
            \delta(g^{\mu\nu})\d^4 x =

        = - \int c { \rho v_\mu \rho v_\nu
            \sqrt{ |\det g| }^2
            \over 2 \rho c \sqrt{ |\det g| } }
             \delta(g^{\mu\nu})\d^4 x =

        = - \int \half \rho v_\mu v_\nu
             \delta(g^{\mu\nu}) \sqrt{ |\det g| } \d^4 x

And the stress energy tensor is:

.. math::

    T_{\mu\nu}
        = - {2\over\sqrt{ |\det g| }}{\delta S_M\over\delta g^{\mu\nu}} =

        = - {2\over\sqrt{ |\det g| }} \left(
                -\half \rho v_\mu v_\nu \sqrt{ |\det g| }
            \right)=

        = \rho v_\mu v_\nu

Now we vary $S_M$ with respect to $x^\mu$:

.. math::

    \delta S_M
        = - \delta \int c \sqrt{\textgoth{p}_\mu \textgoth{p}^\mu} \d^4 x =

        = - \int c {\delta(g^{\mu\nu} \textgoth{p}_\mu \textgoth{p}_\nu)
            \over 2\sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}} \d^4 x =

        = - \int c { 2 g^{\mu\nu} \textgoth{p}_\mu (\delta \textgoth{p}_\nu)
            \over 2\sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}} \d^4 x =

        = - \int c { \textgoth{p}_\mu \over \sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}}
            (\delta \textgoth{p}^\mu) \d^4 x =

        = - \int c { \textgoth{p}_\mu \over \sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}}
            \partial_\nu \left(\textgoth{p}^\nu(\delta x^\mu) - \textgoth{p}^\mu (\delta x^\nu)\right)
                \d^4 x =

        = \int c \partial_\nu \left({ \textgoth{p}_\mu \over
        \sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}}
                \right)
            \left(\textgoth{p}^\nu(\delta x^\mu) - \textgoth{p}^\mu (\delta x^\nu)\right)
                \d^4 x =

        = \int c \left(
            \partial_\nu \left({ \textgoth{p}_\mu \over \sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}} \right)
            -\partial_\mu \left({ \textgoth{p}_\nu \over \sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}} \right)
            \right)
            \textgoth{p}^\nu(\delta x^\mu)
                \d^4 x =

        = \int c \left(
            \nabla_\nu \left({ \textgoth{p}_\mu \over \sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}} \right)
            -\nabla_\mu \left({ \textgoth{p}_\nu \over \sqrt{\textgoth{p}_\alpha \textgoth{p}^\alpha}} \right)
            \right)
            \textgoth{p}^\nu(\delta x^\mu)
                \d^4 x =

        = \int \left( \nabla_\nu v_\mu -\nabla_\mu v_\nu \right)
            \rho v^\nu (\delta x^\mu) \sqrt{ |\det g| }
                \d^4 x =

        = \int \rho (\nabla_\nu v_\mu) v^\nu (\delta x^\mu) \sqrt{ |\det g| }
                \d^4 x

So the equation of motion is the geodesic:

.. math::

    \rho (\nabla_\nu v_\mu) v^\nu = 0

Charged matter has the interaction action:

.. math::

    S_q = -\int \rho_{EM} v^\mu A_\mu \sqrt{ |\det g| } \d^4 x
        = -\int j^\mu A_\mu \sqrt{ |\det g| }\d^4 x
        = -\int \textgoth{j}^\mu A_\mu \d^4 x

where we have introduced the 4-current $j^\mu$ and 4-current density
$\textgoth{j}^\mu$:

.. math::

    j^\mu = \rho_{EM} v^\mu

    \textgoth{j}^\mu = j^\mu \sqrt{ |\det g| }
        = \rho_{EM} v^\mu \sqrt{ |\det g| }

We vary $S_q$ with respect to $x^\mu$:

.. math::

    \delta S_q
        = - \delta \int \textgoth{j}^\mu A_\mu \d^4 x =

        = - \int (\delta \textgoth{j}^\mu) A_\mu \d^4 x =

        = - \int \partial_\nu \left(\textgoth{j}^\nu (\delta x^\mu)
            - \textgoth{j}^\mu (\delta x^\nu)\right) A_\mu \d^4 x =

        = \int \left(\textgoth{j}^\nu (\delta x^\mu)
            - \textgoth{j}^\mu (\delta x^\nu)\right) \partial_\nu A_\mu \d^4 x =

        = \int \textgoth{j}^\nu (\delta x^\mu) (\partial_\nu A_\mu -\partial_\mu A_\nu)
            \d^4 x =

        = \int \rho_{EM} v^\nu  (\partial_\nu A_\mu -\partial_\mu A_\nu)
            (\delta x^\mu) \sqrt{ |\det g| }
            \d^4 x =

        = -\int \rho_{EM} v^\nu  F_{\mu\nu} (\delta x^\mu) \sqrt{ |\det g| }
            \d^4 x

So the combined action $S_M + S_q$ yields:

.. math::

    \rho (\nabla_\nu v_\mu) v^\nu
        -\rho_{EM} v^\nu  F_{\mu\nu} = 0

Varying $S_q$ with respect to $A^\mu$ yields the 4-current $j^\mu = \rho_{EM}
v^\mu$ which ends up on the right hand side of the Maxwell's equations when
varying the $S_{EM}$ action.

Obsolete Section
----------------


This section is obsolete, ideas from it should be polished (sometimes
corrected) and put to other sections.

The problem is, that in general, Christoffel symbols have 40 components and
metrics only 10 and in our case, we cannot find such a metrics, that generates
the Christoffel symbols above. In other words, the spacetime that describes
the Newtonian theory is affine, but not a metric space. The metrics is singular,
and we have one metrics $\diag(-1, 0, 0, 0)$ that describes the time coordinate
and another metrics $\diag(0, 1, 1, 1)$ that describes the spatial coordinates.
We know the affine connection coefficients $\Gamma^\alpha_{\beta\gamma}$, so
that is enough to calculate geodesics and to differentiate vectors and do
everything we need.

However, for me it is still not satisfactory, because I really want to have a
metrics tensor, so that I can easily derive things in exactly the same way as
in general relativity. To do that, we will have to work in the regime $c$ is
finite and only at the end do the limit $c\to\infty$.

We start with Einstein's equations:

.. math::

    R_{\alpha\beta}-\half Rg_{\alpha\beta}={8\pi G\over c^4}T_{\alpha\beta}

or

.. math::

    R_{\alpha\beta}={8\pi G\over c^4}(T_{\alpha\beta}-\half Tg_{\alpha\beta})


.. math::

    R^\alpha{}_\beta={8\pi G\over c^4}(T^\alpha{}_\beta-\half T)

The energy-momentum tensor is

.. math::

    T^{\alpha\beta} = \rho U^\alpha U^\beta

in our approximation $U^i \sim0$ and $U^0 \sim c$, so the only nonzero component
is:

.. math::

    T^{00} = \rho c^2


.. math::

    T = \rho c^2

and

.. math::

    R^i{}_j={8\pi G\over c^4}(-\half \rho c^2)=-{4\pi G\over c^2}\rho


.. math::

    R^0{}_0={8\pi G\over c^4}(\half \rho c^2)={4\pi G\over c^2}\rho

We need to find such a metric tensor, that

.. math::

    R^0{}_0={1\over c^2}\nabla^2\phi

then we get :eq:`grav`.

There are several ways to choose the metrics tensor. We
start
We can always find a coordinate transformation, that converts the metrics to a
diagonal form with only $1$, $0$ and $-1$ on the diagonal. If we want
nondegenerate metrics, we do not accept $0$ (but as it turns out, the metrics
for the Newtonian mechanics *is* degenerated).
Also, it is equivalent if we add a minus to all diagonal elements, e.g. $\diag(1,
1, 1, 1)$ and $\diag(-1, -1, -1, -1)$ are equivalent, so
we are left
with these options only:
signature 4:

.. math::

    g_{\mu\nu}=\diag(1, 1, 1, 1)

signature 2:

.. math::

    g_{\mu\nu}=\diag(-1, 1, 1, 1)


.. math::

    g_{\mu\nu}=\diag(1, -1, 1, 1)


.. math::

    g_{\mu\nu}=\diag(1, 1, -1, 1)


.. math::

    g_{\mu\nu}=\diag(1, 1, 1, -1)

signature 0:

.. math::

    g_{\mu\nu}=\diag(-1, -1, 1, 1)


.. math::

    g_{\mu\nu}=\diag(-1, 1, -1, 1)


.. math::

    g_{\mu\nu}=\diag(-1, 1, 1, -1)

No other possibility exists (up to adding a minus to all elements). We can also
quite easily find coordinate transformations that swap coordinates, i.e. we can
always find a transformation so that we first have only $-1$ and then only $1$
on the diagonal, so we are left with:
signature 4:

.. math::

    g_{\mu\nu}=\diag(1, 1, 1, 1)

signature 2:

.. math::

    g_{\mu\nu}=\diag(-1, 1, 1, 1)

signature 0:

.. math::

    g_{\mu\nu}=\diag(-1, -1, 1, 1)

One possible physical interpretation of the signature 0 metrics is
that we have 2 time coordinates and 2 spatial coordinates. In any case, this
metrics doesn't describe our space (neither Newtonian nor general relativity),
because we really need the spatial coordinates to have the metrics either
$\diag(1, 1, 1)$ or $\diag(-1, -1, -1)$.

So we are left with either (this case will probably not work, but I want to
have an
explicit reason why it doesn't work):

.. math::

    g_{\mu\nu} = \mat{1 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

or (this is the usual special relativity)

.. math::

    g_{\mu\nu} = \mat{-1 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

It turns out, that one option to turn on gravitation is to add the term $-{2\phi\over c^2}\one$ to the
metric tensor, in the first
case:

.. math::

    g_{\mu\nu} = \mat{1-{2\phi\over c^2} & 0 & 0 & 0\cr 0 & 1-{2\phi\over c^2} & 0 & 0\cr 0 & 0 & 1-{2\phi\over c^2} & 0\cr 0 & 0 & 0 & 1-{2\phi\over c^2}\cr}

and second case:

.. math::

    g_{\mu\nu} = \mat{-1-{2\phi\over c^2} & 0 & 0 & 0\cr 0 & 1-{2\phi\over c^2} & 0 & 0\cr 0 & 0 & 1-{2\phi\over c^2} & 0\cr 0 & 0 & 0 & 1-{2\phi\over c^2}\cr}

The second law is derived from the
equation of geodesic:

.. math::

    {\d^2 x^\alpha\over\d\lambda^2} + \Gamma^\alpha_{\beta\gamma} {\d x^\beta\over\d\lambda}{\d x^\gamma\over\d\lambda} = 0

in an equivalent form

.. math::

    {\d U^\alpha\over\d\tau} + \Gamma^\alpha_{\beta\gamma}U^\beta U^\gamma = 0

The only nonzero Christoffel symbols in the first case are (in the expressions
for the Christoffel symbols below, we set $c=1$):

.. math::

    \Gamma^0_{\mu\nu}= \begin{pmatrix}- \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\\- \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0 & 0\\- \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0 & \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0\\- \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0 & 0 & \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\end{pmatrix}

    \Gamma^1_{\mu\nu}= \begin{pmatrix}\frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0 & 0\\- \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\\0 & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0\\0 & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0 & \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\end{pmatrix}

    \Gamma^2_{\mu\nu}= \begin{pmatrix}\frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0 & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0\\0 & \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0\\- \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\\0 & 0 & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\end{pmatrix}

    \Gamma^3_{\mu\nu}= \begin{pmatrix}\frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0 & 0 & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\\0 & \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & 0 & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\\0 & 0 & \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\\- \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 - 2 \phi\left(t,x,y,z\right)}\end{pmatrix}

and in the second case, only $\Gamma^0_{\mu\nu}$ is different:

.. math::

    \Gamma^0_{\mu\nu}= \begin{pmatrix}\frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)} & \frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)} & \frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)} & \frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)}\\\frac{\frac{\partial}{\partial x} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)} & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)} & 0 & 0\\\frac{\frac{\partial}{\partial y} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)} & 0 & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)} & 0\\\frac{\frac{\partial}{\partial z} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)} & 0 & 0 & - \frac{\frac{\partial}{\partial t} \phi\left(t,x,y,z\right)}{1 + 2 \phi\left(t,x,y,z\right)}\end{pmatrix}


Now we assume that $\partial_\mu\phi \sim \phi \ll c^2$, so all $\Gamma^\alpha_{\beta
\gamma}$ are of the same order. Also $|U^i| \ll |U^0|$ and $U^0 = c$, so the only
nonnegligible term is

.. math::

    {\d U^\alpha\over\d\tau} + \Gamma^\alpha_{00}(U^0)^2 = 0

Substituting for the Christoffel symbol we get

.. math::

    {\d U^i\over\d\tau} =-{\delta^{ij}\partial_j{\phi\over c^2}\over1-{2\phi\over c^2}} \, c^2 =-\delta^{ij}(\partial_j\phi)\ \left(1+O\left({\phi\over c^2}\right)\right) =-\delta^{ij}\partial_j\phi  + O\left(\left({\phi\over c^2}\right)^2\right)

and multiplying both sides with $m$:

.. math::

    m{\d U^i\over\d\tau} =-m\partial_j\phi\ \delta^{ij}

which is the second Newton's law. For the zeroth component we get (first case
metric)

.. math::

    m{\d U^0\over\d\tau} =m{\d\phi\over\d\tau}

second case:

.. math::

    m{\d U^0\over\d\tau} =-m{\d\phi\over\d\tau}

Where $mU^0 = p^0$ is the energy of the particle (with respect to this frame
only), this means the energy is conserved unless the gravitational field
depends on time.

To summarize: the Christoffel symbols :eq:`Chris-newton` that we get from the
Newtonian theory contain $c$, which up to this point can be any speed, for
example we can set $c=1\rm\,ms^{-1}$. However, in order to have some metrics
tensor that generates those Christoffel symbols, the only way to do that is by
the metrics

.. math::

    \diag(-1, 1, 1, 1)-{2\phi\over c^2}\one

then calculating the Christoffel symbols. If we neglect the terms of the order
$O\left(\left(\phi\over c^2\right)^2\right)$ and higher, we get the Newtonian
Christoffel symbols :eq:`Chris-newton` that we want. It's clear that in order
to neglect the terms, we must have $|\phi| \ll c^2$, so we must choose $c$
large enough for this to work. To put it plainly, unless $c$ is large, there is
no metrics in our Newtonian spacetime. However for $c$ large, everything is
fine.


Inertial frames
---------------


What is an inertial frame? Inertial frame is such a frame
that doesn't have any fictitious forces. What is a fictitious force?
If we take covariant time derivative of any vector, then fictitious forces are all
the terms with nonzero Christoffel symbols. In other words, nonzero Christoffel symbols
mean that by (partially) differentiating with respect to time, we need to add
additional terms in order to get a proper vector again -- and those terms are
called fictitious forces if we are differentiating the velocity vector.

Inertial frame is a frame without fictitious forces, i.e. with all Christoffel
symbols zero in the whole frame.  This is equivalent to all components of the
Riemann tensor being zero:

.. math::

    R^\alpha{}_{\beta\gamma\delta} = 0

In general, if $R^\alpha{}_{\beta\gamma\delta} \neq 0$ in the whole universe,
then no such frame exists, but one can always achieve that locally, because
one can always find a coordinate transformation so that the Christoffel
symbols are zero locally (e.g. at one point), but unless
$R^\alpha{}_{\beta\gamma\delta} = 0$, the Christoffel symbols will *not*
be zero in the whole frame. So the (local) inertial frame is such a frame that
has zero Christoffel symbols (locally).

What is the metrics of the inertial frame? It is such a metrics, that
$\Gamma^\alpha{}_{\beta\gamma} = 0$. The derivatives
$\partial_\mu\Gamma^\alpha{}_{\beta\gamma}$ however doesn't have to be zero. We
know that taking any of the metrics listed above with $\phi=const$ we get all
the Christoffel symbols zero. So for example these two metrics (one with a plus
sign, the other with a minus sign) have all the Christoffel symbols zero:

.. math::

    g_{\mu\nu} = \mat{\pm c^2 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

Such a metrics corresponds to an inertial frame then.

What are the (coordinate) transformations, that transform from one
inertial frame to another? Those are all transformations that start with an
inertial frame metrics (an example of such a metrics is given above), transform
it using the transformation matrix and the resulting metrics is also inertial.
In particular, let $x^\mu$ be inertial, thus $g_{\mu\nu}$ is an inertial
metrics, then transform to $x'^\mu$ and $g'$:

.. math::

    g'_{\alpha\beta} = {\partial x^\mu\over\partial x'^\alpha} {\partial x^\nu\over\partial x'^\beta} g_{\mu\nu} = \left({\partial x\over\partial x'}\right)^T g \left({\partial x\over\partial x'}\right)

if we denote the transformation matrix by $\Lambda$:

.. math::

    \Lambda^\mu{}_\alpha= {\partial x^\mu\over\partial x'^\alpha}

then the transformation law is:

.. math::

     g' = \Lambda^T g \Lambda

Now let's assume that $g'=g$, i.e. both inertial systems are given by the same
matrix and let's assume this particular form:

.. math::

    g'_{\mu\nu}=g_{\mu\nu} = \mat{\pm c^2 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

(e.g. this covers almost all possible Newtonian metrics tensors).

.. index::
    pair: Lorentz; Group

Lorentz Group
-------------


The Lorentz group is O(3,1), e.g. all matrices satisfying:

.. math::
    :label: ortho

    g =
    \Lambda^T
    g
    \Lambda

with $g=\diag(-c^2, 1, 1, 1)$.
Taking the determinant of :eq:`ortho` we get $(\det\Lambda)^2=1$ or
$\det\Lambda=\pm1$. Writing the 00 component of :eq:`ortho` we get

.. math::

     -c^2 = -c^2(A^0{}_0)^2+(A^0{}_1)^2+(A^0{}_2)^2+(A^0{}_3)^2

or

.. math::

     (A^0{}_0)^2 = 1 + {1\over c^2}\left((A^0{}_1)^2+(A^0{}_2)^2+(A^0{}_3)^2\right)

Thus we can see that either $A^0{}_0\ge1$ (the transformation preserves the
direction of time, orthochronous) or $A^0{}_0\le-1$ (not orthochronous).
Thus we can see that the O(3, 1) group consists of 4 continuous parts, that
are not connected.

First case: elements with $\det\Lambda=1$ and $A^0{}_0\ge1$. Transformations
with $\det\Lambda=1$ form a subgroup and are called SO(3, 1), if they also have
$A^0{}_0\ge1$ (orthochronous), then they also form a subgroup and are called
the proper Lorentz transformations and denoted by ${\rm SO}^+(3, 1)$. They
consists of Lorentz boosts, example in the $x$-direction:

.. math::

    \Lambda^\mu{}_\nu= \mat{ {1\over\sqrt{1-{v^2\over c^2}}}& -{{v\over c^2}\over\sqrt{1-{v^2\over c^2}}} & 0 & 0\cr -{v\over\sqrt{1-{v^2\over c^2}}} & {1\over\sqrt{1-{v^2\over c^2}}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

which in the limit $c\to\infty$ gives

.. math::

    \Lambda^\mu{}_\nu= \mat{ 1 & 0 & 0 & 0\cr -v & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

and spatial rotations:

.. math::

    R_1(\phi)= \mat{ 1 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & \cos\phi & \sin\phi\cr 0 & 0 & -\sin\phi & \cos\phi\cr}


.. math::

    R_2(\phi)= \mat{ 1 & 0 & 0 & 0\cr 0 & \cos\phi & 0 & \sin\phi\cr 0 & 0 & 1 & 0\cr 0 & -\sin\phi & 0  & \cos\phi\cr}


.. math::

    R_3(\phi)= \mat{ 1 & 0 & 0 & 0\cr 0 & \cos\phi & \sin\phi & 0\cr 0 & -\sin\phi & \cos\phi & 0\cr 0 & 0 & 0 & 1\cr}

(More rigorous derivation will be given in a moment.)
It can be shown (see below), that all other elements (improper
Lorentz transformations) of the O(3, 1)
group can be written as products of an element from ${\rm SO}^+(3, 1)$ and an
element of the discrete group:

.. math::

    \{\one,\ P,\ T,\ PT\}

where $P$ is parity (also called space reflection or space inversion):

.. math::

    P= \mat{ 1 & 0 & 0 & 0\cr 0 & -1 & 0 & 0\cr 0 & 0 & -1 & 0\cr 0 & 0 & 0 & -1\cr}

and $T$ is time reversal (also called time inversion):

.. math::

    T= \mat{ -1 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}


Second case: elements with $\det\Lambda=1$ and $A^0{}_0\le-1$. An example of
such an element is $PT$. In general, any product from ${\rm SO}^+(3, 1)$ and
$PT$ belongs here.

Third case: elements with $\det\Lambda=-1$ and $A^0{}_0\ge1$. An example of
such an element is $P$. In general, any product from ${\rm SO}^+(3, 1)$ and
$P$ belongs here.

Fourth case: elements with $\det\Lambda=-1$ and $A^0{}_0\le-1$. An example of
such an element is $T$. In general, any product from ${\rm SO}^+(3, 1)$ and
$T$ belongs here.

Example: where does the reflection around a
single spatial axis $(t, x, y, z)\to(t, -x, y, z)$ belong to? It is the third
case, because the determinant is $\det\Lambda=-1$ and the 00 element is 1.
Written in the matrix form:

.. math::

    \Lambda =
     \mat{ 1 & 0 & 0 & 0\cr 0 & -1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 &
     1\cr} =
     \mat{ 1 & 0 & 0 & 0\cr 0 & -1 & 0 & 0\cr 0 & 0 & -1 & 0\cr 0 & 0 &
     0 & -1\cr}
     \mat{ 1 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & -1 &
     0\cr 0 & 0 & 0 & -1\cr}
     =

     =
     \mat{ 1 & 0 & 0 & 0\cr 0 & -1 & 0 & 0\cr 0 & 0 & -1 & 0\cr 0 & 0 &
     0 & -1\cr}
     \mat{ 1 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr 0 & 0 & \cos\pi &
     \sin\pi\cr 0 & 0 & -\sin\pi & \cos\pi\cr}
     =PR_1(\pi)

So it is constructed using the $R_1$ element from ${\rm SO}^+(3, 1)$ and P
from the discrete group above.

We can now show why the decomposition
${\rm O}(3,1)={\rm SO}^+(3, 1)\times\{\one,\ P,\ T,\ PT\}$ works. Note that
$PT=-\one$. First we show that ${\rm SO}(3,1)={\rm SO}^+(3, 1)\times\{\one,
-\one\}$. This follows from the fact, that all matrices with
$\Lambda^0{}_0\le-1$ can be written using $-\one$ and a matrix with
$\Lambda^0{}_0\ge1$.  All matrices
with $\det \Lambda=-1$ can be constructed from a matrix with $\det\Lambda=1$
(i.e. SO(3, 1)) and
a diagonal matrix with odd number of -1, below we list all of them together
with their construction using time reversal, parity and spatial rotations:

.. math::

    \diag(-1, 1, 1, 1) &= T\\
    \diag(1, -1, 1, 1) &= PR_1(\pi)\\
    \diag(1, 1, -1, 1) &= PR_2(\pi)\\
    \diag(1, 1, 1, -1) &= PR_3(\pi)\\
    \diag(1, -1, -1, -1) &= P\\
    \diag(-1, 1, -1, -1) &= TR_1(\pi)\\
    \diag(-1, -1, 1, -1) &= TR_2(\pi)\\
    \diag(-1, -1, -1, 1) &= TR_3(\pi)\\

But $R_i(\pi)$ belongs to ${\rm SO}^+(3, 1)$, so we just need two extra
elements, $T$ and $P$ to construct all matrices with $\det\Lambda=-1$ using
matrices from SO(3, 1). So to recapitulate, if we start with ${\rm SO}^+(3, 1)$
we need to add the element $PT=-\one$ to construct SO(3, 1) and then we need to
add $P$ and $T$ to construct O(3, 1). Because all other combinations like
$PPT=T$ reduce to just one of $\{\one, P, T, -\one\}$, we are done.

The elements from ${\rm SO}^+(3, 1)$ are proper Lorentz transformations, all
other elements are improper. Now we'd like to construct the proper
Lorentz transformation matrix $A$ explicitly. As said above, all improper
transformations are just proper transformations multiplied by either $P$, $T$
or $PT$, so it is sufficient to construct $A$.

We can always write $A=e^L$, then:

.. math::

    \det A = \det e^L = e^{{\rm Tr}\,L} = 1

so $\Tr L = 0$ and $L$ is a real, traceless matrix. Rewriting :eq:`ortho`:

.. math::

     g = A^T g A


.. math::

     A^{-1} = g^{-1} A^T g


.. math::

     e^{-L} = g^{-1} e^{L^T}g = e^{g^{-1}L^Tg}


.. math::

     -L = g^{-1}L^Tg


.. math::

     -gL = (gL)^T

The matrix $gL$ is thus antisymmetric and the general form of $L$ is then:

.. math::

     L= \mat{ 0 & {L_{01}\over c^2} & {L_{02}\over c^2} & {L_{03}\over c^2}\cr L_{01} & 0 & L_{12} & L_{13}\cr L_{02} & -L_{12} & 0 & L_{23}\cr L_{03} & -L_{13} & -L_{23} & 0\cr}

One can check, that $gL$ is indeed antisymmetric. However, for a better
parametrization, it's better to work with a metric $\diag(-1, 1, 1, 1)$, which
can be achieved by putting $c$ into $(ct, x, y, z)$, or equivalently, to work
with $x^\mu=(t, x, y, z)$ and multiply this by a matrix $C=\diag(c, 1, 1, 1)$
to get $(ct, x, y, z)$. To get a symmetric $\tilde L$, we just have to do
$Cx' = \tilde LCx$,
so to get an unsymmetric $L$ from the symmetric one, we need to do
$C^{-1} \tilde L C$, so we get:

.. math::

     L = C^{-1} \mat{ 0 & \zeta_1 & \zeta_2 & \zeta_3\cr \zeta_1 & 0 & -\varphi_3 & \varphi_2\cr \zeta_2 & \varphi_3 & 0 & -\varphi_1\cr \zeta_3 & -\varphi_2 & \varphi_1 & 0\cr} C = -i\boldsymbol\varphi\cdot{\bf L}-i\boldsymbol\zeta\cdot C^{-1}{\bf M}C

We have parametrized all the proper Lorentz transformations with just 6
parameters
$\zeta_1$,
$\zeta_2$,
$\zeta_3$,
$\varphi_1$,
$\varphi_2$ and
$\varphi_3$. The matrices ${\bf L}$ and ${\bf M}$ are defined as:

.. math::

     L_1=-i\mat{ 0 & 0 & 0 & 0\cr 0 & 0 & 0 & 0\cr 0 & 0 & 0 & 1\cr 0 & 0 & -1 & 0\cr}


.. math::

     L_2=-i\mat{ 0 & 0 & 0 & 0\cr 0 & 0 & 0 & -1\cr 0 & 0 & 0 & 0\cr 0 & 1 & 0 & 0\cr}


.. math::

     L_3=-i\mat{ 0 & 0 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & -1 & 0 & 0\cr 0 & 0 & 0 & 0\cr}


.. math::

     M_1=i\mat{ 0 & 1 & 0 & 0\cr 1 & 0 & 0 & 0\cr 0 & 0 & 0 & 0\cr 0 & 0 & 0 & 0\cr}


.. math::

     M_2=i\mat{ 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 0\cr 1 & 0 & 0 & 0\cr 0 & 0 & 0 & 0\cr}


.. math::

     M_3=i\mat{ 0 & 0 & 0 & 1\cr 0 & 0 & 0 & 0\cr 0 & 0 & 0 & 0\cr 1 & 0 & 0 & 0\cr}

Straightforward calculation shows:

.. math::

    [L_i, L_j] = i\epsilon_{ijk}L_k


.. math::

    [L_i, M_j] = i\epsilon_{ijk}M_k


.. math::

    [M_i, M_j] = -i\epsilon_{ijk}L_k

The first relation corresponds to the commutation relations for angular
momentum, second relation shows that $M$ transforms as a vector under rotations
and the final relation shows that boosts do not in general commute.

We get:

.. math::

     A = e^{-i\boldsymbol\varphi\cdot{\bf L}-i\boldsymbol\zeta\cdot C^{-1}{\bf M}C} = C^{-1}\,e^{-i\boldsymbol\varphi\cdot{\bf L}-i\boldsymbol\zeta\cdot{\bf M}}\,C

As a special case, the rotation around the $z$-axis is given by
$\boldsymbol\varphi=(0, 0, \varphi)$ and $\boldsymbol\zeta=0$:

.. math::

    A= e^{-i\varphi L_3} = \one-L_3^2+iL_3\sin\varphi+L_3^2\cos\varphi= \mat{ 1 & 0 & 0 & 0\cr 0 & \cos\varphi & \sin\varphi & 0\cr 0 & -\sin\varphi & \cos\varphi & 0\cr 0 & 0 & 0 & 1\cr}

The boost in the $x$-direction is $\boldsymbol\varphi=0$ and
$\boldsymbol\zeta=(\zeta, 0, 0)$, e.g.:

.. math::

    A= C^{-1}e^{-i\zeta M_1}C = C^{-1}\left( \one+M_1^2-iM_1\sinh\zeta- M_1^2\cosh\zeta\right)C=


.. math::

     = C^{-1} \mat{ \cosh\zeta & \sinh\zeta & 0 & 0\cr \sinh\zeta & \cosh\zeta & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr} C = \mat{ \cosh\zeta & {1\over c}\sinh\zeta & 0 & 0\cr c\sinh\zeta & \cosh\zeta & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

from the construction, $-\infty<\zeta<\infty$, so we may do the
substitution $\zeta={\rm atanh}\left({v\over c}\right)$, where
$-c<v<c$. The inverse transformation is:

.. math::

    \cosh\zeta={1\over\sqrt{1-{v^2\over c^2}}}


.. math::

    \sinh\zeta={{v\over c}\over\sqrt{1-{v^2\over c^2}}}

and we get the boost given above:

.. math::

    A= \mat{ \cosh\zeta & {1\over c}\sinh\zeta & 0 & 0\cr c\sinh\zeta & \cosh\zeta & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr} = \mat{ {1\over\sqrt{1-{v^2\over c^2}}}& {{v\over c^2}\over\sqrt{1-{v^2\over c^2}}} & 0 & 0\cr {v\over\sqrt{1-{v^2\over c^2}}} & {1\over\sqrt{1-{v^2\over c^2}}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

Depending on the sign of $v$, we can also put a minus sign in front of the
off-diagonal elements.

Adding two boosts together:

.. math::

     A(u)A(v) = \mat{ {1\over\sqrt{1-{u^2\over c^2}}}& -{{u\over c^2}\over\sqrt{1-{u^2\over c^2}}} & 0 & 0\cr -{u\over\sqrt{1-{u^2\over c^2}}} & {1\over\sqrt{1-{u^2\over c^2}}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr} \mat{ {1\over\sqrt{1-{v^2\over c^2}}}& -{{v\over c^2}\over\sqrt{1-{v^2\over c^2}}} & 0 & 0\cr -{v\over\sqrt{1-{v^2\over c^2}}} & {1\over\sqrt{1-{v^2\over c^2}}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr} =


.. math::

    = \mat{ {1\over\sqrt{1-{w^2\over c^2}}}& -{{w\over c^2}\over\sqrt{1-{w^2\over c^2}}} & 0 & 0\cr -{w\over\sqrt{1-{w^2\over c^2}}} & {1\over\sqrt{1-{w^2\over c^2}}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

with

.. math::

    w = {u+v\over1+{uv\over c^2}}

.. index:: O(4) Group

O(4) Group
----------


The group of rotations in 4 dimensions is O(4), e.g. all matrices satisfying:

.. math::
    :label: ortho2

    g =
    \Lambda^T
    g
    \Lambda

with $g=\diag(c^2, 1, 1, 1)$.
Taking the determinant of :eq:`ortho2` we get $(\det\Lambda)^2=1$ or
$\det\Lambda=\pm1$. Writing the 00 component of :eq:`ortho2` we get

.. math::

     c^2 = c^2(A^0{}_0)^2+(A^0{}_1)^2+(A^0{}_2)^2+(A^0{}_3)^2

or

.. math::

     (A^0{}_0)^2 = 1 - {1\over c^2}\left((A^0{}_1)^2+(A^0{}_2)^2+(A^0{}_3)^2\right)

Thus we always have $-1\le A^0{}_0\le1$. That is different to the O(3, 1)
group: the O(4) group consists of only 2 continuous parts, that are not
connected. (The SO(4) part  contains the element $-\one$ though, but one can get
to it continuously, so the group is doubly connected.)

Everything proceeds much like for the O(3, 1) group, so $gL$ is antisymmetric,
but this time $g=\diag(c^2, 1, 1, 1)$, so we get:

.. math::

     L= \mat{ 0 & -{L_{01}\over c^2} & -{L_{02}\over c^2} & -{L_{03}\over c^2}\cr L_{01} & 0 & L_{12} & L_{13}\cr L_{02} & -L_{12} & 0 & L_{23}\cr L_{03} & -L_{13} & -L_{23} & 0\cr}

and so we also have 6 generators, but this time all of them are rotations:

.. math::

     A = C^{-1}\,e^{-i\varphi_a L_a}\,C

with $a=1, 2, 3, 4, 5, 6$. The spatial rotations are the same as for O(3, 1)
and the remaining 3 rotations are $(t,x)$, $(t,y)$ and $(t,z)$ plane rotations.
So for example the $(t,x)$ rotation is:

.. math::

     A= C^{-1} \mat{ \cos\varphi_4 & \sin\varphi_4 & 0 & 0\cr -\sin\varphi_4 & \cos\varphi_4 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr} C = \mat{ \cos\varphi_4 & {1\over c}\sin\varphi_4 & 0 & 0\cr -c\sin\varphi_4 & \cos\varphi_4 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

Now we can do this identification:

.. math::

    \sin\phi_4 = {{v\over c}\over\sqrt{1+({v\over c})^2}}


.. math::

    \cos\phi_4 = {1\over\sqrt{1+({v\over c})^2}}

so we get the Galilean transformation in the limit $c\to\infty$:

.. math::

    A= \mat{ {1\over\sqrt{1+({v\over c})^2}} & {{v\over c^2}\over\sqrt{1+({v\over c})^2}} & 0 & 0\cr -{v\over\sqrt{1+({v\over c})^2}} & {1\over\sqrt{1+({v\over c})^2}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr} \to \mat{ 1 & 0 & 0 & 0\cr -v & 1 & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

Adding two boosts together:

.. math::

     A(u)A(v) = \mat{ {1\over\sqrt{1+{u^2\over c^2}}}& {{u\over c^2}\over\sqrt{1+{u^2\over c^2}}} & 0 & 0\cr -{u\over\sqrt{1+{u^2\over c^2}}} & {1\over\sqrt{1+{u^2\over c^2}}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr} \mat{ {1\over\sqrt{1+{v^2\over c^2}}}& {{v\over c^2}\over\sqrt{1+{v^2\over c^2}}} & 0 & 0\cr -{v\over\sqrt{1+{v^2\over c^2}}} & {1\over\sqrt{1+{v^2\over c^2}}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr} =


.. math::

    = \mat{ {1\over\sqrt{1+{w^2\over c^2}}}& {{w\over c^2}\over\sqrt{1+{w^2\over c^2}}} & 0 & 0\cr -{w\over\sqrt{1+{w^2\over c^2}}} & {1\over\sqrt{1+{w^2\over c^2}}} & 0 & 0\cr 0 & 0 & 1 & 0\cr 0 & 0 & 0 & 1\cr}

with

.. math::

    w = {u+v\over1-{uv\over c^2}}

However, there is one peculiar thing here that didn't exist in the O(3, 1)
case: by adding two velocities less than $c$, for example $u=v=c/2$, we get:

.. math::

    w = {c\over 1-{1\over 4}}={4c\over3}>c

(as opposed to $w = {c\over 1+{1\over 4}}={4c\over5}<c$ in the O(3, 1) case).
So one can get over $c$ easily. By adding $u=v={4c\over3}$ together:

.. math::

    w = {{8c\over 3}\over 1-{16\over 9}}=-{24c\over7}<0

(as opposed to $w = {{8c\over 3}\over 1+{16\over 9}}={24c\over25}>0$ in the
O(3, 1) case). So we can also get to negative speeds easily. One also needs to
be careful with identifying $\cos\phi_4 = {1\over\sqrt{1+({v\over c})^2}}$,
because for $\varphi_4>\pi/2$ we should probably set $\cos\varphi_4 =
-{1\over\sqrt{1+({v\over c})^2}}$. All of this follows directly from the
structure of SO(4), because one can get from $\Lambda^0{}_0>0$ to
$\Lambda^0{}_0<0$ continuously (this corresponds to increasing $\varphi_4$ over
$\pi/2$). In fact, by adding two speeds $u=v>c(\sqrt 2 - 1)$, one always gets
$w>c$. But if $c(\sqrt 2 - 1)\doteq0.414c$ is larger than any speed that we are
concerned about, we are fine.

.. index:: proper time

Proper Time
-----------


Proper time $\tau$ is a time elapsed by (physical) clocks along some (4D)
trajectory.  Coordinate time $t$ is just some time coordinate assigned to each
point in the space and usually one can find some real clocks, that would
measure such a time (many times they are in the infinity). To find a formula
for a proper time (in terms of the coordinate time), we introduce a local
inertial frame at each point of the trajectory -- in this frame, the clocks do
not move, e.g. $x$, $y$, $z$ is
constant (zero) and there is no gravity (this follows from the definition of
the local inertial frame), so the metric is just a Minkowski metric.

For any metrics, $\d s^2$ is invariant:

.. math::

    \d s^2 = g_{\mu\nu} \d x^\mu\d x^\nu

so coming to the local inertial frame, we have $x$, $y$, $z$ constant and we
get:

.. math::

    \d s^2 = g_{00} \d\tau^2

so:

.. math::

    \d\tau=\sqrt{\d s^2\over g_{00}}

since we are still in the local inertial frame (e.g. no gravity), we have
$g_{00}=-c^2$ (depending on which metrics we take it could also be $+c^2$), so:

.. math::

    \d\tau=\sqrt{-{\d s^2\over c^2}}

This formula was derived in the local inertial frame, but the right hand side
is the same in any inertial frame, because $\d s^2$ is invariant and $c$ too.
So in any frame we have:

.. math::

    \d\tau=\sqrt{-{\d s^2\over c^2}} =\sqrt{-{g_{\mu\nu} \d x^\mu\d x^\nu\over c^2}}


We'll explain how to calculate the proper time on the 1971 Hafele and Keating
experiment. They transported cesium-beam atomic clocks around the Earth on
scheduled commercial flights (once flying eastward, once westward) and compared
their reading on return to that of a standard clock at rest on the Earth's
surface.

We'll calculate it with all the metrics discussed above, to see the difference.

Weak Field Metric
~~~~~~~~~~~~~~~~~

Let's start with the metrics:

.. math::

    \d s^2=-\left(1+{2\phi\over c^2}\right)c^2 \d t^2 +\left(1-{2\phi\over c^2}\right)(\d x^2 +\d y^2 +\d z^2)

Then:

.. math::

    \tau_{AB} =\int_A^B\d\tau =\int_A^B\sqrt{-{\d s^2\over c^2}} =\int_A^B\sqrt{\left(1+{2\phi\over c^2}\right)\d t^2 -{1\over c^2}\left(1-{2\phi\over c^2}\right)(\d x^2 +\d y^2 +\d z^2)}=


.. math::

     =\int_A^B\d t\sqrt{\left(1+{2\phi\over c^2}\right) -{1\over c^2}\left(1-{2\phi\over c^2}\right)\left( \left(\d x\over\d t\right)^2 + \left(\d y\over\d t\right)^2 + \left(\d z\over\d t\right)^2\right)}=


.. math::

     =\int_A^B\d t\sqrt{\left(1+{2\phi\over c^2}\right) -{1\over c^2}\left(1-{2\phi\over c^2}\right)|{\bf V}|^2}

where

.. math::

     |{\bf V}|^2= \left(\d x\over\d t\right)^2 + \left(\d y\over\d t\right)^2 + \left(\d z\over\d t\right)^2

is the nonrelativistic velocity. Then we expand the square root into power
series and only keep terms with low powers of $c$:

.. math::

    \tau_{AB} =\int_A^B\d t\sqrt{\left(1+{2\phi\over c^2}\right) -{1\over c^2}\left(1-{2\phi\over c^2}\right)|{\bf V}|^2} =\int_A^B\d t\left(1+{\phi\over c^2}-{1\over 2c^2}|{\bf V}|^2\right)

so

.. math::

    \tau_{AB} =\int_A^B\d t\left(1-{1\over c^2}\left({1\over2}|{\bf V}|^2-\phi\right)\right)


Now let $V_g=V_g(t)$ be the speed of the plane relative to the (rotating) Earth
(positive for the eastbound flights, negative for the westbound ones),
$V_\oplus={2\pi R_\oplus\over24}\,{\rm1\over h}$ the surface speed of the
Earth, then the proper time for the clocks on the surface is:

.. math::

    \tau_\oplus =\int_A^B\d t\left(1-{1\over c^2}\left({1\over2}V_\oplus^2-\phi_\oplus\right) \right)

and for the clocks in the plane

.. math::

    \tau =\int_A^B\d t\left(1-{1\over c^2}\left({1\over2}(V_g+V_\oplus)^2-\phi\right) \right)

then the difference between the proper times is:

.. math::

    \tau-\tau_\oplus=\Delta\tau ={1\over c^2}\int_A^B\d t\left(-{1\over2}(V_g+V_\oplus)^2+\phi +{1\over2}V_\oplus^2-\phi_\oplus\right) ={1\over c^2}\int_A^B\d t\left( \phi-\phi_\oplus-{1\over2}V_g(V_g+2V_\oplus) \right)

but $\phi-\phi_\oplus=g h$, where $h=h(t)$ is the altitude of the plane, so
the final formula is:

.. math::

    \Delta\tau ={1\over c^2}\int_A^B\d t\left( gh-{1\over2}V_g(V_g+2V_\oplus) \right)

Let's evaluate it for typical altitudes and speeds of commercial aircrafts:

.. math::

    R_\oplus=6 378.1{\rm\,km}=6.3781\cdot10^6{\rm\,m}


.. math::

    V_\oplus={2\pi R_\oplus\over24}\,{\rm1\over h} ={2\pi R_\oplus\over24\cdot3600}\,{\rm1\over s} ={2\pi\,6.3781\cdot10^6\over24\cdot3600}{\rm m\over s}=463.83\rm\,{m\over s}


.. math::

    V_g = 870\,{\rm km\over h}=241.67\rm\,{m\over s}


.. math::

    h = 12{\rm\,km}=12000\rm\,m


.. math::

    t = {2\pi R_\oplus\over V_g} = {2\pi\,6.3781\cdot10^6\over 241.67}{\rm\,s} =165824.41{\rm\,s}\approx 46{\rm\,h}


.. math::

    c = 3\cdot10^8\rm\,{m\over s}

For
eastbound flights we get:

.. math::

    \Delta\tau ={t\over c^2} \left( gh-{1\over2}V_g(V_g+2V_\oplus) \right) =-4.344\cdot10^{-8}{\rm\,s}=-43.44{\rm\,ns}

and for westbound flights we get:

.. math::

    \Delta\tau ={t\over c^2} \left( gh-{1\over2}V_g(V_g-2V_\oplus) \right) =3.6964\cdot10^{-7}{\rm\,s}=369.63{\rm\,ns}

By neglecting gravity, one would get:
eastbound flights:

.. math::

    \Delta\tau ={t\over c^2} \left(-{1\over2}V_g(V_g+2V_\oplus) \right) =-260.34{\rm\,ns}

and for westbound flights:

.. math::

    \Delta\tau ={t\over c^2} \left(-{1\over2}V_g(V_g-2V_\oplus) \right) =152.73{\rm\,ns}

By just taking the clocks to the altitude $12\rm\,km$ and staying there for 46
hours (without moving with respect to the inertial frame, e.g. far galaxies), one gets:

.. math::

    \Delta\tau ={ght\over c^2}=216.90{\rm\,ns}


Rotating Disk Metric
~~~~~~~~~~~~~~~~~~~~

The rotating disk metrics is (taking weak field gravitation into account):

.. math::

    \d s^2=-\left(1+{2\phi\over c^2}-{\omega^2\over c^2}(x^2+y^2)\right)c^2 \d t^2 +(\d x^2 +\d y^2 +\d z^2)-2\omega y\,\d x\d t + 2\omega x\,\d y\d t

Then:

.. math::

    \tau_{AB} =\int_A^B\d\tau =\int_A^B\sqrt{-{\d s^2\over c^2}}=


.. math::

     =\int_A^B\sqrt{\left(1+{2\phi\over c^2}-{\omega^2\over c^2}(x^2+y^2)\right) \d t^2 -{1\over c^2}(\d x^2 +\d y^2 +\d z^2) +{2\omega y\over c^2}\,\d x\d t - {2\omega x\over c^2}\,\d y\d t }=


.. math::

     =\int_A^B\d t\sqrt{\left(1+{2\phi\over c^2}-{\omega^2\over c^2}(x^2+y^2)\right) -{1\over c^2}|{\bf V}|^2 +{2\omega y\over c^2}\,{\d x\over\d t} - {2\omega x\over c^2}\,{\d y\over\d t}}

where

.. math::

     |{\bf V}|^2= \left(\d x\over\d t\right)^2 + \left(\d y\over\d t\right)^2 + \left(\d z\over\d t\right)^2

is the nonrelativistic velocity. Then we expand the square root into power
series and only keep terms with low powers of $c$:

.. math::

    \tau_{AB} =\int_A^B\d t\left(1+{\phi\over c^2}-{1\over 2c^2}|{\bf V}|^2 +{\omega y\over c^2}\,{\d x\over\d t} - {\omega x\over c^2}\,{\d y\over\d t} \right)

so

.. math::

    \tau_{AB} =\int_A^B\d t\left(1-{1\over c^2}\left({1\over2}|{\bf V}|^2-\phi -{\omega y}\,{\d x\over\d t} + {\omega x}\,{\d y\over\d t} \right)\right)


Now as before let $V_g=V_g(t)$ be the speed of the plane (relative to the
rotating Earth, e.g. relative to our frame), $V_\oplus={2\pi
R_\oplus\over24}\,{\rm1\over h}$ the surface speed of the Earth, so $\omega
R_\oplus=V_\oplus$. For the clocks on the surface, we have:

.. math::

    x = R_\oplus


.. math::

    y = 0


.. math::

    z = 0

so

.. math::

    {\d x\over\d t}={\d y\over\d t}={\d z\over\d t}=0


.. math::

    |{\bf V}|^2=0

then the proper time for the clocks on the surface is:

.. math::

    \tau_\oplus =\int_A^B\d t\left(1-{1\over c^2}\left(-\phi_\oplus\right) \right)

and for the clocks in the plane we have:

.. math::

    x = (R_\oplus+h)\cos\Omega t


.. math::

    y = (R_\oplus+h)\sin\Omega t


.. math::

    z = 0

where $\Omega$ is defined by $\Omega (R_\oplus+h)=V_g$, so

.. math::

    {\d x\over\d t}=-(R_\oplus+h)\Omega\sin\Omega t


.. math::

    {\d y\over\d t}=(R_\oplus+h)\Omega\cos\Omega t


.. math::

    {\d z\over\d t}=0


.. math::

    |{\bf V}|^2=\Omega^2(R_\oplus+h)^2


.. math::

    \omega y {\d x\over\d t}=-\omega\Omega(R_\oplus+h)^2\sin^2\Omega t


.. math::

    \omega x {\d y\over\d t}=\omega\Omega(R_\oplus+h)^2\cos^2\Omega t

and

.. math::

    \tau =\int_A^B\d t\left(1-{1\over c^2}\left({1\over2} \Omega^2(R_\oplus+h)^2 - \phi +\omega\Omega(R_\oplus+h)^2\right) \right)

then the difference between the proper times is:

.. math::

    \tau-\tau_\oplus=\Delta\tau ={1\over c^2}\int_A^B\d t\left(-{1\over2}\Omega^2(R_\oplus+h)^2 -\omega\Omega(R_\oplus+h)^2 +\phi-\phi_\oplus\right) =


.. math::

     ={1\over c^2}\int_A^B\d t\left( -{1\over2}V_g^2 -V_\oplus V_g \left(1+{h\over R_\oplus}\right) +\phi-\phi_\oplus \right) =


.. math::

     ={1\over c^2}\int_A^B\d t\left( \phi-\phi_\oplus -{1\over2}V_g\left(V_g+2V_\oplus\left(1+{h\over R_\oplus}\right)\right) \right)

but $\phi-\phi_\oplus=g h$, where $h=h(t)$ is the altitude of the plane and we
approximate

.. math::

    \left(1+{h\over R_\oplus}\right)\approx 1 \,,

so
the final formula is the same as before:

.. math::

    \Delta\tau ={1\over c^2}\int_A^B\d t\left( gh-{1\over2}V_g(V_g+2V_\oplus) \right)

Note: for the values above, the bracket $\left(1+{h\over
R_\oplus}\right)^2\doteq1.00377$, so it's effect on the final difference of the
proper times is negligible (e.g. less than $1 \rm\,ns$). The difference is
caused by a slightly vague definition of the speed of the plane, e.g. the
ground speed is a bit different to the speed relative to the rotating Earth
(this depends on how much the atmosphere rotates with the Earth).

Concluding Remarks
~~~~~~~~~~~~~~~~~~

The coordinate time $t$ in both cases above is totally different. One can find
some physical clocks in both cases that measure (e.g. whose proper time is) the
particular coordinate time, but the beauty of the differential geometry
approach is that we don't have to care about this. $t$ is just a coordinate,
that we use to calculate something physical, like a proper time along some
trajectory, which is a frame invariant quantity. In both cases above, we got a
different formulas for the proper time of the surface clocks (and the clocks in
the plane) in terms of the coordinate time (because the coordinate time is
different in both cases), however the difference of the proper times is the
same in both cases:

.. math::

    \Delta\tau ={1\over c^2}\int_A^B\d t\left( gh-{1\over2}V_g(V_g+2V_\oplus) \right)

There is still a slight difference though -- the $t$ here used to evaluate the
integral is different in both cases. To do it correctly, one should take the
total time as measured by any of the clocks and then use the right formula for
the proper time of the particular clock to convert to the particular coordinate
time. However, the difference is small, of the order of nanoseconds, so it's
negligible compared to the total flying time of 46 hours.

FAQ
---


**How does one incorporate the fact, that there are only two possible
transformations, into all of this?**
For more info, see: http://arxiv.org/abs/0710.3398.
Answer: in that article there are
actually three possible transformations, $K<0$ corresponds to O(4), $K>0$ to
O(3, 1) and $K=0$ to either of them in the limit $c\to\infty$.

**What is the real difference between the Newtonian physics and special
relativity?** E.g. how do we derive the Minkowski metrics, how do we know we
need to set $c=const$ and how do we incorporate gravity in it?
Answer:
there are only three possible groups of transformations: O(4), O(3, 1) and a
limit of either for $c\to\infty$. All three provide inequivalent predictions
for high speeds, so we just choose the right one by experiment. It happens to
be the O(3, 1). As to gravity, that can be incorporated in either of them.

Questions Without Answers (Yet)
-------------------------------


How can one reformulate the article http://arxiv.org/abs/0710.3398 into the
language of the O(4) and O(3, 1) groups above? Basically each assumption and
equation must have some counterpart in what we have said above. I'd like to
identify those explicitely.

What are all the possible metrics, that generate the Newtonian Christoffel
symbols?
(Several such are given above, but I want to know all of them) Probable answer:
all metrics, whose inverse reduces to $g^{\mu\nu}=\diag(0, 1, 1, 1)$ in the
limit $c\to\infty$. I would like to have an explicit proof of this though.

What is the role of the different metrics, that generate the same Christoffel
symbols in the limit ($c\to\infty$)? Can one inertial frame be given with one
and another frame with a different form of the metrics (e.g. one with
$g_{00}=c^2$ and the other one with $g_{00}=-c^2$?) Possible answer: there is
no transformation to convert a metrics with signature +4 to signature +2, so
one has to choose one and then all other inertial frames have the same one.

What are all the allowed transformations between inertial frames? If we assume
that the inertial frames are given with one given metrics (see the previous
question), then the answer is: representation of the O(3, 1) group if
$g_{00}=-c^2$ or O(4) group if $g_{00}=c^2$. But if one frame is $g_{00}=-c^2$
and we transform to another frame with $g_{00}=c^2$, then it is not clear what
happens. Possible answer: one has to choose some signature and stick to it,
see also the previous question.

What is the real difference between Newtonian physics and general relativity?
Given our formulation of Newtonian physics using the differential geometry, I
want to know what the physical differences are between all the three theories
are.
