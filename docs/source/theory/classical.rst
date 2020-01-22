.. index:: Newtonian physics, classical mechanics

===================
Classical Mechanics
===================

Rigid Body Rotation
===================

In all the sections below, the rigid body is rotating around
the $\bomega$ axis, so:

.. math::

    {\bf v} = \bomega \times {\bf r}

Kinetic Energy
--------------

The kinetic energy is:

.. math::

    T = \int \half\rho({\bf r}) v^2 \d^3 r =

      = \int \half\rho({\bf r}) {\bf v}\cdot{\bf v} \d^3 r =

      = \int \half\rho({\bf r}) {\bf v}\cdot(\bomega \times {\bf r}) \d^3 r =

      = \int \half\rho({\bf r}) \bomega\cdot({\bf r}\times {\bf v}) \d^3 r =

      = \half \bomega \cdot \int\rho({\bf r}) ({\bf r}\times {\bf v}) \d^3 r =

      = \half \bomega \cdot {\bf L}

where ${\bf L}$ is the total angular momentum:

.. math::

      {\bf L} = \int\rho({\bf r}) ({\bf r}\times {\bf v}) \d^3 r

Angular Momentum
----------------

Total angular momentum is:

.. math::

    {\bf L}
        = \int \rho({\bf r}) ({\bf r} \times {\bf v}) \d^3 r =

        = \int \rho({\bf r}) ({\bf r} \times (\bomega \times {\bf r}))
                \d^3 r=

        = \int \rho({\bf r}) (\bomega r^2 - {\bf r} ({\bf r}
                \cdot \bomega)) \d^3 r =

        = \int \rho({\bf r}) (\one r^2 - {\bf r} {\bf r})
                \d^3 r \cdot \bomega =

        = {\bf I} \cdot \bomega

Where ${\bf I}$ is the moment of inertia tensor:

.. math::

    {\bf I} = \int \rho({\bf r}) (\one r^2 - {\bf r} {\bf r}) \d^3 r

Moment of Inertia
-----------------

The moment of inertia tensor and its components are:

.. math::

    {\bf I} = \int \rho({\bf r}) (\one r^2 - {\bf r} {\bf r}) \d^3 r

    I^{ij} = \int \rho({\bf r}) (\delta^{ij} r_k r^k - r^i r^j) \d^3 r

Let's write $\bomega=\omega {\bf n}$ (where ${\bf n}$ is a unit vector),
then the kinetic energy is:

.. math::

    T = \half \bomega \cdot {\bf L}
      = \half \bomega \cdot {\bf I} \cdot \bomega
      = \half {\bf n} \cdot {\bf I} \cdot {\bf n}\, \omega^2
      = \half I \omega^2

where $I$ is the moment of inertia about the axis of rotation:

.. math::

    I = {\bf n} \cdot {\bf I} \cdot {\bf n} =

      = {\bf n} \cdot \int \rho({\bf r}) (\one r^2 - {\bf r} {\bf r}) \d^3 r
        \cdot {\bf n} =

      = \int \rho({\bf r}) (r^2 - ({\bf r}\cdot {\bf n})^2) \d^3 r

Cylinder
^^^^^^^^

Solid cylinder of radius $R$, height $h$ and mass $m$. We'll use cylindrical
coordinates. First for rotation about the $z$ axis:

.. math::

    V = \pi R^2 h

    {\bf n} = (0, 0, 1)

    {\bf r} = (\rho\cos\phi, \rho\sin\phi, z)

    {\bf r} \cdot {\bf n} = z

    r^2 = \rho^2 + z^2


    I = \int \rho({\bf r}) (r^2 - ({\bf r}\cdot {\bf n})^2) \d^3 r
      = \int {m\over V} (\rho^2+z^2 - z^2) \d^3 r =

      = \int {m\over V} \rho^2 \d^3 r
      = {m\over V} \int_0^{2\pi}\d\phi \int_0^R\d R \int_{-{h\over2}}^{h\over2}
        \d z
         \rho^2 \rho =

      = {m\over V} 2\pi {R^4\over 4} h
      = {m\over \pi R^2 h} 2\pi {R^4\over 4} h
      = \half m R^2

Code::

    >>> from sympy import var, integrate, pi
    >>> var("m V R rho z phi h")
    (m, V, R, rho, z, phi, h)
    >>> I = m/V * integrate(rho**2 * rho, (rho, 0, R), (phi, 0, 2*pi), (z, -h/2, h/2))
    >>> I.subs(V, pi * R**2 * h)
    R**2*m/2


And about the $x$ axis:

.. math::

    {\bf n} = (1, 0, 0)

    {\bf r} = (\rho\cos\phi, \rho\sin\phi, z)

    {\bf r} \cdot {\bf n} = \rho\cos\phi

    r^2 = \rho^2 + z^2


    I = \int \rho({\bf r}) (r^2 - ({\bf r}\cdot {\bf n})^2) \d^3 r
      = \int {m\over V} (\rho^2+z^2 - \rho^2\cos^2\phi) \d^3 r =

      = {m\over V} \int_0^{2\pi}\d\phi \int_0^R\d R \int_{-{h\over2}}^{h\over2}
        \d z (\rho^2+z^2 - \rho^2\cos^2\phi)\rho =

      = {m\over V}\left({\pi R^4 h\over 2}+{\pi R^2 h^3\over 12}
                    -{\pi R^4 h\over 4}\right) =

      = {m\over \pi R^2 h}\left({\pi R^4 h\over 2}+{\pi R^2 h^3\over 12}
                    -{\pi R^4 h\over 4}\right) =

      = {m\over 12} (6R^2 + h^2 - 3R^2) =

      = {m\over 12} (3R^2 + h^2)

Code::

    >>> from sympy import var, integrate, pi, cos
    >>> var("m V R rho z phi h")
    (m, V, R, rho, z, phi, h)
    >>> I = m/V * integrate((rho**2+z**2-rho**2*cos(phi)**2) * rho, (rho, 0, R), (phi, 0, 2*pi), (z, -h/2, h/2))
    >>> I.subs(V, pi * R**2 * h).simplify()
    m*(3*R**2 + h**2)/12

Special cases are a rod of length $h$ (set $R=0$ above) and a thin solid disk
of radius $R$ and mass $m$ (set $h=0$ above).

Sphere
^^^^^^

Solid sphere of radius $R$ and mass $m$. We'll use spherical
coordinates. All axes are equivalent, so we use
rotation about the $z$ axis:

.. math::

    V = {4\over3} \pi R^3

    {\bf n} = (0, 0, 1)

    {\bf r} = (\rho\cos\phi\sin\theta, \rho\sin\phi\sin\theta, \rho\cos\theta)

    {\bf r} \cdot {\bf n} = \rho\cos\theta

    r^2 = \rho^2


    I = \int \rho({\bf r}) (r^2 - ({\bf r}\cdot {\bf n})^2) \d^3 r
      = \int {m\over V} (\rho^2 - \rho^2\cos^2\theta) \d^3 r =

      = {m\over V} \int_0^{2\pi}\d\phi \int_0^R\d R \int_0^\pi
        \d \theta
         \rho^2(1-\cos^2\theta) \rho^2\sin\theta =

      = {m\over V} \int_0^{2\pi}\d\phi \int_0^R\d R \int_0^\pi
        \d \theta
         \rho^4\sin^3\theta =

      = {m\over V} 2\pi {R^5\over 5} {4\over 3} =

      = {m\over V} {8\over 15}\pi R^5
      = {m\over {4\over 3}\pi R^3} {8\over 15}\pi R^5
      = {2\over 5} m R^2

Code::

    >>> from sympy import var, integrate, pi, sin
    >>> var("m V R rho theta phi")
    (m, V, R, rho, theta, phi)
    >>> I = m/V * integrate(rho**4 * sin(theta)**3, (rho, 0, R), (phi, 0, 2*pi), (theta, 0, pi))
    >>> I
    8*pi*R**5*m/(15*V)
    >>> I.subs(V, 4*pi*R**3/3)
    2*R**2*m/5
