"""
Utilities for Geodesic Module

Unit System: M-Units => G = c = M = 1

"""
import warnings

import numpy as np


def _energy(E, q, p, a, mu):
    """
    Utility function to compute Energy of the test particle,
    using `scipy.optimize.fsolve`

    Parameters
    ----------
    E : float
        Variable, which `scipy.optimize.fsolve` solves for
    q : array_like
        Length-3 Array, containing the initial 3-Position
    p: array_like
        Length-3 Array, containing the initial 3-Momentum
    a : float
        Dimensionless Spin Parameter of the Black Hole
        ``0 <= a <= 1``
    mu : float
        Rest Mass of the test particle
        ``mu = 0.`` for Null Geodesics
        ``mu = -1`` for Timelike Geodesics

    Returns
    -------
    float
        Energy of the test particle

    """
    return (
        a ** 4 * (-(E ** 2) + mu ** 2 + 2 * p[0] ** 2)
        + 8 * a * E * p[2] * q[0]
        - a ** 2
        * (
            2 * p[2] ** 2
            - 2 * p[1] ** 2
            + q[0]
            * (
                mu ** 2 * (2 - 3 * q[0])
                - 4 * p[0] ** 2 * (-2 + q[0])
                + E ** 2 * (2 + 3 * q[0])
            )
        )
        + 2
        * q[0]
        * (
            p[1] ** 2 * (-2 + q[0])
            + q[0]
            * (
                p[0] ** 2 * (-2 + q[0]) ** 2
                + q[0] * (mu ** 2 * (-2 + q[0]) - E ** 2 * q[0])
            )
        )
        - (a ** 2 + (-2 + q[0]) * q[0])
        * (
            a ** 2 * (E - mu) * (E + mu) * np.cos(2 * q[1])
            - 2 * p[2] ** 2 * (1 / (np.sin(q[1]) ** 2))
        )
    ) / (
        (a ** 2 + (-2 + q[0]) * q[0])
        * (a ** 2 + 2 * q[0] ** 2 + a ** 2 * np.cos(2 * q[1]))
    )


def _sphToCart(r, th, ph):
    """
    Utility function to convert Spherical Polar
    Coordinates to Cartesian Coordinates

    Parameters
    ----------
    r : float
        r-component of 3-Position
    th : float
        theta-component of 3-Position
    ph : float
        phi-component of 3-Position

    Returns
    -------
    tuple
        3-Tuple, containing the coordinates in cartesian form

    """
    xs = r * np.sin(th) * np.cos(ph)
    ys = r * np.sin(th) * np.sin(ph)
    zs = r * np.cos(th)

    return xs, ys, zs


def _f_vec(ld, y, params):
    """
    Evaluates expressions for the RHS of the dynamic equations,
    from the Kerr Hamiltonian, to be used with ``_verlet_step()``

    Source: `Fuerst and Wu, 2004: <hhttps://www.aanda.org/articles/aa/pdf/2004/36/aa0814.pdf>`_

    Parameters
    ----------
    ld : float
        Affine Parameter, Lambda (Step-size)
    y : numpy.ndarray
        Length-6 array, containing the initial position and momentum of the test particle
    params : array_like
        Length-3 Array, containing - Black Hole Spin Parameter, `a`, Test Particle Energy, `E` and
        Test Particle Rest Mass, `mu`

    Returns
    -------
    yn : numpy.ndarray
        Length-6 array, containing the RHS of the dynamic equations

    """
    yn = np.zeros(shape=y.shape)

    r, th, ph = y[:3]
    pr, pth, pph = y[3:]
    a, E, mu = params

    r2 = r ** 2
    pr2 = pr ** 2
    pth2 = pth ** 2
    pph2 = pph ** 2
    a2 = a ** 2
    E2 = E ** 2
    sint = np.sin(th)
    sint2 = sint ** 2
    cost = np.cos(th)
    sg = r2 + a2 * (cost ** 2)
    dl = r2 - 2 * r + a2
    kappa = pth2 + pph2 * (1 / sint2) + a2 - (E2 * sint2 + mu)

    # Eqs. 21, 22, 27, 30, 31 in Source
    yn[0] = (dl * pr) / sg
    yn[1] = pth / sg
    yn[2] = (a * (-a * pph + 2 * E * r) + pph * (1 / sint2) * dl) / (dl * sg)
    yn[3] = (1 / (sg * dl)) * (
        ((r2 + a2) * mu - kappa) * (r - 1)
        + r * dl * mu
        + 2 * r * (r2 + a2) * E2
        - 2 * a * E * pph
    ) - (2 * pr2 * (r - 1)) / sg
    yn[4] = ((sint * cost) / sg) * ((pph2 / (sint ** 4)) - a2 * (E2 + mu))
    yn[5] = 0.0

    return yn


def _verlet_step(ld, y, params):
    """
    Synchornized Symplectic VerletLeapfrog Integrator
    Advances integration by one step

    Source: `Wikipedia: <https://en.wikipedia.org/wiki/Leapfrog_integration#Algorithm>`_

    Parameters
    ----------
    ld : float
        Affine Parameter, Lambda (Step-size)
    y : numpy.ndarray
        Length-6 array, containing the initial position and momentum of the test particle
    params : array_like
        Length-3 Array, containing - Black Hole Spin Parameter, `a`, Test Particle Energy, `E` and
        Test Particle Rest Mass, `mu`

    Returns
    -------
    y_next : numpy.ndarray
        Length-6 array, containing the values of 3-Position and 3-Momentum,
        corresponding to the integration step

    """
    DIM = 3
    y_next = np.copy(y)

    acc1 = _f_vec(ld, y, params)[3:]

    for i in range(DIM):
        y_next[i] = y[i] + ld * y[i + DIM] + 0.5 * ld * ld * acc1[i]

    acc2 = _f_vec(ld, y_next, params)[3:]

    for i in range(DIM):
        y_next[i + DIM] = y[i + DIM] + 0.5 * ld * (acc1[i] + acc2[i])

    return y_next


def _python_solver(q, p, params, end_lambda, step_size):
    """
    Wrapper to VerletLeapfrog Integrator, defined by ``_verlet_step()``

    This backend is currently in beta and the solver may not be stable for
    certain sets of conditions, e.g. long simulations (`end_lambda > 50.`)
    or high initial radial distances (`position[0] > ~5.`). In these cases,
    or if the output does not seem accurate, it is highly recommended to switch
    to the Julia backend, by setting `julia=True`, in the constructor call to
    ``einsteinpy.geodesic.*``.

    Parameters
    ----------
    q : array_like
        Length-3 Array, containing the initial 3-Position
    p : array_like
        Length-3 Array, containing the initial Covariant 3-Momentum
    params : array_like
        Length-3 Array, containing Black Hole Spin Parameter, `a`, Test Particle Energy, `E` and
        Test Particle Rest Mass, `mu`
    end_lambda : float
        Affine Parameter value, where integration will end
    step_size : float
        Step Size (Fixed)

    Returns
    -------
    lambdas : numpy.ndarray
        Array, containing affine parameter values, where integration was performed
    vecs : numpy.ndarray
        2D Array, containing integrated 3-Positions and Covariant 3-Momenta

    """
    state = np.hstack((q, p))
    a = params[0]

    y = state
    next_step = 0
    lambdas = list()
    vecs = list()

    lambdas.append(next_step)
    vecs.append(list(y))

    outer_event_horizon = 1 + np.sqrt(1 - a ** 2)

    # Integrating and storing results from solver
    while next_step <= end_lambda:
        y = _verlet_step(step_size, y, params)
        next_step += step_size
        lambdas.append(next_step)
        vecs.append(list(y))

        if np.abs(y[0]) <= np.abs(1.01 * outer_event_horizon):
            warnings.warn(
                "Test particle has reached the Event Horizon. ", RuntimeWarning
            )
            break

    lambdas = np.array(lambdas)
    vecs = np.array(vecs)

    return lambdas, vecs
