import astropy.units as u
import numpy as np

from einsteinpy import constant, utils

nonzero_christoffels_list = [
    (0, 0, 1),
    (0, 0, 2),
    (0, 1, 3),
    (0, 2, 3),
    (0, 1, 0),
    (0, 2, 0),
    (0, 3, 1),
    (0, 3, 2),
    (1, 0, 0),
    (1, 1, 1),
    (1, 2, 2),
    (1, 3, 3),
    (2, 0, 0),
    (2, 1, 1),
    (2, 2, 2),
    (2, 3, 3),
    (1, 0, 3),
    (1, 1, 2),
    (2, 0, 3),
    (2, 1, 2),
    (1, 2, 1),
    (1, 3, 0),
    (2, 2, 1),
    (2, 3, 0),
    (3, 0, 1),
    (3, 0, 2),
    (3, 1, 0),
    (3, 1, 3),
    (3, 2, 0),
    (3, 2, 3),
    (3, 3, 1),
    (3, 3, 2),
]  #: Precomputed list of tuples consisting of indices of christoffel symbols which are non-zero in Kerr-Newman Metric


def charge_length_scale(
    Q, c=constant.c.value, G=constant.G.value, Cc=constant.coulombs_const.value
):
    """
    Returns a length scale corrosponding to the Electric Charge Q of the mass

    Parameters
    ----------
    Q : float
        Charge on the massive body
    c : float
        Speed of light. Defaults to 299792458 (SI units)
    G : float
        Gravitational constant. Defaults to 6.67408e-11 (SI units)
    Cc : float
        Coulumb's constant. Defaults to 8.98755e9 (SI units)

    Returns
    -------
    float
        returns (coulomb's constant^0.5)*(Q/c^2)*G^0.5

    """
    return (Q / (c ** 2)) * np.sqrt(G * Cc)


def rho(r, theta, a):
    """
    Returns the value sqrt(r^2 + a^2 * cos^2(theta)).
    Specific to Boyer-Lindquist coordinates

    Parameters
    ----------
    r : float
        Component r in vector
    theta : float
        Component theta in vector
    a : float
        Any constant

    Returns
    -------
    float
        The value sqrt(r^2 + a^2 * cos^2(theta))
    
    """
    return np.sqrt((r ** 2) + ((a * np.cos(theta)) ** 2))


def delta(
    r, M, a, Q, c=constant.c.value, G=constant.G.value, Cc=constant.coulombs_const.value
):
    """
    Returns the value r^2 - Rs * r + a^2
    Specific to Boyer-Lindquist coordinates

    Parameters
    ----------
    r : float
        Component r in vector
    M : float
        Mass of the massive body
    a : float
        Any constant
    Q : float
        Charge on the massive body
    c : float
        Speed of light
    G : float
        Gravitational constant
    Cc : float
        Coulomb's constant

    Returns
    -------
    float
        The value r^2 - Rs * r + a^2 + Rq^2
    
    """
    Rs = utils.schwarzschild_radius_dimensionless(M, c, G)
    return (r ** 2) - (Rs * r) + (a ** 2) + (charge_length_scale(Q, c, G, Cc) ** 2)


def metric(
    r,
    theta,
    M,
    a,
    Q,
    c=constant.c.value,
    G=constant.G.value,
    Cc=constant.coulombs_const.value,
):
    """
    Returns the Kerr-Newman Metric

    Parameters
    ----------
    
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    M : float
        Mass of the massive body
    a : float
        Black Hole spin factor
    Q : float
        Charge on the massive body
    c : float
        Speed of light
    G : float
        Gravitational constant
    Cc : float
        Coulomb's constant
    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4)

    """
    m = np.zeros((4, 4), dtype=float)
    rh2, dl = rho(r, theta, a) ** 2, delta(r, M, a, Q, c, G, Cc)
    c2 = c ** 2
    # set the diagonal/off-diagonal terms of metric
    m[0, 0] = (dl - ((a * np.sin(theta)) ** 2)) / (rh2)
    m[1, 1] = -rh2 / (dl * c2)
    m[2, 2] = -rh2 / c2
    m[3, 3] = (
        (((a * np.sin(theta)) ** 2) * dl - ((r ** 2 + a ** 2) ** 2))
        * (np.sin(theta) ** 2)
        / (rh2 * c2)
    )
    m[0, 3] = m[3, 0] = (
        -a * (np.sin(theta) ** 2) * (dl - (r ** 2) - (a ** 2)) / (rh2 * c)
    )
    return m


def metric_inv(
    r,
    theta,
    M,
    a,
    Q,
    c=constant.c.value,
    G=constant.G.value,
    Cc=constant.coulombs_const.value,
):
    """
    Returns the inverse of Kerr-Newman Metric

    Parameters
    ----------
    
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    M : float
        Mass of the massive body
    a : float
        Black Hole spin factor
    Q : float
        Charge on the massive body
    c : float
        Speed of light
    G : float
        Gravitational constant
    Cc : float
        Coulomb's constant
    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4)

    """
    return np.linalg.inv(metric(r, theta, M, a, Q, c, G, Cc))


def dmetric_dx(
    r,
    theta,
    M,
    a,
    Q,
    c=constant.c.value,
    G=constant.G.value,
    Cc=constant.coulombs_const.value,
):
    """
    Returns differentiation of each component of Kerr-Newman metric tensor w.r.t. t, r, theta, phi

    Parameters
    ----------
    
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    M : float
        Mass of the massive body
    a : float
        Black Hole spin factor
    Q : float
        Charge on the massive body
    c : float
        Speed of light
    G : float
        Gravitational constant
    Cc : float
        Coulomb's constant
    
    Returns
    -------
    dmdx : ~numpy.array
        Numpy array of shape (4,4,4)
        dmdx[0], dmdx[1], dmdx[2] & dmdx[3] is differentiation of metric w.r.t. t, r, theta & phi respectively

    """
    Rs = utils.schwarzschild_radius_dimensionless(M, c, G)
    dmdx = np.zeros((4, 4, 4), dtype=float)
    rh2, dl = rho(r, theta, a) ** 2, delta(r, M, a, Q, c, G, Cc)
    c2 = c ** 2
    # metric is invariant on t & phi
    # differentiation of metric wrt r
    def due_to_r():
        nonlocal dmdx
        drh2dr = 2 * r
        dddr = 2 * r - Rs
        dmdx[1, 0, 0] = (dddr * rh2 - drh2dr * (dl - (a * np.sin(theta)) ** 2)) / (
            rh2 ** 2
        )
        dmdx[1, 1, 1] = (-1 / (c2 * (dl ** 2))) * (drh2dr * dl - dddr * rh2)
        dmdx[1, 2, 2] = -drh2dr / c2
        dmdx[1, 3, 3] = ((np.sin(theta) ** 2) / (c2 * (rh2 ** 2))) * (
            (
                (((a * np.sin(theta)) ** 2) * dddr - 4 * (r ** 3) - 4 * (r * (a ** 2)))
                * rh2
            )
            - (drh2dr * (((a * np.sin(theta)) ** 2) * dl - ((r ** 2 + a ** 2) ** 2)))
        )
        dmdx[1, 0, 3] = dmdx[1, 3, 0] = (
            (-a) * (np.sin(theta) ** 2) / (c * (rh2 ** 2))
        ) * ((dddr - 2 * r) * rh2 - drh2dr * (dl - r ** 2 - a ** 2))

    # differentiation of metric wrt theta
    def due_to_theta():
        nonlocal dmdx
        drh2dth = -2 * (a ** 2) * np.cos(theta) * np.sin(theta)
        dmdx[2, 0, 0] = (
            (-2 * (a ** 2) * np.sin(theta) * np.cos(theta)) * rh2
            - drh2dth * (dl - ((a * np.sin(theta)) ** 2))
        ) / (rh2 ** 2)
        dmdx[2, 1, 1] = -drh2dth / (c2 * dl)
        dmdx[2, 2, 2] = -drh2dth / c2
        dmdx[2, 3, 3] = (1 / (c2 * (rh2 ** 2))) * (
            (
                (
                    (4 * (a ** 2) * (np.sin(theta) ** 3) * np.cos(theta) * dl)
                    - (2 * np.sin(theta) * np.cos(theta) * ((r ** 2 + a ** 2) ** 2))
                )
                * rh2
            )
            - (
                drh2dth
                * (((a * np.sin(theta)) ** 2) * dl - ((r ** 2 + a ** 2) ** 2))
                * (np.sin(theta) ** 2)
            )
        )
        dmdx[2, 0, 3] = dmdx[2, 3, 0] = (
            (-a * (dl - r ** 2 - a ** 2)) / (c * (rh2 ** 2))
        ) * (
            (2 * np.sin(theta) * np.cos(theta) * rh2) - (drh2dth * (np.sin(theta) ** 2))
        )

    due_to_r()
    due_to_theta()
    return dmdx


def christoffels(
    r,
    theta,
    M,
    a,
    Q,
    c=constant.c.value,
    G=constant.G.value,
    Cc=constant.coulombs_const.value,
):
    """
    Returns the 3rd rank Tensor containing Christoffel Symbols for Kerr-Newman Metric

    Parameters
    ----------
    
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    M : float
        Mass of the massive body
    a : float
        Black Hole spin factor
    Q : float
        Charge on the massive body
    c : float
        Speed of light
    G : float
        Gravitational constant
    Cc : float
        Coulomb's constant
    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4,4)

    """
    invg = metric_inv(r, theta, M, a, Q, c, G, Cc)
    dmdx = dmetric_dx(r, theta, M, a, Q, c, G, Cc)
    chl = np.zeros(shape=(4, 4, 4), dtype=float)
    for _, k, l in nonzero_christoffels_list[0:4]:
        val1 = dmdx[l, 0, k] + dmdx[k, 0, l]
        val2 = dmdx[l, 3, k] + dmdx[k, 3, l]
        chl[0, k, l] = chl[0, l, k] = 0.5 * (invg[0, 0] * (val1) + invg[0, 3] * (val2))
        chl[3, k, l] = chl[3, l, k] = 0.5 * (invg[3, 0] * (val1) + invg[3, 3] * (val2))
    for i, k, l in nonzero_christoffels_list[8:16]:
        chl[i, k, l] = 0.5 * (
            invg[i, i] * (dmdx[l, i, k] + dmdx[k, i, l] - dmdx[i, k, l])
        )
    for i, k, l in nonzero_christoffels_list[16:20]:
        chl[i, k, l] = chl[i, l, k] = 0.5 * (
            invg[i, i] * (dmdx[l, i, k] + dmdx[k, i, l] - dmdx[i, k, l])
        )
    return chl


def em_potential(
    r,
    theta,
    a,
    Q,
    M,
    c=constant.c.value,
    G=constant.G.value,
    Cc=constant.coulombs_const.value,
):
    """
    Returns a 4-d vector(for each component of 4-d space-time) containing the electromagnetic potential around a Kerr-Newman body

    Parameters
    ----------
    
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    a : float
        Black Hole spin factor
    Q : float
        Charge on the massive body
    M : float
        Mass of the massive body
    c : float
        Speed of light
    G : float
        Gravitational constant
    Cc : float
        Coulomb's constant
    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,)
    
    """
    rq, rh2, c2 = charge_length_scale(Q, c, G, Cc), rho(r, theta, a) ** 2, c ** 2
    vec = np.zeros((4,), dtype=float)
    vec[0] = r * rq / rh2
    vec[3] = ((-c2) / (rh2 * G * M)) * a * r * rq * (np.sin(theta) ** 2)
    return vec


def maxwell_tensor_covariant(
    r,
    theta,
    a,
    Q,
    M,
    c=constant.c.value,
    G=constant.G.value,
    Cc=constant.coulombs_const.value,
):
    """
    Returns a 2nd rank Tensor containing Maxwell Tensor with lower indices for Kerr-Newman Metric

    Parameters
    ----------
    
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    a : float
        Black Hole spin factor
    Q : float
        Charge on the massive body
    M : float
        Mass of the massive body
    c : float
        Speed of light
    G : float
        Gravitational constant
    Cc : float
        Coulomb's constant
    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4)

    """
    c2 = c ** 2
    m = np.zeros((4, 4), dtype=float)
    rh2, rq = rho(r, theta, a) ** 2, charge_length_scale(Q, c, G, Cc)
    drh2dr, drh2dth = 2 * r, -2 * (a ** 2) * np.cos(theta) * np.sin(theta)
    # set the tensor terms
    m[0, 1] = ((-rq) / (rh2 ** 2)) * (rh2 - drh2dr * r)
    m[0, 2] = r * rq * drh2dth / (rh2 ** 2)
    m[3, 1] = (c2 * a * rq * (np.sin(theta) ** 2) / (G * M * (rh2 ** 2))) * (
        rh2 - r * drh2dr
    )
    m[3, 2] = (c2 * a * rq * r / (G * M * (rh2 ** 2))) * (
        (2 * np.sin(theta) * np.cos(theta) * rh2) - (drh2dth * (np.sin(theta) ** 2))
    )
    for i, j in [(0, 1), (0, 2), (3, 1), (3, 2)]:
        m[j, i] = -m[i, j]
    return m


def maxwell_tensor_contravariant(
    r,
    theta,
    a,
    Q,
    M,
    c=constant.c.value,
    G=constant.G.value,
    Cc=constant.coulombs_const.value,
):
    """
    Returns a 2nd rank Tensor containing Maxwell Tensor with upper indices for Kerr-Newman Metric

    Parameters
    ----------
    
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    a : float
        Black Hole spin factor
    Q : float
        Charge on the massive body
    M : float
        Mass of the massive body
    c : float
        Speed of light
    G : float
        Gravitational constant
    Cc : float
        Coulomb's constant
    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4)

    """
    mcov = maxwell_tensor_covariant(r, theta, a, Q, M, c, G, Cc)
    ginv = metric_inv(r, theta, M, a, Q, c, G, Cc)
    # contravariant F = contravariant g X covariant F X transpose(contravariant g)
    # but g is symettric
    return np.matmul(np.matmul(ginv, mcov), ginv)


@u.quantity_input(mass=u.kg, Q=u.C)
def kerrnewman_time_velocity(pos_vec, vel_vec, mass, a, Q):
    """
    Velocity of coordinate time wrt proper metric

    Parameters
    ----------
    pos_vector : ~numpy.array
        Vector with r, theta, phi components in SI units
    vel_vector : ~numpy.array
        Vector with velocities of r, theta, phi components in SI units
    mass : ~astropy.units.kg
        Mass of the body
    a : float
        Any constant
    Q : ~astropy.units.C
        Charge on the massive body

    Returns
    -------
    ~astropy.units.one
        Velocity of time

    """
    _scr = utils.schwarzschild_radius(mass).value
    Qc = Q.to(u.C)
    g = metric(pos_vec[0], pos_vec[1], mass.value, a, Qc.value)
    A = g[0, 0]
    B = 2 * g[0, 3]
    C = (
        g[1, 1] * (vel_vec[0] ** 2)
        + g[2, 2] * (vel_vec[1] ** 2)
        + g[3, 3] * (vel_vec[2] ** 2)
        - 1
    )
    D = (B ** 2) - (4 * A * C)
    vt = (B + np.sqrt(D)) / (2 * A)
    return vt * u.one
