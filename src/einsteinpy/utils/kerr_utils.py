import astropy.units as u
import numpy as np

from einsteinpy import constant
from einsteinpy.utils import schwarzschild_radius as scr

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
]  #: Precomputed list of tuples consisting of indices of christoffel symbols which are non-zero in Kerr Metric


def sigma(r, theta, a):
    """
    Returns the value r^2 + a^2 * cos^2(theta)
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
        The value r^2 + a^2 * cos^2(theta)
    
    """
    return (r ** 2) + ((a * np.cos(theta)) ** 2)


def delta(r, Rs, a):
    """
    Returns the value r^2 - Rs * r + a^2
    Specific to Boyer-Lindquist coordinates

    Parameters
    ----------
    r : float
        Component r in vector
    Rs : float
        Schwarzschild radius
    a : float
        Any constant

    Returns
    -------
    float
        The value r^2 - Rs * r + a^2
    
    """
    return (r ** 2) - (Rs * r) + (a ** 2)


def metric(c, r, theta, Rs, a):
    """
    Returns the Kerr Metric

    Parameters
    ----------
    c : float
        Speed of light
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    Rs : float
        Schwarzschild Radius
    a : float
        Black Hole spin factor

    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4)

    """
    m = np.zeros(shape=(4, 4), dtype=float)
    sg, dl = sigma(r, theta, a), delta(r, Rs, a)
    c2 = c ** 2
    # set the diagonal/off-diagonal terms of metric
    m[0][0] = 1 - (Rs * r / sg)
    m[1][1] = (sg / dl) * (-1 / c2)
    m[2][2] = -1 * sg / c2
    m[3][3] = (
        (-1 / c2)
        * ((r ** 2) + (a ** 2) + (Rs * r * (np.sin(theta) ** 2) * ((a ** 2) / sg)))
        * (np.sin(theta) ** 2)
    )
    m[0][3] = Rs * r * a * (np.sin(theta) ** 2) / (sg * c)
    m[3][0] = m[0][3]
    return m


def metric_inv(c, r, theta, Rs, a):
    """
    Returns the inverse of Kerr Metric

    Parameters
    ----------
    c : float
        Speed of light
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    Rs : float
        Schwarzschild Radius
    a : float
        Black Hole spin factor

    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4)

    """
    m = metric(c, r, theta, Rs, a)
    return np.linalg.inv(m)


def dmetric_dx(c, r, theta, Rs, a):
    """
    Returns differentiation of each component of Kerr metric tensor w.r.t. t, r, theta, phi

    Parameters
    ----------
    c : float
        Speed of light
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    Rs : float
        Schwarzschild Radius
    a : float
        Black Hole spin factor

    Returns
    -------
    dmdx : ~numpy.array
        Numpy array of shape (4,4,4)
        dmdx[0], dmdx[1], dmdx[2] & dmdx[3] is differentiation of metric w.r.t. t, r, theta & phi respectively

    """
    dmdx = np.zeros(shape=(4, 4, 4), dtype=float)
    sg, dl = sigma(r, theta, a), delta(r, Rs, a)
    c2 = c ** 2

    # metric is invariant on t & phi
    # differentiation of metric wrt r
    def due_to_r():
        nonlocal dmdx
        dsdr = 2 * r
        dddr = 2 * r - Rs
        tmp = (Rs * (sg - r * dsdr) / sg) * (1 / sg)
        dmdx[1][0][0] = -1 * tmp
        dmdx[1][1][1] = (-1 / c2) * (dsdr - (sg * (dddr / dl))) / dl
        dmdx[1][2][2] = (-1 / c2) * dsdr
        dmdx[1][3][3] = (
            (-1 / c2)
            * (2 * r + (a ** 2) * (np.sin(theta) ** 2) * tmp)
            * (np.sin(theta) ** 2)
        )
        dmdx[1][0][3] = (1 / c) * (a * (np.sin(theta) ** 2) * tmp)
        dmdx[1][3][0] = dmdx[1][0][3]

    # differentiation of metric wrt theta
    def due_to_theta():
        nonlocal dmdx
        dsdth = -2 * (a ** 2) * np.cos(theta) * np.sin(theta)
        tmp = (-1 / sg) * Rs * r * dsdth / sg
        dmdx[2][0][0] = -1 * tmp
        dmdx[2][1][1] = (-1 / c2) * (dsdth / dl)
        dmdx[2][2][2] = (-1 / c2) * dsdth
        dmdx[2][3][3] = (-1 / c2) * (
            2 * np.sin(theta) * np.cos(theta) * ((r ** 2) + (a ** 2))
            + tmp * (a ** 2) * (np.sin(theta) ** 4)
            + (4 * (np.sin(theta) ** 3) * np.cos(theta) * (a ** 2) * r * Rs / sg)
        )
        dmdx[2][0][3] = (a / c) * (
            (np.sin(theta) ** 2) * tmp
            + (2 * np.sin(theta) * np.cos(theta) * Rs * r / sg)
        )
        dmdx[2][3][0] = dmdx[2][0][3]

    due_to_r()
    due_to_theta()
    return dmdx


def christoffels(c, r, theta, Rs, a):
    """
    Returns the 3rd rank Tensor containing Christoffel Symbols for Kerr Metric

    Parameters
    ----------
    c : float
        Speed of light
    r : float
        Distance from the centre
    theta : float
        Angle from z-axis
    Rs : float
        Schwarzschild Radius
    a : float
        Black Hole spin factor

    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4,4)

    """
    invg = metric_inv(c, r, theta, Rs, a)
    dmdx = dmetric_dx(c, r, theta, Rs, a)
    chl = np.zeros(shape=(4, 4, 4), dtype=float)
    for _, k, l in nonzero_christoffels_list[0:4]:
        val1 = dmdx[l][0][k] + dmdx[k][0][l]
        val2 = dmdx[l][3][k] + dmdx[k][3][l]
        chl[0][k][l] = chl[0][l][k] = 0.5 * (invg[0][0] * (val1) + invg[0][3] * (val2))
        chl[3][k][l] = chl[3][l][k] = 0.5 * (invg[3][0] * (val1) + invg[3][3] * (val2))
    for i, k, l in nonzero_christoffels_list[8:12]:
        dmdxVal = dmdx[l]
        val1 =  0.5 * (
            invg[i][i] * (dmdxVal[i][k] + dmdxVal[i][l] - dmdx[i][k][l])
        )
        val2 =  0.5 * (
            invg[i+1][i+1] * (dmdxVal[i+1][k] + dmdxVal[i+1][l] - dmdx[i+1][k][l])
        )

        chl[i][k][l] = val1
        chl[i+1][k][l] = val2
        
    for i, k, l in nonzero_christoffels_list[16:20]:
        chl[i][k][l] = chl[i][l][k] = 0.5 * (
            invg[i][i] * (dmdx[l][i][k] + dmdx[k][i][l] - dmdx[i][k][l])
        )
    return chl


@u.quantity_input(mass=u.kg)
def kerr_time_velocity(pos_vec, vel_vec, mass, a):
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

    Returns
    -------
    ~astropy.units.one
        Velocity of time

    """
    _scr = scr(mass).value
    g = metric(constant.c.value, pos_vec[0], pos_vec[1], _scr, a)
    A = g[0][0]
    B = 2 * g[0][3]
    C = (
        g[1][1] * (vel_vec[0] ** 2)
        + g[2][2] * (vel_vec[1] ** 2)
        + g[3][3] * (vel_vec[2] ** 2)
        - 1
    )
    D = (B ** 2) - (4 * A * C)
    vt = (B + np.sqrt(D)) / (2 * A)
    return vt * u.one


def nonzero_christoffels():
    """
    Returns a list of tuples consisting of indices of christoffel symbols which are non-zero in Kerr Metric computed in real-time.

    Returns
    -------
    list
        List of tuples
        each tuple (i,j,k) represent christoffel symbol with i as upper index and j,k as lower indices.
    """
    # Below is the code for algorithmically calculating the indices of nonzero christoffel symbols in Kerr Metric.
    invg = np.zeros(shape=(4, 4), dtype=bool)
    dmdx = np.zeros(shape=(4, 4, 4), dtype=bool)
    invg[3][0] = invg[0][3] = True
    dmdx[1][3][0] = dmdx[1][0][3] = True
    dmdx[2][0][3] = dmdx[2][3][0] = True
    # Code Climate cyclomatic complexity hack ; Instead of "for i in range(4)"
    invg[0][0] = invg[1][1] = invg[2][2] = invg[3][3] = True
    dmdx[1][0][0] = dmdx[1][1][1] = dmdx[1][2][2] = dmdx[1][3][3] = True
    dmdx[2][0][0] = dmdx[2][1][1] = dmdx[2][2][2] = dmdx[2][3][3] = True
    # hack ends
    chl = np.zeros(shape=(4, 4, 4), dtype=bool)
    tmp = np.array([i for i in range(4 ** 3)])
    vcl = list()
    for t in tmp:
        i = int(t / (4 ** 2)) % 4
        k = int(t / 4) % 4
        l = t % 4
        for m in range(4):
            chl[i][k][l] |= invg[i][m] & (dmdx[l][m][k] | dmdx[k][m][l] | dmdx[m][k][l])
        if chl[i][k][l]:
            vcl.append((i, k, l))
    return vcl
