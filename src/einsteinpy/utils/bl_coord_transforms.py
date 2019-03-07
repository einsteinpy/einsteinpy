import numpy as np


def CartesianToBL_pos(pos_vec, a):
    """
    Function to convert Cartesian to Boyer-Lindquits coordinates (position vector)

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having x, y, z coordinates in SI units(m)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array with r, theta, phi in (m, rad, rad)

    """
    r_vec = np.zeros(shape=(3,), dtype=float)
    w = float(np.sum(np.square(pos_vec[:3])) - (a ** 2))
    r_vec[0] = np.sqrt(
        0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (pos_vec[2] ** 2))))
    )
    r_vec[1] = np.arccos(pos_vec[2] / r_vec[0])
    r_vec[2] = np.arctan2(pos_vec[1], pos_vec[0])
    return r_vec


def CartesianToBL_vel(pos_vec, vel_vec, a):
    """
    Function to convert velocity in Cartesian coordinates to velocity in Boyer-Lindquist coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having x, y, z coordinates in SI units(m)
    vel_vec : ~numpy.array
        3-length numpy array with vx, vy, vz in (m/s)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array having v_r, v_theta, v_phi in (m/s, rad/s, rad/s)
    
    """
    v_vec = np.zeros(shape=(3,), dtype=float)
    w = float(np.sum(np.square(pos_vec[:3])) - (a ** 2))
    dw_dt = float(2 * np.sum(np.multiply(pos_vec, vel_vec)))
    r = np.sqrt(0.5 * (w + np.sqrt((w ** 2) + (4 * (a ** 2) * (pos_vec[2] ** 2)))))
    v_vec[0] = (1 / (2 * r)) * (
        (dw_dt / 2)
        + (
            (w * dw_dt + 4 * (a ** 2) * pos_vec[2] * vel_vec[2])
            / (2 * np.sqrt((w ** 2) + (4 * (a ** 2) * (pos_vec[2] ** 2))))
        )
    )
    v_vec[1] = (-1 / np.sqrt(1 - np.square(pos_vec[2] / r))) * (
        (vel_vec[2] * r - v_vec[0] * pos_vec[2]) / (r ** 2)
    )
    v_vec[2] = (1 / (1 + np.square(pos_vec[1] / pos_vec[0]))) * (
        (vel_vec[1] * pos_vec[0] - vel_vec[0] * pos_vec[1]) / (pos_vec[0] ** 2)
    )
    return v_vec


def BLToCartesian_pos(pos_vec, a):
    """
    Function to convert Boyer-Lindquist coordinates to Cartesian coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having r, theta, phi coordinates in SI units(m, rad, rad)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array with x, y, z in m.

    """
    r_vec = np.zeros(shape=(3,), dtype=float)
    tmp = np.sqrt(pos_vec[0] ** 2 + a ** 2) * np.sin(pos_vec[1])
    r_vec[0] = tmp * np.cos(pos_vec[2])
    r_vec[1] = tmp * np.sin(pos_vec[2])
    r_vec[2] = pos_vec[0] * np.cos(pos_vec[1])
    return r_vec


def BLToCartesian_vel(pos_vec, vel_vec, a):
    """
    Function to convert velocities in Boyer-Lindquist coordinates to velocities in Cartesian coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having r, theta, phi coordinates in SI units(m, rad, rad)
    vel_vec : ~numpy.array
        3-length numpy array with V_r, V_theta, V_phi (m/s, rad/s, rad/s)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array having vx, vy, vz in SI units(m/s)

    """
    v_vec = np.zeros(shape=(3,), dtype=float)
    tmp = np.sqrt(pos_vec[0] ** 2 + a ** 2)
    v_vec[0] = (
        (pos_vec[0] * vel_vec[0] * np.sin(pos_vec[1]) * np.cos(pos_vec[2]) / tmp)
        + (tmp * np.cos(pos_vec[1]) * np.cos(pos_vec[2]) * vel_vec[1])
        - (tmp * np.sin(pos_vec[1]) * np.sin(pos_vec[2]) * vel_vec[2])
    )
    v_vec[1] = (
        (pos_vec[0] * vel_vec[0] * np.sin(pos_vec[1]) * np.sin(pos_vec[2]) / tmp)
        + (tmp * np.cos(pos_vec[1]) * np.sin(pos_vec[2]) * vel_vec[1])
        + (tmp * np.sin(pos_vec[1]) * np.cos(pos_vec[2]) * vel_vec[2])
    )
    v_vec[2] = (vel_vec[0] * np.cos(pos_vec[1])) - (
        pos_vec[0] * np.sin(pos_vec[1]) * vel_vec[1]
    )
    return v_vec
