import astropy.units as u
import numpy as np
import coord_transforms as sph
import bl_coord_transforms as bl


def BLToSpherical_pos(pos_vec, a):
    """
    Function to convert Boyer-Lindquist coordinates to Spherical coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having r, theta, phi coordinates in SI units(m, rad, rad)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array with r, theta, phi (m, rad, rad)

    """
    r_vec = bl.BLToCartesian_pos(pos_vec, a)
    r_vec = sph.CartesianToSpherical_pos(r_vec)

    return r_vec

def BLToSpherical_vel(pos_vec, vel_vec, a):
    """
    Function to convert velocities in Boyer-Lindquist coordinates to velocities in Spherical coordinates

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
        3-length numpy array with V_r, V_theta, V_phi (m/s, rad/s, rad/s)

    """
    r_vec = bl.BLToCartesian_pos(pos_vec, a)

    v_vec = bl.BLToCartesian_vel(pos_vec, vel_vec, a)
    v_vec = sph.CartesianToSpherical_vel(r_vec, v_vec)

    return v_vec


def SphericalToBL_pos(pos_vec, a):
    """
    Function to convert Spherical coordinates to Boyer-Lindquist coordinates

    Parameters
    ----------
    pos_vec : ~numpy.array
        3-length numpy array having r, theta, phi coordinates in SI units(m, rad, rad)
    a : float
        Any constant

    Returns
    -------
    ~numpy.array
        3-length numpy array with r, theta, phi in (m, rad, rad)

    """
    r_vec = sph.SphericalToCartesian_pos(pos_vec)
    r_vec = bl.CartesianToBL_pos(r_vec, a)

    return r_vec

def SphericalToBL_vel(pos_vec, vel_vec, a):
    """
    Function to convert velocities in Spherical coordinates to velocities in Boyer-Lindquist coordinates

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
        3-length numpy array having v_r, v_theta, v_phi in (m/s, rad/s, rad/s)

    """
    r_vec = sph.SphericalToCartesian_pos(pos_vec)

    v_vec = sph.SphericalToCartesian_vel(pos_vec, vel_vec)
    v_vec = bl.CartesianToBL_vel(r_vec, v_vec, a)

    return v_vec

def S2BL_units(pos_vec, vel_vec, a):
    """
    Function to convert Spherical to Boyer-Lindquist Coordinates along with handling units

    Parameters
    ----------
    pos_vec : list
        list of 3 position coordinates along with appropriate units
        [r, theta, phi]
        (u.m, u.rad, u.rad)
    vel_vec : list
        list of 3 velocity coordinates along with appropriate units
        [vr, vtheta, vphi]
        (u.m/u.s, u.rad/u.s, u.rad/u.s)
    a : float
        Any constant

    Returns
    -------
    tuple
        consisting of 2 lists
        ([r, theta, phi], [vr, vtheta, vphi]) in units
        ([u.m, u.rad, u.rad], [u.m/u.s, u.rad/u.s, u.rad/u.s])

    """
    npos_vec, nvel_vec = sph.S2C_units(pos_vec, vel_vec)
    npos_vec, nvel_vec = bl.C2BL_units(npos_vec, nvel_vec, a)

    return (npos_vec, nvel_vec)

def BL2S_units(pos_vec, vel_vec, a):
    """
    Function to convert Boyer-Lindquist to Spherical Coordinates along with handling units

    Parameters
    ----------
    pos_vec : list
        list of 3 position coordinates along with appropriate units
        [r, theta, phi]
        (u.m, u.rad, u.rad)
    vel_vec : list
        list of 3 velocity coordinates along with appropriate units
        [vr, vtheta, vphi]
        (u.m/u.s, u.rad/u.s, u.rad/u.s)
    a : float
        Any constant

    Returns
    -------
    tuple
        consisting of 2 lists
        ([r, theta, phi], [vr, vtheta, vphi]) in units
        ([u.m, u.rad, u.rad], [u.m/u.s, u.rad/u.s, u.rad/u.s])

    """
    npos_vec, nvel_vec = bl.BL2C_units(pos_vec, vel_vec, a)
    npos_vec, nvel_vec = sph.C2S_units(npos_vec, nvel_vec)

    return (npos_vec, nvel_vec)