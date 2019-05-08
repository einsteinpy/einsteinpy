import astropy.units as u
import numpy as np

from einsteinpy import constant


@u.quantity_input(mass=u.kg)
def schwarzschild_radius(mass):
    """
    Schwarzschild radius is the radius defining the event horizon of a
    Schwarzschild black hole. It is characteristic radius associated with every
    quantity of mass.

    Parameters
    ----------
    mass : ~astropy.units.kg

    Returns
    -------
    r : ~astropy.units.m
        Schwarzschild radius for a given mass

    """
    M = mass.to(u.kg)
    num = 2 * constant.G * M
    deno = constant.c ** 2
    return num / deno


@u.quantity_input(mass=u.kg)
def time_velocity(pos_vec, vel_vec, mass):
    """
    Velocity of time calculated from einstein's equation.
    See http://www.physics.usu.edu/Wheeler/GenRel/Lectures/GRNotesDecSchwarzschildGeodesicsPost.pdf

    Parameters
    ----------
    pos_vector : ~numpy.array
        Vector with r, theta, phi components in SI units
    vel_vector : ~numpy.array
        Vector with velocities of r, theta, phi components in SI units
    mass : ~astropy.units.kg
        Mass of the body

    Returns
    -------
    ~astropy.units.one
        Velocity of time

    """
    # this function considers SI units only
    a = schwarzschild_radius(mass).value
    c = constant.c.value
    num1 = (1 / ((c ** 2) * (1 + (a / pos_vec[0])))) * (vel_vec[0] ** 2)
    num2 = ((pos_vec[0] ** 2) / (c ** 2)) * (vel_vec[1] ** 2)
    num3 = (
        ((pos_vec[0] ** 2) / (c ** 2)) * (np.sin(pos_vec[1]) ** 2) * (vel_vec[2] ** 2)
    )
    deno = 1 + (a / pos_vec[0])
    time_vel_squared = (1 + num1 + num2 + num3) / deno
    time_vel = np.sqrt(time_vel_squared)
    return time_vel * u.one


def metric(c, r, theta, Rs):
    """
    Returns the Schwarzschild Metric

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

    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4)

    """
    m = np.zeros(shape=(4, 4), dtype=float)
    tmp, c2 = 1.0 - (Rs / r), c ** 2
    m[0, 0] = tmp
    m[1, 1] = -1.0 / (tmp * c2)
    m[2, 2] = -1 * (r ** 2) / c2
    m[3, 3] = -1 * ((r * np.sin(theta)) ** 2) / c2
    return m


def christoffels(c, r, theta, Rs):
    """
    Returns the 3rd rank Tensor containing Christoffel Symbols for Schwarzschild Metric

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

    Returns
    -------
    ~numpy.array
        Numpy array of shape (4,4,4)

    """
    chl = np.zeros(shape=(4, 4, 4), dtype=float)
    c2 = c ** 2
    chl[1, 0, 0] = 0.5 * Rs * (r - Rs) * c2 / (r ** 3)
    chl[1, 1, 1] = 0.5 * Rs / (Rs * r - r ** 2)
    chl[1, 2, 2] = Rs - r
    chl[1, 3, 3] = (Rs - r) * (np.sin(theta) ** 2)
    chl[0, 0, 1] = chl[0, 1, 0] = -chl[1, 1, 1]
    chl[2, 2, 1] = chl[2, 1, 2] = chl[3, 3, 1] = chl[3, 1, 3] = 1 / r
    chl[2, 3, 3] = -np.cos(theta) * np.sin(theta)
    chl[3, 3, 2] = chl[3, 2, 3] = 1 / np.tan(theta)
    return chl
