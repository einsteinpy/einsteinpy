import astropy.units as u

from einsteinpy import constant


@u.quantity_input(mass=u.kg)
def schwarzschild_radius(mass):
    """
    Schwarzschild radius is the radius defining the event horizon of a
    Schwarzschild black hole. It is characteristic radius associated with every
    quantity of mass.

    Parameters
    ----------
    mass : float
    """
    num = 2 * constant.G * mass
    deno = constant.c ** 2
    return num / deno
