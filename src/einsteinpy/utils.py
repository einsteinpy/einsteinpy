import astropy.units as u

from . import constant


@u.quantity_input(mass=u.kg)
def schwarzschild_radius(mass):
    """
    Schwarzschild radius
    """
    num = 2 * constant.G * mass
    deno = constant.c**2
    return num/deno
