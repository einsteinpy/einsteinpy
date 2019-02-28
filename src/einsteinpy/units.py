import astropy.units as u

from einsteinpy.constant import G, c


def astro_dist(mass):
    """
    Function for turning distance into astronomical perspective.
    """
    value = (G.value * mass) / (c.value ** 2)
    astro_dist = u.def_unit("astro-m(GM/c2)", u.m / value)
    return astro_dist


def astro_sec(mass):
    """
    Function for turning time into astronomical perspective.
    """
    value = (G.value * mass) / (c.value ** 3)
    astro_sec = u.def_unit("astro-sec(GM/c3)", u.s / value)
    return astro_sec
