"""Example data.
"""
import astropy.units as u
import numpy as np

from einsteinpy.bodies import Body
from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.geodesic import Geodesic


def perihelion():
    Attractor = Body(name="BH", mass=6e24 * u.kg, parent=None)
    sph_obj = SphericalDifferential(
        130 * u.m,
        np.pi / 2 * u.rad,
        -np.pi / 8 * u.rad,
        0 * u.m / u.s,
        0 * u.rad / u.s,
        1900 * u.rad / u.s,
    )
    Object = Body(differential=sph_obj, parent=Attractor)
    geodesic = Geodesic(body=Object, time=0 * u.s, end_lambda=0.002, step_size=5e-8)
    return geodesic
