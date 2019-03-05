import astropy.units as u
import numpy as np
from astropy import constants

h = constants.h  # planck's constant
h_red = constants.h / (2 * np.pi)  # Reduced planck's constant
c = constants.c  # speed of light
G = constants.G  # Gravitational Constant

Cosmo_Const_base = constants.Constant(
    "lambda",
    "Cosmological Constant",
    2.036e-35,
    "1 / s2",
    0.000081e-35,
    "Wikipedia",
    system="si",
)

Cosmo_Const = Cosmo_Const_base.value * (u.s ** -2)
