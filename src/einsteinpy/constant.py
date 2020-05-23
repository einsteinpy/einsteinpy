import astropy.units as u
import numpy as np
from astropy import constants

c = constants.c
G = constants.G
eps0 = constants.eps0
coulombs_const = 1 / (4 * np.pi * eps0)

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

Solar_Mass = 1.9891e30 * u.kg
R_sun = 695510 * u.km
