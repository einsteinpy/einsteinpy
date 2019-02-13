import astropy.units as u
from astropy import constants

c = constants.c.value * (u.m / u.s)

G = constants.G.value * ((u.m ** 3) / (u.kg * u.s ** 2))

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
