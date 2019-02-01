from astropy import constants
from astropy import units as u

c = constants.c.value * (u.m / u.s)

G = constants.G.value * ((u.m ** 3) / (u.kg * u.s**2))

Cosmo_Const = constants.Constant(
    'lambda',
    'Cosmological Constant',
    2.036e-35,
    '1 / s2',
    0.000081e-35,
    "Wikipedia",
    system='si'
).value * (u.s ** -2)
