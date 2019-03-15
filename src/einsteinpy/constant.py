import astropy.units as u
from astropy import constants

c = constants.c.value * (u.m / u.s)
G = constants.G.value * ((u.m ** 3) / (u.kg * u.s ** 2))
pi_by_2 = 1.5707963267948965579989817342720925807952880859375

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
