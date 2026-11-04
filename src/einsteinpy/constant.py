import numpy as np
from astropy import constants, units as u
from astropy.cosmology import Planck18

__all__ = [
    "c",
    "G",
    "eps0",
    "coulombs_const",
    "Cosmo_Const",
    "Solar_Mass",
    "R_sun",
]

c = constants.c
G = constants.G
eps0 = constants.eps0
coulombs_const = 1 / (4 * np.pi * eps0)

# Lambda c^2 = 3 H0^2 Omega_Lambda, from the Planck 2018 cosmology in astropy.
# The uncertainty comes from H0 = 67.66 +/- 0.42 km/s/Mpc and
# Omega_Lambda = 0.6889 +/- 0.0056 (Planck 2018 VI, Table 2, last column),
# added in quadrature, which ignores their correlation.
_Lambda = (3 * Planck18.H0**2 * Planck18.Ode0).to_value(u.s**-2)
_Lambda_rel_err = np.hypot(2 * 0.42 / 67.66, 0.0056 / 0.6889)

Cosmo_Const_base = constants.Constant(
    "lambda",
    "Cosmological Constant",
    _Lambda,
    "1 / s2",
    _Lambda * _Lambda_rel_err,
    "Planck 2018 (astropy.cosmology.Planck18)",
    system="si",
)

Cosmo_Const = Cosmo_Const_base.value * (u.s**-2)

# IAU 2015 nominal values, through astropy
Solar_Mass = constants.M_sun.to(u.kg)
R_sun = constants.R_sun.to(u.km)
