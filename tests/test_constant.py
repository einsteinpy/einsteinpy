import astropy.units as u
import numpy as np
import pytest
from astropy import constants
from numpy.testing import assert_allclose

from einsteinpy.constant import Cosmo_Const, Cosmo_Const_base, R_sun, Solar_Mass


def test_Cosmo_Const_returns_correct_value_units():
    cnst = Cosmo_Const
    assert_allclose(cnst.value, 9.9359e-36, rtol=1e-4)
    assert isinstance(cnst.unit, u.core.CompositeUnit)
    assert isinstance(cnst, u.quantity.Quantity)


def test_Cosmo_Const_matches_Planck_2018_in_inverse_square_metres():
    # Planck 2018 quotes Lambda = 1.1056e-52 m^-2, which is Lambda / c^2 here
    assert u.allclose(Cosmo_Const / constants.c**2, 1.1056e-52 * u.m**-2, rtol=1e-3)


def test_Cosmo_Const_gives_Planck_2018_de_Sitter_Hubble_rate():
    # In the dark energy era H = sqrt(Lambda / 3) = H0 sqrt(Omega_Lambda)
    H = np.sqrt(Cosmo_Const / 3)
    H_expected = 67.66 * np.sqrt(0.6889) * u.km / u.s / u.Mpc
    assert u.allclose(H, H_expected, rtol=1e-3)


def test_Cosmo_const_is_astropy_constant_with_unit_and_uncert():
    cnst = Cosmo_Const_base
    assert isinstance(cnst, constants.Constant)
    assert_allclose(cnst.uncertainty / cnst.value, 0.0148, rtol=0.01)
    assert cnst.unit == u.Unit("1 / s2")


def test_Cosmo_Const_has_correct_metadata():
    cnst = Cosmo_Const_base
    assert cnst.name == "Cosmological Constant"
    assert cnst.system == "si" and cnst.abbrev == "lambda"
    assert cnst.reference == "Planck 2018 (astropy.cosmology.Planck18)"


def test_solar_values_match_IAU_2015_nominal():
    # IAU 2015 Resolution B3: GM_sun = 1.3271244e20 m^3 / s^2, R_sun = 695700 km
    rs = 2 * constants.G * Solar_Mass / constants.c**2
    assert u.allclose(
        rs, 2 * 1.3271244e20 * u.m**3 / u.s**2 / constants.c**2, rtol=1e-7
    )
    assert u.allclose(R_sun, 695700 * u.km)
    assert Solar_Mass.unit == u.kg and R_sun.unit == u.km
