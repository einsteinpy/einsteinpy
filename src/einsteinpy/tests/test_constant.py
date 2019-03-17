import astropy.units as u
import numpy as np
import pytest
from astropy import constants
from numpy.testing import assert_allclose

from einsteinpy.constant import Cosmo_Const, Cosmo_Const_base


def test_Cosmo_Const_returns_correct_value_units():
    cnst = Cosmo_Const
    assert cnst.value == 2.036e-35
    assert isinstance(cnst.unit, u.core.CompositeUnit)
    assert isinstance(cnst, u.quantity.Quantity)


def test_Cosmo_const_is_astropy_constant_with_unit_and_uncert():
    cnst = Cosmo_Const_base
    assert isinstance(cnst, constants.Constant)
    assert cnst.uncertainty == 8.1e-40
    assert cnst.unit == u.Unit("1 / s2")


def test_Cosmo_Const_has_correct_metadata():
    cnst = Cosmo_Const_base
    assert cnst.name == "Cosmological Constant"
    assert cnst.system == "si" and cnst.abbrev == "lambda"
    assert cnst.reference == "Wikipedia"
