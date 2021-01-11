import pytest

from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.utils import scalar_factor, scalar_factor_derivative


@pytest.mark.parametrize(
    "t, era, tuning_param, sf, sfd",
    [
        (
            10 * u.s,
            "md",
            1.2,
            5.569906600335334,
            0.3713271066890223,
        ),
        (
            10 * u.h,
            "rd",
            2.,
            379.4733192202055,
            0.005270462766947299,
        ),
        (
            10 * u.yr,
            "ded",
            1.,
            1.0000000008221144,
            2.6051231598190325e-18,
        ),
    ],
)
def test_scalar_factor_and_derivative(t, era, tuning_param, sf, sfd):
    assert_allclose(scalar_factor(t, era, tuning_param), sf, rtol=1e-6, atol=1e-6)
    assert_allclose(scalar_factor_derivative(t, era, tuning_param), sfd, rtol=1e-6, atol=1e-6)


def test_scalar_factor_and_derivative_raise_ValueError():
    with pytest.raises(ValueError):
        sf = scalar_factor(1 * u.s, "era")
    
    with pytest.raises(ValueError):
        sf = scalar_factor_derivative(1 * u.s, "era")
