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
            1.0000000005743104,
            1.819879425905945e-18,
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


def test_scalar_factor_dark_energy_era_grows_at_de_Sitter_rate():
    # a'/a in the dark energy era is H = sqrt(Lambda / 3), about 56 km/s/Mpc
    # for Planck 2018. The parametrized test above cannot see this, because
    # its atol is much larger than a'.
    t = 10 * u.yr
    H = scalar_factor_derivative(t, "ded") / scalar_factor(t, "ded")
    H_expected = (67.66 * 0.6889**0.5 * u.km / u.s / u.Mpc).to_value(u.s**-1)
    assert_allclose(H, H_expected, rtol=1e-3, atol=0)
