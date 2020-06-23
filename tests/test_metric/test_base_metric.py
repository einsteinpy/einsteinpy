import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy import constant

from einsteinpy.metric import BaseMetric, Kerr, KerrNewman, Schwarzschild

_c = constant.c.value


def test_str_repr():
    """
    Tests, if the ``__str__`` and ``__repr__`` messages match

    """
    ms = Schwarzschild(M=1e22)
    mk = Kerr(coords="BL", M=1e22, a=0.5)
    mkn = KerrNewman(coords="BL", M=1e22, a=0.5, Q=0.)

    assert str(ms) == repr(ms)
    assert str(mk) == repr(mk)
    assert str(mkn) == repr(mkn)


def test_r_ks_raises_NotImplementedError():
    """
    Tests, if ``r_ks()`` raises a NotImplementedError, when invoked

    """
    try:
        BaseMetric.r_ks(1., 1., 1., 0.8)
        assert False

    except NotImplementedError:
        assert True


def dummy_met(x_vec):
    """
    Dummy Metric function for ``test_perturbation()``

    """
    return 1.5 * np.ones((4, 4), dtype=float)


def dummy_perturb(x_vec):
    """
    Dummy Perturbation function for ``test_perturbation()``

    """
    return 1.9 * np.ones((4, 4), dtype=float)


def test_perturbation():
    """
    Tests, if the ``pertubation`` is added correctly

    """
    met = BaseMetric(coords="BL", M=1e22, a=0.75, Q=1., metric_cov=dummy_met, perturbation=dummy_perturb)

    met_calc = 3.4 * np.ones((4, 4), dtype=float)

    assert_allclose(met.metric_covariant(np.ones(4, dtype=float)), met_calc, rtol=1e-10)


@pytest.mark.parametrize("a", [1.1, -0.1])
def test_alpha_raises_error(a):
    """
    Tests, if AssertionError is raised for invalid 'a' inputs

    """
    try:
        BaseMetric.alpha(a, 3e20)
        assert False
    except ValueError:
        assert True


def test_compare_metrics_under_limits():
    """
    Tests if KerrNewman Metric reduces to Kerr Metric, in the limit Q -> 0 and to Schwarzschild \
    Metric, in the limits, a -> 0 & Q -> 0

    """
    r, theta = 99.9, 5 * np.pi / 6
    M = 6.73317655e26

    x_vec = np.array([0., r, theta, 0.])
    ms = Schwarzschild(M=M)
    mk = Kerr(coords="BL", M=M, a=0.5)
    mk0 = Kerr(coords="BL", M=M, a=0.)
    mkn = KerrNewman(coords="BL", M=M, a=0.5, Q=0.)
    mkn0 = KerrNewman(coords="BL", M=M, a=0., Q=0.)

    ms_mat = ms.metric_covariant(x_vec)
    mk_mat = mk.metric_covariant(x_vec)
    mk0_mat = mk0.metric_covariant(x_vec)
    mkn_mat = mkn.metric_covariant(x_vec)
    mkn0_mat = mkn0.metric_covariant(x_vec)

    assert_allclose(ms_mat, mk0_mat, rtol=1e-8)
    assert_allclose(mk_mat, mkn_mat, rtol=1e-8)
    assert_allclose(mkn0_mat, ms_mat, rtol=1e-8)


def test_compare_kerr_kerrnewman_metric():
    """
    Tests, if covariant & contravariant forms of Kerr & Kerr-Newman metrics match, when Q -> 0

    """
    r, theta, M, a = 0.1, 4 * np.pi / 5, 1e23, 0.99
    x_vec = np.array([0., r, theta, 0.])

    mk = Kerr(coords="BL", M=M, a=a)
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=0.)
    mk_contra = mk.metric_contravariant(x_vec)
    mkn_contra = mkn.metric_contravariant(x_vec)
    mk_cov = mk.metric_covariant(x_vec)
    mkn_cov = mkn.metric_covariant(x_vec)

    assert_allclose(mk_contra, mkn_contra, rtol=1e-10)
    assert_allclose(mk_cov, mkn_cov, rtol=1e-10)


def test_singularities_for_uncharged_nonrotating_case():
    """
    Tests, if all metric singularities match up across Schwarzschild, Kerr & Kerr-Newman \
    spacetimes, when a -> 0 & Q -> 0, subject to choice to coordinates (here, Spherical and BL)

    """
    theta = np.pi / 4
    M, a, Q = 5e27, 0., 0.

    ms = Schwarzschild(M=M)
    mk = Kerr(coords="BL", M=M, a=a)
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q)
    mssing = ms.singularities()
    mksing = mk.singularities()
    mknsing = mkn.singularities()

    mssinglist = [
        mssing["inner_ergosphere"],
        mssing["inner_horizon"],
        mssing["outer_horizon"],
        mssing["outer_ergosphere"]
    ]

    mksinglist = [
        mksing["inner_ergosphere"](theta),
        mksing["inner_horizon"],
        mksing["outer_horizon"],
        mksing["outer_ergosphere"](theta)
    ]

    mknsinglist = [
        mknsing["inner_ergosphere"](theta),
        mknsing["inner_horizon"],
        mknsing["outer_horizon"],
        mknsing["outer_ergosphere"](theta)
    ]

    scr = ms.schwarzschild_radius(M)

    assert_allclose(mssinglist, mksinglist, rtol=1e-4, atol=0.0)
    assert_allclose(mksinglist, mknsinglist, rtol=1e-4, atol=0.0)
    assert_allclose(mknsinglist, mssinglist, rtol=1e-4, atol=0.0)
    assert_allclose(mksinglist[2], scr, rtol=1e-4, atol=0.0)


def test_deprecation_warning_for_calculate_trajectory():
    """
    Tests, if a Deprecation Warning is shown, when accessing calculate_trajectory \
    for all metric classes
    """
    M, a, Q = 5e27, 0., 0.
    ms = Schwarzschild(M=M)
    mk = Kerr(coords="BL", M=M, a=a)
    mkn = KerrNewman(coords="BL", M=M, a=a, Q=Q)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        ms.calculate_trajectory()
        mk.calculate_trajectory()
        mkn.calculate_trajectory()

        assert len(w) == 3  # 3 warnings to be shown
        assert issubclass(w[-1].category, DeprecationWarning)
