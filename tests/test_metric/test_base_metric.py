import warnings

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy import constant

from einsteinpy.coordinates import SphericalDifferential, BoyerLindquistDifferential

from einsteinpy.metric import BaseMetric, Kerr, KerrNewman, Schwarzschild
from einsteinpy.utils.exceptions import CoordinateError

_c = constant.c.value


@pytest.fixture
def sph():
    return SphericalDifferential(
        e0=10000.0 * u.s,
        e1=130.0 * u.m,
        e2=np.pi / 2 * u.rad,
        e3=-np.pi / 8 * u.rad,
        u0=0.0 * u.m / u.s,
        u1=0.0 * u.rad / u.s,
        u2=0.0 * u.rad / u.s,
    )


@pytest.fixture
def bl():
    return BoyerLindquistDifferential(
        e0=10000.0 * u.s,
        e1=130.0 * u.m,
        e2=np.pi / 2 * u.rad,
        e3=-np.pi / 8 * u.rad,
        u0=0.0 * u.m / u.s,
        u1=0.0 * u.rad / u.s,
        u2=0.0 * u.rad / u.s,
    )


def test_str_repr(sph, bl):
    """
    Tests, if the ``__str__`` and ``__repr__`` messages match

    """
    M = 1e22 * u.kg
    a = 0.5 * u.one
    Q = 0. * u.C
    ms = Schwarzschild(coords=sph, M=M)
    mk = Kerr(coords=bl, M=M, a=a)
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)

    assert str(ms) == repr(ms)
    assert str(mk) == repr(mk)
    assert str(mkn) == repr(mkn)


@pytest.mark.parametrize(
    "M, a, Q, q",
    [
        (1e22 * u.m, 0.3 * u.m, 90, 200 * u.C),
        (1e22, 0.3, 78 * u.C, 20)
    ],
)
def test_wrong_or_no_units_in_init(M, a, Q, q):
    """
    Tests, if wrong or no units are flagged as error, while instantiation

    """
    sph = SphericalDifferential(
        e0=10000.0 * u.s,
        e1=130.0 * u.m,
        e2=np.pi / 2 * u.rad,
        e3=-np.pi / 8 * u.rad,
        u0=0.0 * u.m / u.s,
        u1=0.0 * u.rad / u.s,
        u2=0.0 * u.rad / u.s,
    )

    bl = BoyerLindquistDifferential(
        e0=10000.0 * u.s,
        e1=130.0 * u.m,
        e2=np.pi / 2 * u.rad,
        e3=-np.pi / 8 * u.rad,
        u0=0.0 * u.m / u.s,
        u1=0.0 * u.rad / u.s,
        u2=0.0 * u.rad / u.s,
    )

    try:
        bm = BaseMetric(coords=sph, M=M, a=a, Q=Q)
        ms = Schwarzschild(coords=sph, M=M)
        mk = Kerr(coords=bl, M=M, a=a)
        mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q, q=q)

        assert False

    except (u.UnitsError, TypeError):
        assert True


@pytest.fixture
def dummy_input(sph, bl):
    M = 5e27 * u.kg
    a = 0. * u.one
    Q = 0. * u.C

    ms = Schwarzschild(coords=bl, M=M)
    mk = Kerr(coords=sph, M=M, a=a)
    mkn = KerrNewman(coords=sph, M=M, a=a, Q=Q)

    x_vec_sph = sph.position()
    x_vec_bl = bl.position()

    return ms, mk, mkn, x_vec_sph, x_vec_bl


def test_coordinate_mismatch0(dummy_input):
    """
    Tests, if CoordinateError is raised, in case of coordinate system mismatch, \
    between metric coordinates and supplied coordinates

    """
    ms, mk, mkn, x_vec_sph, x_vec_bl = dummy_input

    def test_ms_cov(ms):
        ms_cov = ms.metric_covariant(x_vec_sph)

    def test_mk_cov(mk):
        mk_cov = mk.metric_covariant(x_vec_bl)

    def test_mkn_cov(mkn):
        mkn_cov = mkn.metric_covariant(x_vec_bl)

    def test_ms_con(ms):
        ms_cov = ms.metric_contravariant(x_vec_sph)

    with pytest.raises(CoordinateError):
        test_ms_cov(ms)

    with pytest.raises(CoordinateError):
        test_mk_cov(mk)

    with pytest.raises(CoordinateError):
        test_mkn_cov(mkn)

    with pytest.raises(CoordinateError):
        test_ms_con(ms)


def test_coordinate_mismatch1(dummy_input):
    """
    Tests, if CoordinateError is raised, in case of coordinate system mismatch, \
    between metric coordinates and supplied coordinates

    """
    ms, mk, mkn, x_vec_sph, x_vec_bl = dummy_input

    def test_ms_chl(ms):
        ms_chl = ms.christoffels(x_vec_sph)

    def test_mk_chl(mk):
        mk_chl = mk.christoffels(x_vec_bl)

    def test_mkn_chl(mkn):
        mkn_chl = mkn.christoffels(x_vec_bl)

    def test_ms_fvec(ms):
        ms_fvec = ms.f_vec(0., x_vec_sph)

    with pytest.raises(CoordinateError):
        test_ms_chl(ms)

    with pytest.raises(CoordinateError):
        test_mk_chl(mk)

    with pytest.raises(CoordinateError):
        test_mkn_chl(mkn)

    with pytest.raises(CoordinateError):
        test_ms_fvec(ms)


def test_coordinate_mismatch2(dummy_input):
    """
    Tests, if CoordinateError is raised, in case of coordinate system mismatch, \
    between metric coordinates and supplied coordinates

    """
    ms, mk, mkn, x_vec_sph, x_vec_bl = dummy_input

    def test_mk_con(mk):
        mk_cov = mk.metric_contravariant(x_vec_bl)

    def test_mkn_con(mkn):
        mkn_cov = mkn.metric_contravariant(x_vec_bl)

    def test_mk_fvec(mk):
        ms_fvec = mk.f_vec(0., x_vec_bl)

    def test_mkn_fvec(mkn):
        ms_fvec = mkn.f_vec(0., x_vec_bl)

    with pytest.raises(CoordinateError):
        test_mk_con(mk)

    with pytest.raises(CoordinateError):
        test_mkn_con(mkn)

    with pytest.raises(CoordinateError):
        test_mk_fvec(mk)

    with pytest.raises(CoordinateError):
        test_mkn_fvec(mkn)


def dummy_met(x_vec):
    """
    Dummy Metric function for ``test_perturbation()``

    """
    # This isn't a physical metric
    # Just checking the matrix addition
    return 1.5 * np.ones((4, 4), dtype=float)


def dummy_perturb(x_vec):
    """
    Dummy Perturbation function for ``test_perturbation()``

    """
    # This isn't a physical metric/perturbation
    # Just checking the matrix addition
    return 1.9 * np.ones((4, 4), dtype=float)


def test_perturbation(sph):
    """
    Tests, if the ``pertubation`` is added correctly

    """
    met = BaseMetric(
        coords=sph,
        M=1e22 * u.kg,
        a=0.75 * u.one,
        Q=1. * u.C,
        metric_cov=dummy_met,
        perturbation=dummy_perturb
    )

    met_calc = 3.4 * np.ones((4, 4), dtype=float)

    assert_allclose(met.metric_covariant(np.ones(4, dtype=float)), met_calc, rtol=1e-10)


def test_unperturbed_metric_covariant(sph):
    """
    Tests, if unperturbed metric is returned

    """
    met = BaseMetric(coords=sph, M=1e22 * u.kg, a=0.75 * u.one, Q=1. * u.C, metric_cov=dummy_met, perturbation=None)

    met_calc = dummy_met(np.ones(4, dtype=float))

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


def test_compare_metrics_under_limits(sph, bl):
    """
    Tests if KerrNewman Metric reduces to Kerr Metric, in the limit Q -> 0 and to Schwarzschild \
    Metric, in the limits, a -> 0 & Q -> 0

    """
    M = 6.73317655e26 * u.kg
    a1, a2 = 0.5 * u.one, 0. * u.one
    Q = 0. * u.C

    ms = Schwarzschild(coords=sph, M=M)
    mk = Kerr(coords=bl, M=M, a=a1)
    mk0 = Kerr(coords=bl, M=M, a=a2)
    mkn = KerrNewman(coords=bl, M=M, a=a1, Q=Q)
    mkn0 = KerrNewman(coords=bl, M=M, a=a2, Q=Q)

    x_vec_sph = sph.position()
    x_vec_bl = bl.position()

    ms_mat = ms.metric_covariant(x_vec_sph)
    mk_mat = mk.metric_covariant(x_vec_bl)
    mk0_mat = mk0.metric_covariant(x_vec_bl)
    mkn_mat = mkn.metric_covariant(x_vec_bl)
    mkn0_mat = mkn0.metric_covariant(x_vec_bl)

    assert_allclose(ms_mat, mk0_mat, rtol=1e-8)
    assert_allclose(mk_mat, mkn_mat, rtol=1e-8)
    assert_allclose(mkn0_mat, ms_mat, rtol=1e-8)


def test_compare_kerr_kerrnewman_metric(bl):
    """
    Tests, if covariant & contravariant forms of Kerr & Kerr-Newman metrics match, when Q -> 0

    """
    M = 1e23 * u.kg
    a = 0.99 * u.one
    Q = 0. * u.C

    x_vec = bl.position()

    mk = Kerr(coords=bl, M=M, a=a)
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
    mk_contra = mk.metric_contravariant(x_vec)
    mkn_contra = mkn.metric_contravariant(x_vec)
    mk_cov = mk.metric_covariant(x_vec)
    mkn_cov = mkn.metric_covariant(x_vec)

    assert_allclose(mk_contra, mkn_contra, rtol=1e-10)
    assert_allclose(mk_cov, mkn_cov, rtol=1e-10)


def test_singularities_for_uncharged_nonrotating_case(sph, bl):
    """
    Tests, if all metric singularities match up across Schwarzschild, Kerr & Kerr-Newman \
    spacetimes, when a -> 0 & Q -> 0, subject to choice to coordinates (here, Spherical and BL)

    """
    M = 5e27 * u.kg
    a = 0. * u.one
    Q = 0. * u.C

    theta = np.pi / 4

    ms = Schwarzschild(coords=sph, M=M)
    mk = Kerr(coords=bl, M=M, a=a)
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)
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

    scr = ms.schwarzschild_radius(M.value)

    assert_allclose(mssinglist, mksinglist, rtol=1e-4, atol=0.0)
    assert_allclose(mksinglist, mknsinglist, rtol=1e-4, atol=0.0)
    assert_allclose(mknsinglist, mssinglist, rtol=1e-4, atol=0.0)
    assert_allclose(mksinglist[2], scr, rtol=1e-4, atol=0.0)


def test_singularities_raises_CoordinateError(sph, bl):
    """
    Tests, if ``singularities()`` raises a CoordinateError, \
    when there is a coordinate mismatch between supplied coordinate \
    object and the coordinate object, metric was instantiated with

    """
    M = 5e27 * u.kg
    a = 0. * u.one
    Q = 0. * u.C

    theta = np.pi / 4

    ms = Schwarzschild(coords=bl, M=M)
    mk = Kerr(coords=sph, M=M, a=a)
    mkn = KerrNewman(coords=sph, M=M, a=a, Q=Q)

    def mssing(ms):
        mssing = ms.singularities()

    def mksing(mk):
        mksing = mk.singularities()

    def mknsing(mkn):
        mknsing = mkn.singularities()

    with pytest.raises(CoordinateError):
        mssing(ms)

    with pytest.raises(CoordinateError):
        mksing(mk)

    with pytest.raises(CoordinateError):
        mknsing(mkn)


def test_deprecation_warning_for_calculate_trajectory(sph, bl):
    """
    Tests, if a Deprecation Warning is shown, when accessing calculate_trajectory \
    for all metric classes

    """
    M, a, Q = 5e27 * u.kg, 0. * u.one, 0. * u.C
    ms = Schwarzschild(coords=sph, M=M)
    mk = Kerr(coords=bl, M=M, a=a)
    mkn = KerrNewman(coords=bl, M=M, a=a, Q=Q)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        ms.calculate_trajectory()
        mk.calculate_trajectory()
        mkn.calculate_trajectory()

        assert len(w) == 3  # 3 warnings to be shown
        assert issubclass(w[-1].category, DeprecationWarning)
