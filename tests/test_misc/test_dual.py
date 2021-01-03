import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.misc.dual import Dual, _deriv, _diff_g, _jacobian_g
from einsteinpy.geodesic.utils import _sch, _kerr, _kerrnewman


def test_str_repr_Dual():
    x = Dual(1., 2.)

    assert str(x) == repr(x) == "Dual(1.0, 2.0)"


def test_Dual1():
    x = Dual(1., 2.)
    y = Dual(2., 3.)
    c = 2.

    add_ = x + y
    assert_allclose([add_.a, add_.b], [3., 5.], atol=1e-8, rtol=1e-8)

    radd_ = x.__radd__(c)
    assert_allclose([radd_.a, radd_.b], [3., 2.], atol=1e-8, rtol=1e-8)

    sub_ = x - y
    assert_allclose([sub_.a, sub_.b], [-1., -1.], atol=1e-8, rtol=1e-8)

    rsub_ = x.__rsub__(c)
    assert_allclose([rsub_.a, rsub_.b], [1., -2.], atol=1e-8, rtol=1e-8)

    mul_ = x * y
    assert_allclose([mul_.a, mul_.b], [2., 7.], atol=1e-8, rtol=1e-8)

    rmul_ = x.__rmul__(c)
    assert_allclose([rmul_.a, rmul_.b], [2., 4.], atol=1e-8, rtol=1e-8)


def test_Dual2():
    x = Dual(1., 2.)
    y = Dual(2., 3.)
    c = 2.

    add_ = x + c
    assert_allclose([add_.a, add_.b], [3., 2.], atol=1e-8, rtol=1e-8)

    radd_ = x.__radd__(y)
    assert_allclose([radd_.a, radd_.b], [3., 5.], atol=1e-8, rtol=1e-8)

    sub_ = x - c
    assert_allclose([sub_.a, sub_.b], [-1., 2.], atol=1e-8, rtol=1e-8)

    rsub_ = x.__rsub__(y)
    assert_allclose([rsub_.a, rsub_.b], [1., 1.], atol=1e-8, rtol=1e-8)

    mul_ = x * c
    assert_allclose([mul_.a, mul_.b], [2., 4.], atol=1e-8, rtol=1e-8)

    rmul_ = x.__rmul__(y)
    assert_allclose([rmul_.a, rmul_.b], [2., 7.], atol=1e-8, rtol=1e-8)


def test_Dual3():
    x = Dual(1., 2.)
    y = Dual(2., 3.)
    c = 2.

    truediv_ = x / y
    assert_allclose([truediv_.a, truediv_.b], [.5, .25], atol=1e-8, rtol=1e-8)

    rtruediv1_ = x / c
    assert_allclose([rtruediv1_.a, rtruediv1_.b], [.5, 1.], atol=1e-8, rtol=1e-8)

    rtruediv2_ = c / x
    assert_allclose([rtruediv2_.a, rtruediv2_.b], [2., -4.], atol=1e-8, rtol=1e-8)

    assert x == x

    assert x != y

    pow_ = y ** 2
    assert_allclose([pow_.a, pow_.b], [4., 12.], atol=1e-8, rtol=1e-8)

    neg_ = -y
    assert_allclose([neg_.a, neg_.b], [-2., -3.], atol=1e-8, rtol=1e-8)


def test_Dual4():
    x = Dual(1., 2.)
    y = Dual(2., 3.)
    c = 2.

    truediv_ = x / c
    assert_allclose([truediv_.a, truediv_.b], [.5, 1.], atol=1e-8, rtol=1e-8)

    rtruediv1_ = x.__rtruediv__(y)
    assert_allclose([rtruediv1_.a, rtruediv1_.b], [2., -1.], atol=1e-8, rtol=1e-8)

    x = Dual(0., 2.)
    y = Dual(0., 4.)

    truediv_ = x / y
    assert_allclose([truediv_.a, truediv_.b], [.5, 0.], atol=1e-8, rtol=1e-8)

    rtruediv1_ = x.__rtruediv__(y)
    assert_allclose([rtruediv1_.a, rtruediv1_.b], [2., 0.], atol=1e-8, rtol=1e-8)


def test_Dual5():
    x = Dual(1., 2.)
    y = Dual(2., 3.)
    c = 2.

    sin = x.sin()
    assert_allclose([sin.a, sin.b], [0.8414709848078965, 1.0806046117362795], atol=1e-8, rtol=1e-8)

    cos = y.cos()
    assert_allclose([cos.a, cos.b], [-0.4161468365471424, -2.727892280477045], atol=1e-8, rtol=1e-8)

    tan = x.tan()
    assert_allclose([tan.a, tan.b], [1.557407724654902, 6.851037641629518], atol=1e-8, rtol=1e-8)

    log = y.log()
    assert_allclose([log.a, log.b], [0.6931471805599453, 1.5], atol=1e-8, rtol=1e-8)

    exp = x.exp()
    assert_allclose([exp.a, exp.b], [2.718281828459045, 5.43656365691809], atol=1e-8, rtol=1e-8)


def test_deriv1():
    df = _deriv(lambda y: np.cos(y) * np.sin(y) ** 2, 1.)

    assert_allclose(df, -0.10452774015707361, atol=1e-8, rtol=1e-8)


def test_deriv2():
    df = _deriv(lambda y: 0., 0.)

    assert_allclose(df, 0., atol=1e-8, rtol=1e-8)


@pytest.mark.parametrize(
    "g, g_prms, coords, indices, wrt, expected",
    [
        (
            _sch,
            (),
            [0., 2.5, np.pi / 6, np.pi / 2],
            (2, 2),
            2,
            0.
        ),
        (
            _sch,
            (),
            [0., 2.5, np.pi / 6, np.pi / 2],
            (2, 2),
            1,
            -0.128
        ),
        (
            _kerr,
            (0.9,),
            [0., 25, np.pi / 2, 0.],
            (0, 3),
            1,
            1.542519202468211e-05
        ),
        (
            _kerr,
            (0.9,),
            [0., 25, np.pi / 2, 0.],
            (0, 3),
            3,
            0.
        ),
        (
            _kerrnewman,
            (0.5, 0.1,),
            [0., 5.5, np.pi / 4, 0.],
            (3, 0),
            0,
            0.
        ),
        (
            _kerrnewman,
            (0.5, 0.1,),
            [0., 5.5, np.pi / 4, 0.],
            (3, 0),
            2,
            -7.631639775075298e-05
        ),
    ],
)
def test_diff_g(g, g_prms, coords, indices, wrt, expected):
    diff = _diff_g(g, g_prms, coords, indices, wrt)

    assert_allclose(diff, expected, atol=1e-8, rtol=1e-8)


@pytest.mark.parametrize(
    "g, g_prms, coords, wrt, expected",
    [
        (
            _sch,
            (),
            [0., 2.5, np.pi / 6, np.pi / 2],
            2,
            0.
        ),
        (
            _kerr,
            (0.9,),
            [0., 25, np.pi / 2, 0.],
            1,
            0.0037790856763836173
        ),
        (
            _kerrnewman,
            (0.5, 0.1,),
            [0., 5.5, np.pi / 4, 0.],
            0,
            0.
        ),
    ],
)
def test_jacobian_g(g, g_prms, coords, wrt, expected):
    diff = _jacobian_g(g, g_prms, coords, wrt)[0, 0]

    assert_allclose(diff, expected, atol=1e-8, rtol=1e-8)


def test_ValueError1():
    try:
        d = _diff_g(
            g=_sch,
            g_prms=(),
            coords=[0., 2.5, np.pi / 6, np.pi / 2],
            indices=(2, 2, 2),
            wrt=2
        )

        assert False

    except ValueError:
        assert True


def test_ValueError2():
    try:
        d = _diff_g(
            g=_sch,
            g_prms=(),
            coords=[0., 2.5, np.pi / 6, np.pi / 2],
            indices=(2, 2),
            wrt=10
        )

        assert False

    except ValueError:
        assert True
