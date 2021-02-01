import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.utils.dual import DualNumber, _deriv, _diff_g, _jacobian_g
from einsteinpy.geodesic.utils import _sch, _kerr, _kerrnewman


def test_str_repr_DualNumber():
    x = DualNumber(1.0, 2.0)

    assert str(x) == repr(x) == "DualNumber(1.0, 2.0)"


def test_DualNumber1():
    x = DualNumber(1.0, 2.0)
    y = DualNumber(2.0, 3.0)
    c = 2.0

    add_ = x + y
    assert_allclose([add_.val, add_.deriv], [3.0, 5.0], atol=1e-8, rtol=1e-8)

    radd_ = x.__radd__(c)
    assert_allclose([radd_.val, radd_.deriv], [3.0, 2.0], atol=1e-8, rtol=1e-8)

    sub_ = x - y
    assert_allclose([sub_.val, sub_.deriv], [-1.0, -1.0], atol=1e-8, rtol=1e-8)

    rsub_ = x.__rsub__(c)
    assert_allclose([rsub_.val, rsub_.deriv], [1.0, -2.0], atol=1e-8, rtol=1e-8)

    mul_ = x * y
    assert_allclose([mul_.val, mul_.deriv], [2.0, 7.0], atol=1e-8, rtol=1e-8)

    rmul_ = x.__rmul__(c)
    assert_allclose([rmul_.val, rmul_.deriv], [2.0, 4.0], atol=1e-8, rtol=1e-8)


def test_DualNumber2():
    x = DualNumber(1.0, 2.0)
    y = DualNumber(2.0, 3.0)
    c = 2.0

    add_ = x + c
    assert_allclose([add_.val, add_.deriv], [3.0, 2.0], atol=1e-8, rtol=1e-8)

    radd_ = x.__radd__(y)
    assert_allclose([radd_.val, radd_.deriv], [3.0, 5.0], atol=1e-8, rtol=1e-8)

    sub_ = x - c
    assert_allclose([sub_.val, sub_.deriv], [-1.0, 2.0], atol=1e-8, rtol=1e-8)

    rsub_ = x.__rsub__(y)
    assert_allclose([rsub_.val, rsub_.deriv], [1.0, 1.0], atol=1e-8, rtol=1e-8)

    mul_ = x * c
    assert_allclose([mul_.val, mul_.deriv], [2.0, 4.0], atol=1e-8, rtol=1e-8)

    rmul_ = x.__rmul__(y)
    assert_allclose([rmul_.val, rmul_.deriv], [2.0, 7.0], atol=1e-8, rtol=1e-8)


def test_DualNumber3():
    x = DualNumber(1.0, 2.0)
    y = DualNumber(2.0, 3.0)
    c = 2.0

    truediv_ = x / y
    assert_allclose([truediv_.val, truediv_.deriv], [0.5, 0.25], atol=1e-8, rtol=1e-8)

    rtruediv1_ = x / c
    assert_allclose(
        [rtruediv1_.val, rtruediv1_.deriv], [0.5, 1.0], atol=1e-8, rtol=1e-8
    )

    rtruediv2_ = c / x
    assert_allclose(
        [rtruediv2_.val, rtruediv2_.deriv], [2.0, -4.0], atol=1e-8, rtol=1e-8
    )

    assert x == x

    assert x != y

    pow_ = y ** 2
    assert_allclose([pow_.val, pow_.deriv], [4.0, 12.0], atol=1e-8, rtol=1e-8)

    neg_ = -y
    assert_allclose([neg_.val, neg_.deriv], [-2.0, -3.0], atol=1e-8, rtol=1e-8)


def test_DualNumber4():
    x = DualNumber(1.0, 2.0)
    y = DualNumber(2.0, 3.0)
    c = 2.0

    truediv_ = x / c
    assert_allclose([truediv_.val, truediv_.deriv], [0.5, 1.0], atol=1e-8, rtol=1e-8)

    rtruediv1_ = x.__rtruediv__(y)
    assert_allclose(
        [rtruediv1_.val, rtruediv1_.deriv], [2.0, -1.0], atol=1e-8, rtol=1e-8
    )

    x = DualNumber(0.0, 2.0)
    y = DualNumber(0.0, 4.0)

    truediv_ = x / y
    assert_allclose([truediv_.val, truediv_.deriv], [0.5, 0.0], atol=1e-8, rtol=1e-8)

    rtruediv1_ = x.__rtruediv__(y)
    assert_allclose(
        [rtruediv1_.val, rtruediv1_.deriv], [2.0, 0.0], atol=1e-8, rtol=1e-8
    )


def test_DualNumber5():
    x = DualNumber(1.0, 2.0)
    y = DualNumber(2.0, 3.0)
    c = 2.0

    sin = x.sin()
    assert_allclose(
        [sin.val, sin.deriv],
        [0.8414709848078965, 1.0806046117362795],
        atol=1e-8,
        rtol=1e-8,
    )

    cos = y.cos()
    assert_allclose(
        [cos.val, cos.deriv],
        [-0.4161468365471424, -2.727892280477045],
        atol=1e-8,
        rtol=1e-8,
    )

    tan = x.tan()
    assert_allclose(
        [tan.val, tan.deriv],
        [1.557407724654902, 6.851037641629518],
        atol=1e-8,
        rtol=1e-8,
    )

    log = y.log()
    assert_allclose(
        [log.val, log.deriv], [0.6931471805599453, 1.5], atol=1e-8, rtol=1e-8
    )

    exp = x.exp()
    assert_allclose(
        [exp.val, exp.deriv],
        [2.718281828459045, 5.43656365691809],
        atol=1e-8,
        rtol=1e-8,
    )


def test_deriv1():
    df = _deriv(lambda y: np.cos(y) * np.sin(y) ** 2, 1.0)

    assert_allclose(df, -0.10452774015707361, atol=1e-8, rtol=1e-8)


def test_deriv2():
    df = _deriv(lambda y: 0.0, 0.0)

    assert_allclose(df, 0.0, atol=1e-8, rtol=1e-8)


@pytest.mark.parametrize(
    "g, g_prms, coords, indices, wrt, expected",
    [
        (_sch, (), [0.0, 2.5, np.pi / 6, np.pi / 2], (2, 2), 2, 0.0),
        (_sch, (), [0.0, 2.5, np.pi / 6, np.pi / 2], (2, 2), 1, -0.128),
        (_kerr, (0.9,), [0.0, 25, np.pi / 2, 0.0], (0, 3), 1, 1.542519202468211e-05),
        (_kerr, (0.9,), [0.0, 25, np.pi / 2, 0.0], (0, 3), 3, 0.0),
        (
            _kerrnewman,
            (
                0.5,
                0.1,
            ),
            [0.0, 5.5, np.pi / 4, 0.0],
            (3, 0),
            0,
            0.0,
        ),
        (
            _kerrnewman,
            (
                0.5,
                0.1,
            ),
            [0.0, 5.5, np.pi / 4, 0.0],
            (3, 0),
            2,
            -7.631639775075298e-05,
        ),
    ],
)
def test_diff_g(g, g_prms, coords, indices, wrt, expected):
    diff = _diff_g(g, g_prms, coords, indices, wrt)

    assert_allclose(diff, expected, atol=1e-8, rtol=1e-8)


@pytest.mark.parametrize(
    "g, g_prms, coords, wrt, expected",
    [
        (_sch, (), [0.0, 2.5, np.pi / 6, np.pi / 2], 2, 0.0),
        (_kerr, (0.9,), [0.0, 25, np.pi / 2, 0.0], 1, 0.0037790856763836173),
        (
            _kerrnewman,
            (
                0.5,
                0.1,
            ),
            [0.0, 5.5, np.pi / 4, 0.0],
            0,
            0.0,
        ),
    ],
)
def test_jacobian_g(g, g_prms, coords, wrt, expected):
    diff = _jacobian_g(g, g_prms, coords, wrt)[0, 0]

    assert_allclose(diff, expected, atol=1e-8, rtol=1e-8)


def test_ValueError1():
    with pytest.raises(ValueError):
        d = _diff_g(
            g=_sch,
            g_prms=(),
            coords=[0.0, 2.5, np.pi / 6, np.pi / 2],
            indices=(2, 2, 2),
            wrt=2,
        )


def test_ValueError2():
    with pytest.raises(ValueError):
        d = _diff_g(
            g=_sch,
            g_prms=(),
            coords=[0.0, 2.5, np.pi / 6, np.pi / 2],
            indices=(2, 2),
            wrt=10,
        )
