import numpy as np
import pytest
from numpy.testing import assert_allclose

from einsteinpy.geodesic.utils import _P, _sch, _kerr, _kerrnewman


@pytest.mark.parametrize(
    "g, g_prms, q, p, time_like, expected",
    [
        (
            _sch,
            (),
            [0., 2.5, np.pi / 6, np.pi / 2],
            [0.1, 2., 2.],
            True,
            [-0.9167333309092672, 0.1, 2., 2.]
        ),
        (
            _kerr,
            (0.9,),
            [0., 25, np.pi / 2, 0.],
            [0., 0, 2.427],
            False,
            [-0.09333038545829885, 0., 0, 2.427]
        ),
        (
            _kerrnewman,
            (0.5, 0.1,),
            [0., 5.5, np.pi / 4, 0.],
            [0.1, -0.2, -4.],
            True,
            [-1.1220908695019767, 0.1, -0.2, -4.]
        ),
    ],
)
def test_P(g, g_prms, q, p, time_like, expected):
    P = _P(g, g_prms, q, p, time_like)

    assert_allclose(P, expected, atol=1e-8, rtol=1e-8)


def _generic(x_vec, *params):
    """
    A contravariant metric with non-zero off-diagonal terms in every slot,
    used to check that ``_P`` does not assume Kerr-like sparsity

    """
    r, th = x_vec[1], x_vec[2]

    g = np.zeros(shape=(4, 4))

    tmp = 1.0 - (2 / r)
    g[0, 0] = -1 / tmp
    g[1, 1] = tmp
    g[2, 2] = 1 / (r**2)
    g[3, 3] = 1 / ((r * np.sin(th)) ** 2)

    g[0, 1] = g[1, 0] = params[0]
    g[0, 2] = g[2, 0] = params[0] / 2
    g[1, 2] = g[2, 1] = params[0] / 4
    g[1, 3] = g[3, 1] = params[0] / 8
    g[2, 3] = g[3, 2] = params[0] / 16

    return g


@pytest.mark.parametrize(
    "g, g_prms, q, p, time_like",
    [
        (_sch, (), [0., 2.5, np.pi / 6, np.pi / 2], [0.1, 2., 2.], True),
        (_kerr, (0.9,), [0., 25, np.pi / 2, 0.], [0., 0, 2.427], False),
        (_kerrnewman, (0.5, 0.1,), [0., 5.5, np.pi / 4, 0.], [0.1, -0.2, -4.], True),
        (_generic, (0.1,), [0., 6., np.pi / 2, 0.], [0.12, 0.04, 3.5], True),
        (_generic, (0.1,), [0., 6., np.pi / 2, 0.], [0.12, 0.04, 3.5], False),
        (_generic, (0.3,), [0., 8., np.pi / 3, 0.], [-0.2, 0.5, 1.5], True),
    ],
)
def test_P_satisfies_normalization(g, g_prms, q, p, time_like):
    """
    ``_P`` must return a 4-Momentum satisfying
    :math:`g^{\\mu\\nu} p_\\mu p_\\nu = -m^2`, with :math:`m = 1` for
    time-like and :math:`m = 0` for null geodesics

    """
    P = _P(g, g_prms, q, p, time_like)
    guu = np.array(g(q, *g_prms), dtype=float)

    assert_allclose(P[1:], p, atol=1e-8, rtol=1e-8)
    assert_allclose(P @ guu @ P, -int(time_like), atol=1e-10)


@pytest.mark.parametrize(
    "x",
    [
        (
            [0., 2.5, np.pi / 6, np.pi / 2],
        ),
        (
            [0., 25, np.pi / 2, 0.],
        ),
        (
            [0., 5.5, np.pi / 4, 0.],
        ),
    ],
)
def test_metrics(x):
    x = x[0]

    # a = Q = 0.
    s = _sch(x).astype(float)
    k = _kerr(x, 0.).astype(float)
    kn = _kerrnewman(x, 0., 0.).astype(float)
    assert_allclose(s, k, atol=1e-8, rtol=1e-8)
    assert_allclose(k, kn, atol=1e-8, rtol=1e-8)
    assert_allclose(kn, s, atol=1e-8, rtol=1e-8)

    # Non-zero Spin
    k = _kerr(x, 0.4).astype(float)
    kn = _kerrnewman(x, 0.4, 0.).astype(float)
    assert_allclose(k, kn, atol=1e-8, rtol=1e-8)
