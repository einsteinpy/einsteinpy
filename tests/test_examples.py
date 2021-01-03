import numpy as np
from unittest.mock import patch
from numpy.testing import assert_allclose

import einsteinpy


def Timelike_test(metric, metric_params, position, momentum, steps, delta):
    partial_lambdas = np.array([0, 1, 2, 3])
    partial_vecs = np.array(
        [
            [
                1.03068069e00,
                3.99997742e01,
                9.58510673e-02,
                2.44928680e-15,
                -9.79146613e-01,
                -3.87363444e-04,
                5.62571365e-19,
                3.83405000e00,
            ],
            [
                2.06136190e00,
                3.99987445e01,
                1.91699898e-01,
                2.44924485e-15,
                -9.79146730e-01,
                -9.48048366e-04,
                1.12515772e-18,
                3.83404999e00,
            ],
            [
                3.09204405e00,
                3.99971561e01,
                2.87545822e-01,
                2.44918275e-15,
                -9.79146663e-01,
                -1.21478762e-03,
                1.68776376e-18,
                3.83405001e00,
            ],
            [
                4.12272819e00,
                3.99949454e01,
                3.83388599e-01,
                2.44909661e-15,
                -9.79146558e-01,
                -1.80682765e-03,
                2.25040630e-18,
                3.83405000e00,
            ],
        ]
    )

    return partial_lambdas, partial_vecs


def test_precession_attr():
    with patch("einsteinpy.examples.Timelike", new=Timelike_test):
        p = einsteinpy.examples.precession()
        L = p[1][0, -1]

        assert p[0].shape[0] == p[1].shape[0]
        assert p[1].shape[1] == 8
        assert_allclose(p[1][:, 6], 0, atol=1e-12, rtol=1e-12)
        assert_allclose(p[1][:, 7], L, atol=1e-6, rtol=1e-6)
