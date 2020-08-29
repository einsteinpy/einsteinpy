import numpy as np
import pytest

from unittest.mock import patch
from numpy.testing import assert_allclose

import einsteinpy
import einsteinpy_geodesics


def solveSystem_test(q, p, params, end_lambda, step_size):
    partial_lambdas = np.array([0.0000e+00, 5.0000e-01, 1.0000e+00, 1.9990e+03, 1.9995e+03, 2.0000e+03])
    partial_vecs = np.array(
        [[ 4.00000000e+01,  0.00000000e+00,  2.44929360e-15, 0.00000000e+00,  0.00000000e+00,  3.83405000e+00],
        [ 3.99999197e+01, -4.79255517e-02,  2.44929044e-15, 2.17126283e-04, -2.81284902e-19,  3.83405000e+00],
        [ 3.99996789e+01, -9.58508493e-02,  2.44928097e-15, 4.34253495e-04, -5.62571255e-19,  3.83405000e+00],
        [-6.00084085e+00, -2.52448871e+00,  3.98636802e-16, -2.33673001e-01, -9.49738949e-15,  3.83405000e+00],
        [-6.18465627e+00, -2.28046143e+00,  4.03625070e-16, -2.36036036e-01, -9.50787744e-15,  3.83405000e+00],
        [-6.35745272e+00, -2.03236148e+00,  4.08689510e-16, -2.38256302e-01, -9.51810736e-15,  3.83405000e+00]]
    )

    return partial_lambdas, partial_vecs

def test_mock_julia_backend():
    q0 = [2.5, np.pi / 2, 0.]
    p0 = [0., 0., -8.5]
    a = 0.9
    end_lambda = 10.
    step_size = 0.005

    with patch('einsteinpy.geodesic.geodesic.solveSystem', new=solveSystem_test):
        geod = einsteinpy.geodesic.Geodesic(
            position=q0,
            momentum=p0,
            a=a,
            end_lambda=end_lambda,
            step_size=step_size,
            time_like=False,
            return_cartesian=True,
            julia=True
        )

        traj = geod.trajectory
        L = traj[1][0, 5]

        assert traj[0].shape[0] == traj[1].shape[1]
        assert traj[1].shape[1] == 6
        assert_allclose(traj[1][:, 4], 0, atol=1e-12, rtol=1e-12)
        assert_allclose(traj[1][:, 5], L, atol=1e-12, rtol=1e-12)
