import numpy as np
import astropy.units as u
from numpy.testing import assert_allclose
import pytest

from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.coordinates import CartesianDifferential
from einsteinpy.coordinates import Spherical
from einsteinpy.metric import ReissnerNordstrom
from einsteinpy.utils import CoordinateError


def test_metric_reissnernordstrom():
    M = 2.0e30 * u.kg
    Q = 1.0e30 * u.C
    sph = Spherical(
        t=0. * u.s,
        r=1e6 * u.m,
        theta=4 * np.pi / 5 * u.rad,
        phi=0. * u.rad
    )

    # define expected metric values
    exp_metric = np.zeros(shape=(4, 4), dtype=float)

    exp_metric[0, 0] = 8.387172744E+31

    exp_metric[1, 1] = -1.071583007E-15

    exp_metric[2, 2] = -1E+12

    exp_metric[3, 3] = -3.454915028E+11

    ms = ReissnerNordstrom(coords=sph, M=M, Q=Q)

    metric = ms.metric_covariant(sph.position())

    for i in range(3):
        for j in range(3):
            assert(np.isclose(metric[i, j], metric[j, i],rtol=1e-8))

    assert_allclose(exp_metric, metric, rtol=1e-8)


def test_christoffels_reissnernordstrom():
    M = 2.0e30 * u.kg
    Q = 1.0e30 * u.C
    sph = Spherical(
        t=0. * u.s,
        r=1e6 * u.m,
        theta=4 * np.pi / 5 * u.rad,
        phi=0. * u.rad
    )

    # define expected christoffel values
    chl = np.zeros(shape=(4, 4, 4), dtype=float)

    # \Gamma(t,t,r) = \Gamma(t,r,t)
    chl[0, 0, 1] = -1.00000000E-06
    chl[0, 1, 0] = -1.00000000E-06

    # \Gamma(r,t,t)
    chl[1, 0, 0] = -7.82689973E+40

    # \Gamma(r,r,r)
    chl[1, 1, 1] = 1.00000000E-06

    # \Gamma(r,th,th)
    chl[1, 2, 2] = -9.33198822E+20

    # \Gamma(r,phi,phi)
    chl[1, 3, 3] = -3.22412264E+20

    # \Gamma(th,th,r)=\Gamma(th,r,th)
    chl[2, 2, 1] = 1.00000000E-06
    chl[2, 1, 2] = 1.00000000E-06

    # \Gamma(th,phi,phi)
    chl[2, 3, 3] = 4.75528258E-01

    # \Gamma(phi,phi,r)
    chl[3, 3, 1] = 1.00000000E-06
    chl[3, 1, 3] = 1.00000000E-06

    # \Gamma(phi,phi,th)
    chl[3, 3, 2] = -1.37638192E+00
    chl[3, 2, 3] = -1.37638192E+00

    ms = ReissnerNordstrom(coords=sph, M=M, Q=Q)

    christoffels = ms._christoffels(sph.position())

    for i in range(4):
        for j in range(4):
            for k in range(4):
                assert(np.isclose(christoffels[i, j, k],christoffels[i, k, j],rtol=1e-8))

    assert_allclose(chl, christoffels, rtol=1e-8)


def test_f_vec_s_reissnernordstrom():
    M = 2.0e30 * u.kg
    Q = 1.0e30 * u.C
    sph = SphericalDifferential(
        t=0. * u.s,
        r=1e6 * u.m,
        theta=4 * np.pi / 5 * u.rad,
        phi=0. * u.rad,
        v_r=1. * u.m / u.s,
        v_th=1. * u.rad / u.s,
        v_p=2e6 * u.rad / u.s
    )
    f_vec_expected = np.array(
        [
            0.000128363308154944, 1, 1, 2.00000000e+06,
            2.56726616310E-10, 2.57929819230E+33, -1.90211303259E+12, 5.50552368188E+06
        ]
    )

    ms = ReissnerNordstrom(coords=sph, M=M, Q=Q)
    state = np.hstack((sph.position(), sph.velocity(ms)))

    f_vec = ms._f_vec(0., state)

    assert isinstance(f_vec, np.ndarray)
    assert_allclose(f_vec_expected, f_vec, rtol=1e-8)


def test_raiseerror_reissnernordstrom():
    M = 2.0e30 * u.kg
    Q = 1.0e30 * u.C
    crds = CartesianDifferential(
        t=0. * u.s,
        x=1 * u.m,
        y=5 * u.m,
        z=0. * u.m,
        v_x=1. * u.m / u.s,
        v_y=1. * u.m / u.s,
        v_z=2 * u.m / u.s
    )

    with pytest.raises(CoordinateError):
        ms = ReissnerNordstrom(coords=crds, M=M, Q=Q)
        metric = ms.metric_covariant(crds.position())

    with pytest.raises(CoordinateError):
        ms = ReissnerNordstrom(coords=crds, M=M, Q=Q)
        chs = ms._christoffels(crds.position())

    with pytest.raises(CoordinateError):
        ms = ReissnerNordstrom(coords=crds, M=M, Q=Q)
        state = np.array([1, 1, 1, 1, 1, 1, 1, 1])
        f_vec = ms._f_vec(0., state)
