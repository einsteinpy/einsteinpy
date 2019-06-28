from unittest import mock

from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy.hypersurface import SchwarzschildEmbedding


def test_gradient():
    M = 4e24 * u.kg
    r = 1200
    obj = SchwarzschildEmbedding(M)
    grad = obj.gradient(r)
    grad_test_value = 13.7267
    assert_allclose(grad_test_value, grad, 1e-4)


def test_radial_coord():
    M = 5.972e24 * u.kg
    r = 170000
    obj = SchwarzschildEmbedding(M)
    R = obj.radial_coord(r)
    R_test_value = 170059.7514
    assert_allclose(R_test_value, R, 1e-4)


def test_get_values():
    M = 5.972e24 * u.g
    obj = SchwarzschildEmbedding(M)
    x_axis, y_axis = obj.get_values(100)
    x_axis = x_axis[:10]
    y_axis = y_axis[:10]
    x_axis_test_value = [
        0.31031436748282454,
        0.3103211992691816,
        0.31033806032413,
        0.3103646908134154,
        0.3104008399665836,
        0.31044626567977485,
        0.3105007341394409,
        0.3105640194657007,
        0.31063590337413993,
        0.3107161748549419,
    ]
    y_axis_test_value = [
        0,
        0.000996913602540167,
        0.001987250685704889,
        0.002971124505173177,
        0.003948645664953687,
        0.004919922187844408,
        0.005885059584344937,
        0.006844160919989087,
        0.00779732688107651,
        0.008744655838791766,
    ]
    assert_allclose(x_axis_test_value, x_axis, 1e-4)
    assert_allclose(y_axis_test_value, y_axis, 1e-4)


def test_get_values_surface():
    M = 5.972e24 * u.g
    obj = SchwarzschildEmbedding(M)
    X, Y, Z = obj.get_values_surface(100)
    X = X[1][:10].tolist()
    Y = Y[1][:10].tolist()
    Z = Z[1].tolist()
    X_test_value = [
        0.3085563745933694,
        0.3085631676762918,
        0.3085799332099953,
        0.3086064128322384,
        0.30864235719323163,
        0.3086875255606821,
        0.30874168544564184,
        0.30880461224788464,
        0.3088760889196209,
        0.3089559056464467,
    ]
    Y_test_value = [
        0.03298439576620688,
        0.03298512194059198,
        0.03298691416089414,
        0.03298974480842477,
        0.03299358722791822,
        0.03299841568531133,
        0.03300420532774708,
        0.03301093214566557,
        0.03301857293685553,
        0.03302710527234805,
    ]
    Z_test_value = 0.000996913602540167
    assert_allclose(X_test_value, X, 1e-4)
    assert_allclose(Y_test_value, Y, 1e-4)
    assert_allclose(Z_test_value, Z, 1e-4)
