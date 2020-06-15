import numpy as np

from einsteinpy.coordinates.utils import four_position, stacked_vec
from einsteinpy.geodesic import Geodesic
from einsteinpy.metric import Schwarzschild


def perihelion():
    # ToDo: Use dedicated coordinate class, after refactoring `coordinate`
    M = 6e24  # Mass
    t = 10000  # Coordinate Time (has no effect in this case - Schwarzschild)
    x_vec = np.array([130.0, 0.0, 0.0])  # 3-Pos
    v_vec = np.array([0.0, 0.0, 1900.0])  # 3-Vel

    # Schwarzschild Metric Object
    ms_cov = Schwarzschild(M=M)
    # Getting Position 4-Vector
    x_4vec = four_position(t, x_vec)
    # Calculating Schwarzschild Metric at x_4vec
    ms_cov_mat = ms_cov.metric_covariant(x_4vec)
    # Getting stacked initial vector
    init_vec = stacked_vec(ms_cov_mat, t, x_vec, v_vec, time_like=True)

    # Calculating Geodesic
    geod = Geodesic(metric=ms_cov, init_vec=init_vec, end_lambda=0.002, step_size=5e-8)

    return geod
