import numpy as np

from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.geodesic import Geodesic
from einsteinpy.metric import Schwarzschild


def perihelion():
    """
    An example to showcase the usage of the various modules in ``einsteinpy``. \
    Here, we assume a Schwarzschild spacetime and obtain & plot the apsidal precession of \
    test particle orbit in it.

    Returns
    -------
    geod: ~einsteinpy.geodesic.Geodesic
        Geodesic defining test particle trajectory

    """
    # Mass of the black hole in SI
    M = 6e24

    # Defining the initial coordinates of the test particle
    # in SI
    sph = SphericalDifferential(
        t=10000.0,
        r=130.0,
        theta=np.pi / 2,
        phi=-np.pi / 8,
        v_r=0.0,
        v_th=0.0,
        v_p=1900.0,
    )

    # Schwarzschild Metric Object
    ms = Schwarzschild(M=M)

    # Getting the initial state of the test particle
    state = sph.state(metric=ms, time_like=True)

    # Calculating Geodesic
    geod = Geodesic(metric=ms, state=state, end_lambda=0.002, step_size=5e-8)

    return geod
