import numpy as np
from astropy import units as u

from einsteinpy.coordinates import SphericalDifferential
from einsteinpy.geodesic import Timelike
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
    M = 6e24 * u.kg

    # Defining the initial coordinates of the test particle
    # in SI
    sph = SphericalDifferential(
        t=10000.0 * u.s,
        r=130.0 * u.m,
        theta=np.pi / 2 * u.rad,
        phi=-np.pi / 8 * u.rad,
        v_r=0.0 * u.m / u.s,
        v_th=0.0 * u.rad / u.s,
        v_p=1900.0 * u.rad / u.s,
    )

    # Schwarzschild Metric Object
    ms = Schwarzschild(coords=sph, M=M)

    # Calculating Geodesic
    geod = Timelike(metric=ms, coords=sph, end_lambda=0.002, step_size=5e-8)

    return geod
