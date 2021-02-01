import numpy as np

from einsteinpy.geodesic import Timelike

__all__ = ["precession"]


def precession():
    """
    An example to showcase the usage of the various modules in ``einsteinpy``.
    Here, we assume a Schwarzschild spacetime and obtain a test particle orbit, that
    shows apsidal precession.

    Returns
    -------
    geod: ~einsteinpy.geodesic.Timelike
        Timelike Geodesic, defining test particle trajectory

    """
    # Defining initial conditions
    metric = "Schwarzschild"
    position = [40.0, np.pi / 2, 0.0]
    momentum = [0.0, 0.0, 3.83405]

    # Calculating Geodesic
    geod = Timelike(
        metric=metric,
        metric_params=(),
        position=position,
        momentum=momentum,
        steps=5500,
        delta=1.0,
    )

    return geod
