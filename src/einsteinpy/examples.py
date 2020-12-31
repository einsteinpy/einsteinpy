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
    position = [40.0, np.pi / 2, 0.0]
    momentum = [0.0, 0.0, 3.83405]
    spin = 0.0

    # Calculating Geodesic
    geod = Timelike(
        position=position,
        momentum=momentum,
        a=spin,
        end_lambda=2000.0,
        step_size=0.5,
        return_cartesian=True,
        julia=True,
    )

    return geod
