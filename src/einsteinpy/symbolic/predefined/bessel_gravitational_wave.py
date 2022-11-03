from sympy import Function, Symbol, besselj, cos, diag, exp, sin, sqrt, symbols

from einsteinpy.symbolic import constants
from einsteinpy.symbolic.metric import MetricTensor


def BesselGravitationalWave(C=symbols("C")):
    """
    Exact gravitational wave solution without diffraction.
    Class. Quantum Grav., 16:L75-78, 1999.
    D. Kramer.

    An exact solution describing an axisymmetric gravitational wave propagating in the
    z-direction in closed form. This solution to Einstein's vacuum field equations has the
    remarkable property that the curvature invariants decrease monotonically with increasing radial
    distance from the axis and vanish at infinity. The solution is regular at the symmetry axis.

    Parameters
    ----------
    C : ~sympy.core.basic.Basic or int or float
        Constant for Bessel metric, the choice of the constant is not really relavent for details see the paper. Defaults to 'C'.
    """
    coords = symbols("t rho phi z")
    t, rho, ph, z = coords

    # Useful helper functions, these wrap the Bessel functions that are used. C is some constant here
    U = C * besselj(rho, 0) * cos(t)
    K = (
        (1 / 2)
        * (C**2)
        * rho
        * (
            (rho * ((besselj(rho, 0) ** 2) + (besselj(rho, 1) ** 2)))
            - (2 * besselj(rho, 0) * besselj(rho, 1) * (cos(t) ** 2))
        )
    )

    # define the metric
    metric = diag(
        -1 * exp(-2 * U) * exp(2 * K),
        exp(-2 * U) * exp(2 * K),
        exp(-2 * U) * (rho**2),
        exp(2 * U),
    ).tolist()
    return MetricTensor(metric, coords, "ll", name="BesselGravitationalWaveMetric")
