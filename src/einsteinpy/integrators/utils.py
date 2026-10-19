"""
Utilities for Integration Module

"""
import numpy as np

from einsteinpy.utils.dual import _jacobian_g


def _PartHamFlow(g, g_prms, q, p, wrt):
    """
    Partial Hamiltonian Flow computed from the Metric

    Parameters
    ----------
    g : callable
        Metric Function
    g_prms : array_like
        Tuple of parameters to pass to the metric
        E.g., ``(a,)`` for Kerr
    q : array_like
        Initial 4-Position
    p : array_like
        Initial 3-Momentum
    wrt : int
        Coordinate, with respect to which, the derivative
        will be calculated
        Takes values from ``[0, 1, 2, 3]``

    Returns
    -------
    float
        Partial Hamiltonian Flow

    References
    ----------
    .. [1] Christian, Pierre and Chan, Chi-Kwan;
        "FANTASY: User-Friendly Symplectic Geodesic Integrator
        for Arbitrary Metrics with Automatic Differentiation";
        `2021 ApJ 909 67 <https://doi.org/10.3847/1538-4357/abdc28>`__

    """
    return _jacobian_g(g, g_prms, q, wrt) @ p @ p


def _flow_A(g, g_prms, q1, p1, q2, p2, delta=0.5):
    """
    Overall flow of Hamiltonian, :math:`H_A`

    Parameters
    ----------
    g : callable
        Metric Function
    g_prms : array_like
        Tuple of parameters to pass to the metric
        E.g., ``(a,)`` for Kerr
    q1 : array_like
        First copy of 4-Position
    p1 : array_like
        First copy of 3-Momentum
    q2 : array_like
        Second copy of 4-Position
    p2 : array_like
        Second copy of 3-Momentum
    delta : float
        Initial integration step-size
        Defaults to ``0.5``

    Returns
    -------
    float
        Hamiltonian Flow for :math:`H_A`

    References
    ----------
    .. [1] Christian, Pierre and Chan, Chi-Kwan;
        "FANTASY: User-Friendly Symplectic Geodesic Integrator
        for Arbitrary Metrics with Automatic Differentiation";
        `2021 ApJ 909 67 <https://doi.org/10.3847/1538-4357/abdc28>`__

    """
    dH1 = [0.5 * (_PartHamFlow(g, g_prms, q1, p2, i)) for i in range(4)]
    dp1 = np.array(dH1)
    p1_next = p1 - delta * dp1

    dH2 = g(q1, *g_prms) @ p2
    dq2 = np.array(dH2)
    q2_next = q2 + delta * dq2

    return q2_next, p1_next


def _flow_B(g, g_prms, q1, p1, q2, p2, delta=0.5):
    """
    Overall flow of Hamiltonian, :math:`H_B`

    Parameters
    ----------
    g : callable
        Metric Function
    g_prms : array_like
        Tuple of parameters to pass to the metric
        E.g., ``(a,)`` for Kerr
    q1 : array_like
        First copy of 4-Position
    p1 : array_like
        First copy of 3-Momentum
    q2 : array_like
        Second copy of 4-Position
    p2 : array_like
        Second copy of 3-Momentum
    delta : float
        Initial integration step-size
        Defaults to ``0.5``

    Returns
    -------
    float
        Hamiltonian Flow for :math:`H_B`

    References
    ----------
    .. [1] Christian, Pierre and Chan, Chi-Kwan;
        "FANTASY: User-Friendly Symplectic Geodesic Integrator
        for Arbitrary Metrics with Automatic Differentiation";
        `2021 ApJ 909 67 <https://doi.org/10.3847/1538-4357/abdc28>`__

    """
    dH2 = [
        0.5
        * (
            _PartHamFlow(
                g,
                g_prms,
                q2,
                p1,
                i,
            )
        )
        for i in range(4)
    ]
    dp2 = np.array(dH2)
    p2_next = p2 - delta * dp2

    dH1 = g(q2, *g_prms) @ p1
    dq1 = np.array(dH1)
    q1_next = q1 + delta * dq1

    return q1_next, p2_next


def _flow_mixed(q1, p1, q2, p2, delta=0.5, omega=1.0):
    """
    Mixed flow of Hamiltonian, :math:`\\tilde{H}`

    Parameters
    ----------
    q1 : array_like
        First copy of 4-Position
    p1 : array_like
        First copy of 3-Momentum
    q2 : array_like
        Second copy of 4-Position
    p2 : array_like
        Second copy of 3-Momentum
    delta : float
        Initial integration step-size
        Defaults to ``0.5``
    omega : float
        Coupling for Hamiltonian Flows
        Defaults to ``1.0``

    Returns
    -------
    float
        Hamiltonian Flow for :math:`\\tilde{H}`

    References
    ----------
    .. [1] Christian, Pierre and Chan, Chi-Kwan;
        "FANTASY: User-Friendly Symplectic Geodesic Integrator
        for Arbitrary Metrics with Automatic Differentiation";
        `2021 ApJ 909 67 <https://doi.org/10.3847/1538-4357/abdc28>`__

    """
    q_sum = q1 + q2
    q_dif = q1 - q2
    p_sum = p1 + p2
    p_dif = p1 - p2
    cos = np.cos(2.0 * omega * delta)
    sin = np.sin(2.0 * omega * delta)

    q1_next = 0.5 * (q_sum + (q_dif) * cos + (p_dif) * sin)
    p1_next = 0.5 * (p_sum + (p_dif) * cos - (q_dif) * sin)
    q2_next = 0.5 * (q_sum - (q_dif) * cos - (p_dif) * sin)
    p2_next = 0.5 * (p_sum - (p_dif) * cos + (q_dif) * sin)

    return q1_next, p1_next, q2_next, p2_next


def _Z(order):
    """
    Returns the constants for Yoshida Triple Jump.
    Used to compose higher order (even) integrators.

    References
    ----------
    .. [1] Yoshida, Haruo,
        "Construction of higher order symplectic integrators";
         Physics Letters A, vol. 150, no. 5-7, pp. 262-268, 1990.
        `DOI: <https://doi.org/10.1016/0375-9601(90)90092-3>`__

    """
    n = (order - 2) / 2
    x = 2 ** (1 / (2 * n + 1))
    Z0 = -x / (2 - x)
    Z1 = 1 / (2 - x)

    return Z0, Z1
