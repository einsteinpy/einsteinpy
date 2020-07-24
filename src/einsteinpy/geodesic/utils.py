import numpy as np

from einsteinpy.ijit import jit


# @jit(parallel=True)
def _calc_state(alpha, beta, r, theta, a, E):
    """
    Jit-ed funtion for calculating state

    Source: RAPTOR (?????)

    """
    # Useful expressions
    sint, cost, tant = np.sin(theta), np.cos(theta), np.tan(theta)
    sint2, cost2 = sint ** 2, cost ** 2
    r2, a2 = r ** 2, a ** 2
    sg = r2 + a2 * cost2
    dl = r2 + a2 - 2.0 * r

    # Constants of Motion
    L = -alpha * E * sint
    Q = (E ** 2) * (beta ** 2 + cost2 * (alpha ** 2 - 1.0))

    # Co- & Contra-variant wave vectors
    k_d = np.zeros(4)
    k_u = np.zeros(4)

    k_d[0] = -E
    k_d[2] = np.sign(beta) * np.sqrt(np.abs(Q - (L / tant) ** 2 + cost2 * (E ** 2)))
    k_d[3] = L

    # Relevant metric components, required to build k_u
    g_dd_11 = sg / dl
    g_uu_00 = -(((r2 + a2) ** 2) - (dl * a2 * sint2)) / (sg * dl)
    g_uu_03 = -2.0 * a * r / (sg * dl)
    g_uu_33 = (dl - a2 * sint2) / (sg * dl * sint2)
    g_uu_22 = 1.0 / sg

    # Building k_u
    k_u[0] = g_uu_00 * k_d[0] + g_uu_03 * k_d[3]
    k_u[3] = g_uu_33 * k_d[3] + g_uu_03 * k_d[0]
    k_u[2] = g_uu_22 * k_d[2]
    k_u[1] = np.sqrt((-k_u[0] * k_d[0] - k_u[2] * k_d[2] - k_u[3] * k_d[3]) / g_dd_11)

    # Constructing State Vector
    x_u = np.array([0.0, r, theta, np.pi / 2])  # Position Vector

    state = np.hstack((x_u, k_u))

    return state


# @jit(parallel=True)
def _christoffels_M(a, r, th):
    """
    Jit-ed funtion for calculating Christoffel Symbols \
    in Kerr Metric in M-Units & Boyer-Lindquist Coordinates

    Source: RAPTOR (?????)

    """
    # Useful expressions
    sint, cost = np.sin(th), np.cos(th)
    sint2, cost2 = sint ** 2, cost ** 2
    r2, a2 = r ** 2, a ** 2
    sg = r2 + a2 * cost2
    dl = r2 + a2 - 2 * r
    sg2, sg3 = sg ** 2, sg ** 3

    chl = np.zeros((4, 4, 4))

    # Non-zero components
    chl[0, 0, 1] = (r2 + a2) / (sg2 * dl) * (2.0 * r2 - sg)
    chl[0, 0, 2] = -2.0 * a2 * r * sint * cost / sg2
    chl[0, 1, 3] = -a * sint2 / (sg * dl) * (2.0 * r2 / sg * (r2 + a2) + r2 - a2)
    chl[0, 2, 3] = 2.0 * a2 * a * r * sint2 * sint * cost / sg2
    chl[1, 0, 0] = dl / sg3 * (2.0 * r2 - sg)
    chl[1, 0, 3] = -a * dl * sint2 / sg3 * (2.0 * r2 - sg)
    chl[1, 1, 1] = (1.0 - r) / dl + r / sg
    chl[1, 1, 2] = -a2 * sint * cost / sg
    chl[1, 2, 2] = -r * dl / sg
    chl[1, 3, 3] = -dl * sint2 / sg * (r - a2 * sint2 / sg2 * (2.0 * r2 - sg))
    chl[2, 0, 0] = -2.0 * a2 * r * sint * cost / sg3
    chl[2, 0, 3] = 2.0 * a * r * (r2 + a2) * sint * cost / sg3
    chl[2, 1, 1] = a2 * sint * cost / (sg * dl)
    chl[2, 1, 2] = r / sg
    chl[2, 3, 3] = (
        -sint
        * cost
        / sg3
        * ((r2 + a2) * (((r2 + a2) ** 2) - dl * a2 * sint2) - sg * dl * a2 * sint2)
    )
    chl[3, 0, 1] = a / (sg2 * dl) * (2.0 * r2 - sg)
    chl[3, 0, 2] = -2.0 * a * r * cost / (sg2 * sint)
    chl[3, 1, 3] = r / sg - a2 * sint2 / (sg * dl) * (r - 1.0 + 2.0 * r2 / sg)
    chl[3, 2, 3] = cost / sint * (1.0 + 2.0 * a2 * r * sint2 / sg2)
    chl[2, 2, 2] = chl[1, 1, 2]

    return chl


# @jit(parallel=True)
def _f_geod_M(chl, vals, vec):
    """
    Jit-ed funtion for calculating ``f_vec`` in \
    Kerr Metric in M-Units & Boyer-Lindquist Coordinates

    Source: RAPTOR (?????)

    """

    vals[:4] = vec[4:]

    vals[4] = -2.0 * (
        chl[0, 0, 1] * vec[4] * vec[5]
        + chl[0, 0, 2] * vec[4] * vec[6]
        + chl[0, 1, 3] * vec[5] * vec[7]
        + chl[0, 2, 3] * vec[6] * vec[7]
    )
    vals[5] = -1.0 * (
        chl[1, 0, 0] * vec[4] * vec[4]
        + 2 * chl[1, 0, 3] * vec[4] * vec[7]
        + chl[1, 1, 1] * vec[5] * vec[5]
        + 2 * chl[1, 1, 2] * vec[5] * vec[6]
        + chl[1, 2, 2] * vec[6] * vec[6]
        + chl[1, 3, 3] * vec[7] * vec[7]
    )
    vals[6] = -1.0 * (
        chl[2, 0, 0] * vec[4] * vec[4]
        + 2 * chl[2, 0, 3] * vec[4] * vec[7]
        + chl[2, 1, 1] * vec[5] * vec[5]
        + 2 * chl[2, 1, 2] * vec[5] * vec[6]
        + chl[2, 2, 2] * vec[6] * vec[6]
        + chl[2, 3, 3] * vec[7] * vec[7]
    )
    vals[7] = -2.0 * (
        chl[3, 0, 1] * vec[4] * vec[5]
        + chl[3, 0, 2] * vec[4] * vec[6]
        + chl[3, 1, 3] * vec[5] * vec[7]
        + chl[3, 2, 3] * vec[6] * vec[7]
    )

    return vals
