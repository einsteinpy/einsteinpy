import numpy as np

"""
Utilities for Geodesic Module

Unit System: M-Units => G = c = 1
Coordinate System: Boyer-Lindquist Coordinates

"""

# Not used
def _sphToCart(r, th, p):
    """
    Converts Spherical Polar Coordinates into Cartesian Coordinates

    """
    xs = r * np.sin(th) * np.cos(p)
    ys = r * np.sin(th) * np.sin(p)
    zs = r * np.cos(th)

    return xs, ys, zs


# Not used
def _norm_U(y, a):
    """
    Returns norm of 4-Velocity

    """
    u = y[4:]
    r, th = y[1], y[2]

    g_dd = _g_dd(r, th, a)

    return u @ g_dd @ u


def _v0(g_dd, v1, v2, v3):
    """
    Returns Timelike component of 4-Velocity for Null Geodesics

    """
    g = g_dd

    # Defining coefficients for quadratic equation
    A = g[0, 0]
    B = 2 * (g[0, 1] * v1 + g[0, 2] * v2 + g[0, 3] * v3)
    C = (
        (g[1, 1] * np.square(v1) + g[2, 2] * np.square(v2) + g[3, 3] * np.square(v3))
        + 2 * v1 * (g[1, 2] * v2 + g[1, 3] * v3)
        + 2 * v2 * g[2, 3] * v3
    )
    D = np.square(B) - (4 * A * C)

    v_t = (-B + np.sqrt(D)) / (2 * A)

    return v_t


def _g_dd(r, th, a):
    """
    Returns covariant Kerr metric in Boyer-Lindquist Coordinates \
    and M-Units (G = c = M = 1)

    """
    r2, a2 = np.square(r), np.square(a)
    sint2, cost2 = np.square(np.sin(th)), np.square(np.cos(th))

    sg = r2 + a2 * cost2
    dl = r2 - 2 * r + a2

    g_dd = np.zeros((4, 4), dtype=np.float64)

    g_dd[0, 0] = -1 + 2 * r / sg
    g_dd[1, 1] = sg / dl
    g_dd[2, 2] = sg
    g_dd[3, 3] = (dl + (2 * r * (a2 + r2)) / sg) * sint2
    g_dd[0, 3] = g_dd[3, 0] = (-2 * a * r * sint2) / sg

    return g_dd


# Not used
def _g_uu(r, th, a):
    """
    Returns contravariant Kerr metric in \
    Boyer-Lindquist Coordinates \
    and M-Units (G = c = M = 1)

    """
    g_dd = g_dd(r, th, a)
    g_uu = np.linalg.inv(g_dd)

    return g_uu


# For use with _f_vec
def _gamma_udd(r, th, a):
    """
    Returns Christoffel Symbols in Boyer-Lindquist Coordinates \
    and M-Units (G = c = M = 1)

    """
    sg, dl = sigma(r, th, a), delta(r, a)
    sg2 = np.square(sg)
    r2, a2 = np.square(r), np.square(a)
    r4, a4 = np.power(r, 4), np.power(a, 4)
    sint, cost = np.sin(th), np.cos(th)
    sint2, cost2 = np.square(sint), np.square(cost)
    sint4, cost4 = np.power(sint, 4), np.power(cost, 4)
    sin2t, cos2t = np.sin(2 * th), np.cos(2 * th)

    denom = (
        2 * a2 * (2 + dl) * r2 * cost2
        + a4 * dl * cost4
        + r2 * ((-2 + r) * r2 * r + a2 * (4 + r2) - 8 * a2 * cos2t)
    )

    gamma_udd = np.zeros((4, 4, 4))

    gamma_udd[0, 1, 0] = (
        (r2 - a2 * cost2)
        * (a4 + 3 * a2 * (-2 + r) * r + 2 * r4 + a2 * (a2 + r * (6 + r)) * cos2t)
        / (2 * sg * denom)
    )

    gamma_udd[0, 2, 0] = -(
        a2
        * r
        * cost
        * (a4 * 2 * (-8 + r) * r2 * r + a2 * r * (-14 + 3 * r) + a2 * dl * cos2t)
        * sint
    ) / (sg * denom)

    gamma_udd[0, 3, 1] = (
        2 * a * (-r2 * (a2 + 3 * r2) + a2 * (a2 - r2) * cost2) * sint2
    ) / denom

    gamma_udd[0, 3, 2] = (4 * a2 * a * dl * r * cost * sint2 * sint) / denom

    gamma_udd[1, 0, 0] = -(dl * (-r2 + a2 * cost2)) / (sg2 * sg)

    gamma_udd[1, 1, 1] = ((a2 - r) * r - a2 * (-1 + r) * cost2) / (dl * sg)

    gamma_udd[1, 2, 1] = -(a2 * cost * sint) / sg

    gamma_udd[1, 2, 2] = -(dl * r) / sg

    gamma_udd[1, 3, 0] = (2 * a * dl * (-r2 + a2 * cost2) * sint2) / (sg2 * sg)

    gamma_udd[1, 3, 3] = -(
        dl
        * (
            -a2 * r2
            + r4 * r
            + a2 * (a2 + r2 + 2 * r2 * r) * cost2
            + a4 * (-1 + r) * cost4
        )
        * sint2
    ) / (sg2 * sg)

    gamma_udd[2, 0, 0] = -(2 * a2 * r * cost * sint) / (sg2 * sg)

    gamma_udd[2, 1, 1] = (a2 * cost * sint) / (dl * sg)

    gamma_udd[2, 2, 1] = r / sg

    gamma_udd[2, 2, 2] = -(a2 * cost * sint) / sg

    gamma_udd[2, 3, 0] = (2 * a * r * (a2 + r2) * sin2t) / (sg2 * sg)

    gamma_udd[2, 3, 3] = (
        cost * sint * (-sg2 * dl - sg * (2 * r * (a2 + r2)))
        - 2 * a2 * r * (a2 + r2) * sint2
    ) / (sg2 * sg)

    gamma_udd[3, 1, 0] = (2 * a * r2 - 2 * a2 * a * cost2) / denom

    gamma_udd[3, 2, 0] = -(4 * a * dl * r * (cost / sint)) / denom

    gamma_udd[3, 3, 1] = (
        (
            2 * a2 * r2 * r
            - a2 * r4
            - 2 * r4 * r2
            + r4 * r2 * r
            - a2 * r * (2 * a2 + r2 * (2 + 3 * r - 3 * r2)) * cost2
        )
        + (a4 * (a2 + r * (2 - 2 * r + 3 * r2)) * cost4)
        + (a4 * a2 * (-1 + r) * cost4 * cost2)
        - (8 * a2 * r2 * r * sint2)
        + (2 * a4 * r * np.square(sin2t))
    ) / (sg * denom)

    gamma_udd[3, 3, 2] = (
        (cost / sint)
        * (
            (a4 * r2 * (4 + 3 * a2 - 6 * r + 3 * r2) * cost4)
            + (a4 * a2 * dl * cost4 * cost2)
            + (
                r2
                * (
                    (-2 + r) * r4 * r
                    + a4 * (6 + r)
                    + a2 * r2 * (2 + r + r2)
                    - a2 * (6 + r) * (a2 + r2) * cos2t
                )
            )
            + (
                a2
                * r
                * cost2
                * (
                    a2 * r * (-4 + 3 * r2)
                    + r2 * r * (4 - 6 * r + 3 * r2)
                    + 2 * a2 * (a2 + r2) * sint2
                )
            )
        )
        / (sg * denom)
    )

    gamma_udd[0, 0, 1] = gamma_udd[0, 1, 0]
    gamma_udd[0, 0, 2] = gamma_udd[0, 2, 0]
    gamma_udd[0, 1, 3] = gamma_udd[0, 3, 1]
    gamma_udd[0, 2, 3] = gamma_udd[0, 3, 2]
    gamma_udd[1, 1, 2] = gamma_udd[1, 2, 1]
    gamma_udd[1, 0, 3] = gamma_udd[1, 3, 0]
    gamma_udd[2, 1, 2] = gamma_udd[2, 2, 1]
    gamma_udd[2, 0, 3] = gamma_udd[2, 3, 0]
    gamma_udd[3, 0, 1] = gamma_udd[3, 1, 0]
    gamma_udd[3, 0, 2] = gamma_udd[3, 2, 0]
    gamma_udd[3, 1, 3] = gamma_udd[3, 3, 1]
    gamma_udd[3, 2, 3] = gamma_udd[3, 3, 2]

    return gamma_udd


# For use with _gamma_udd
def _f_vec(vec, a):
    """
    Returns RHS of the Geodesic Equations
    
    """
    vals = np.zeros(shape=vec.shape, dtype=vec.dtype)
    # chl = _gamma_udd(vec[1], vec[2])

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


# Standalone
def _f_vec_geod(t, vec, a):
    """
    Returns RHS of Geodesic Equations
    
    """
    vals = np.zeros(shape=vec.shape, dtype=vec.dtype)
    r, th = vec[1], vec[2]
    tdot = vec[4]
    rdot = vec[5]
    thdot = vec[6]
    pdot = vec[7]

    r2, a2 = np.square(r), np.square(a)
    r4, a4 = np.power(r, 4), np.power(a, 4)
    sint, cost = np.sin(th), np.cos(th)
    sint2, cost2 = np.square(sint), np.square(cost)
    cost4 = np.power(cost, 4)
    sin2t, cos2t, cos4t = np.sin(2 * th), np.cos(2 * th), np.cos(4 * th)

    vals[:4] = vec[4:]

    vals[4] = (
        (
            -2
            * a2
            * r
            * thdot
            * (
                -2
                * (
                    a4
                    + 2 * (-8 + r) * r2 * r
                    + a2 * r * (-14 + 3 * r)
                    + a2 * (a2 + (-2 + r) * r) * cos2t
                )
                * sin2t
                * tdot
                + 8
                * a
                * (a2 + (-2 + r) * r)
                * cost
                * (a2 + 2 * r2 + a2 * cos2t)
                * sin2t
                * sint
                * pdot
            )
        )
        + rdot
        * (
            (
                3 * a4 * a2
                - 6 * a4 * r
                + 3 * a4 * r2
                + 24 * a2 * r2 * r
                - 8 * a2 * r4
                - 8 * r4 * r2
                + 4 * a2 * (a4 + a2 * r2 - 6 * r2 * r) * cos2t
                + a4 * (a2 + r * (6 + r)) * cos4t
            )
            * tdot
            - 16
            * a
            * (
                (-r4 * (a2 + 3 * r2) + a4 * (a2 - r2) * cost4) * sint2
                - a2 * r4 * (np.square(sin2t))
            )
            * pdot
        )
    ) / (
        4
        * (r2 + a2 * cost2)
        * (
            2 * a2 * r2 * (2 + a2 - 2 * r + r2) * cost2
            + a4 * (a2 + (-2 + r) * r) * cost4
            + r2 * ((-2 + r) * r2 * r + a2 * (4 + r2) - 8 * a2 * cos2t)
        )
    )

    vals[5] = (1 / (np.power((r2 + a2 * cost2), 3))) * (
        (
            (
                (np.square((r2 + a2 * cost2)))
                * (r * (-a2 + r) + a2 * (-1 + r) * cost2)
                * (np.square(rdot))
            )
            / (a2 + (-2 + r) * r)
        )
        + ((a2 + (-2 + r) * r) * (-r2 + a2 * cost2) * (np.square(tdot)))
        + (2 * a2 * cost * (np.square((r2 + a2 * cost2))) * sint * rdot * thdot)
        + (r * (a2 + (-2 + r) * r) * (np.square((r2 + a2 * cost2))) * (np.square(thdot)))
        - (4 * a * (a2 + (-2 + r) * r) * (-r2 + a2 * cost2) * sint2 * tdot * pdot)
        + (
            (a2 + (-2 + r) * r)
            * (
                -a2 * r2
                + r4 * r
                + a2 * (a2 + r2 + 2 * r2 * r) * cost2
                + a4 * (-1 + r) * cost4
            )
            * sint2
            * (np.square(pdot))
        )
    )

    vals[6] = (1 / (np.power((r2 + a2 * cost2), 3))) * (
        (
            -(a2 * cost * (np.square((r2 + a2 * cost2))) * sint * (np.square(rdot)))
            / (a2 + (-2 + r) * r)
        )
        + (a2 * r * sin2t * (np.square(tdot)))
        - (2 * r * (np.square((r2 + a2 * cost2))) * rdot * thdot)
        + (a2 * cost * (np.square((r2 + a2 * cost2))) * sint * (np.square(thdot)))
        - (4 * a * r * (a2 + r2) * sin2t * tdot * pdot)
        - (
            cost
            * sint
            * (
                -(np.square((r2 + a2 * cost2)))
                * (a2 - 2 * r + r2 + ((2 * r * (a2 + r2)) / (r2 + a2 * cost2)))
                - 2 * a2 * r * (a2 + r2) * sint2
            )
            * (np.square(pdot))
        )
    )

    vals[7] = -(
        (
            2
            * (
                thdot
                * (
                    (
                        -2
                        * a
                        * r
                        * (a2 + (-2 + r) * r)
                        * (a2 + 2 * r2 + a2 * cos2t)
                        * (cost / sint)
                        * tdot
                    )
                    + (
                        (
                            (
                                a2
                                * r2
                                * (a2 * (-4 + 3 * r2) + r2 * (4 - 6 * r + 3 * r2))
                                * cost2
                            )
                            + (a4 * r2 * (4 + 3 * a2 - 6 * r + 3 * r2) * cost4)
                            + (a4 * a2 * (a2 + (-2 + r) * r) * cost4 * cost2)
                            + (
                                r2
                                * (
                                    (-2 + r) * r4 * r
                                    + a4 * (6 + r)
                                    + a2 * r2 * (2 + r + r2)
                                    - a2 * (6 + r) * (a2 + r2) * cos2t
                                )
                            )
                        )
                        * (cost / sint)
                        + (2 * a4 * r * (a2 + r2) * cost2 * cost * sint)
                    )
                    * pdot
                )
                + (
                    rdot
                    * (
                        ((2 * a * r4 - 2 * a4 * a * cost4) * tdot)
                        + (
                            2 * a2 * r2 * r
                            - a2 * r4
                            - 2 * r4 * r2
                            + r4 * r2 * r
                            - a2 * r * (2 * a2 + r2 * (2 + 3 * r - 3 * r2)) * cost2
                            + a4 * (a2 + r * (2 - 2 * r + 3 * r2)) * cost4
                            + a4 * a2 * (-1 + r) * cost4 * cost2
                            - 8 * a2 * r2 * r * sint2
                            + 2 * a4 * r * (np.square(sin2t))
                        )
                        * pdot
                    )
                )
            )
        )
        / (
            (r2 + a2 * cost2)
            * (
                2 * a2 * r2 * (2 + a2 - 2 * r + r2) * cost2
                + a4 * (a2 + (-2 + r) * r) * cost4
                + r2 * ((-2 + r) * r2 * r + a2 * (4 + r2) - 8 * a2 * cos2t)
            )
        )
    )

    return vals
