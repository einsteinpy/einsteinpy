def v1(g_cov_mat, v0, v2, v3):
    """
        OUTDATED - ?????
        Utility function to return first space-like of 4-Velocity
        Assumes a (+, -, -, -) Metric Signature

        Parameters
        ----------
        g_cov_mat : ~numpy.ndarray
            Matrix, containing Covariant Metric \
            Tensor values, in same coordinates as ``v_vec``
            Numpy array of shape (4,4)
        v0 : float
            First component of 4-Velocity
        v2 : float
            Third component of 4-Velocity
        v3 : float
            Fourth component of 4-Velocity

        Returns
        -------
        float
            First spacelike component of 4-Velocity/4-Momentum

        """
    g = g_cov_mat
    # fac = 0, so removed from C
    # Defining coefficients for quadratic equation
    A = g[1, 1]
    B = 2 * (g[0, 1] * v0 + g[1, 2] * v2 + g[1, 3] * v3)
    C = (
        (g[0, 0] * v0 ** 2 + g[2, 2] * v2 ** 2 + g[3, 3] * v3 ** 2)
        + 2 * v0 * (g[0, 2] * v2 + g[0, 3] * v3)
        + 2 * v2 * g[2, 3] * v3
    )
    D = (B ** 2) - (4 * A * C)

    v1 = (-B + np.sqrt(D)) / (2 * A)

    return v1
