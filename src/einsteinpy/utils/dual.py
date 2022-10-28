import numpy as np


class DualNumber:
    """
    Numbers of the form, :math:`a + b\\epsilon`, where
    :math:`\\epsilon^2 = 0` and :math:`\\epsilon \\ne 0`.
    Their addition and multiplication properties make them
    suitable for Automatic Differentiation (AD).
    EinsteinPy uses AD for solving Geodesics in arbitrary spacetimes.
    This module is based on [1]_.

    References
    ----------
    .. [1] Christian, Pierre and Chan, Chi-Kwan;
        "FANTASY: User-Friendly Symplectic Geodesic Integrator
        for Arbitrary Metrics with Automatic Differentiation";
        `2021 ApJ 909 67 <https://doi.org/10.3847/1538-4357/abdc28>`__

    """

    def __init__(self, val, deriv):
        """
        Constructor

        Parameters
        ----------
        val : float
            Value
        deriv : float
            Directional Derivative

        """
        self.val = float(val)
        self.deriv = float(deriv)

    def __str__(self):
        return f"DualNumber({self.val}, {self.deriv})"

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        if isinstance(other, DualNumber):
            return DualNumber(self.val + other.val, self.deriv + other.deriv)

        return DualNumber(self.val + other, self.deriv)

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, DualNumber):
            return DualNumber(self.val - other.val, self.deriv - other.deriv)

        return DualNumber(self.val - other, self.deriv)

    def __rsub__(self, other):
        if isinstance(other, DualNumber):
            return DualNumber(other.val - self.val, other.deriv - self.deriv)

        return DualNumber(other, 0) - self

    def __mul__(self, other):
        if isinstance(other, DualNumber):
            return DualNumber(
                self.val * other.val, self.deriv * other.val + self.val * other.deriv
            )

        return DualNumber(self.val * other, self.deriv * other)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, DualNumber):
            if self.val == 0 and other.val == 0:
                return DualNumber(self.deriv / other.deriv, 0.0)

            return DualNumber(
                self.val / other.val,
                (self.deriv * other.val - self.val * other.deriv) / (other.val**2),
            )

        return DualNumber(self.val / other, self.deriv / other)

    def __rtruediv__(self, other):
        if isinstance(other, DualNumber):
            if self.val == 0 and other.val == 0:
                return DualNumber(other.deriv / self.deriv, 0.0)

            return DualNumber(
                other.val / self.val,
                (other.deriv * self.val - other.val * self.deriv) / (self.val**2),
            )

        return DualNumber(other, 0).__truediv__(self)

    def __eq__(self, other):
        return (self.val == other.val) and (self.deriv == other.deriv)

    def __ne__(self, other):
        return not (self == other)

    def __neg__(self):
        return DualNumber(-self.val, -self.deriv)

    def __pow__(self, power):
        return DualNumber(
            self.val**power, self.deriv * power * self.val ** (power - 1)
        )

    def sin(self):
        return DualNumber(np.sin(self.val), self.deriv * np.cos(self.val))

    def cos(self):
        return DualNumber(np.cos(self.val), -self.deriv * np.sin(self.val))

    def tan(self):
        return np.sin(self) / np.cos(self)

    def log(self):
        return DualNumber(np.log(self.val), self.deriv / self.val)

    def exp(self):
        return DualNumber(np.exp(self.val), self.deriv * np.exp(self.val))


def _deriv(func, x):
    """
    Calculates first (partial) derivative of ``func`` at ``x``

    Parameters
    ----------
    func : callable
        Function to differentiate
    x : float
        Point, at which, the derivative will be evaluated

    Returns
    _______
    float
        First partial derivative of ``func`` at ``x``

    """
    funcdual = func(DualNumber(x, 1.0))

    if isinstance(funcdual, DualNumber):
        return funcdual.deriv

    return 0.0


def _diff_g(g, g_prms, coords, indices, wrt):
    """
    Computes derivative of metric elements

    Parameters
    ----------
    g : callable
        Metric (Contravariant) Function
    g_prms : array_like
        Tuple of parameters to pass to the metric
        E.g., ``(a,)`` for Kerr
    coords : array_like
        4-Position
    indices : array_like
        2-tuple, containing indices, indexing a metric
        element, whose derivative will be calculated
    wrt : int
        Coordinate, with respect to which, the derivative
        will be calculated
        Takes values from ``[0, 1, 2, 3]``

    Returns
    -------
    float
        Value of derivative of metric element at ``coords``

    Raises
    ------
    ValueError
        If ``wrt`` is not in [1, 2, 3, 4]
        or ``len(indices) != 2``

    """
    if wrt not in [0, 1, 2, 3]:
        raise ValueError(f"wrt takes values from [0, 1, 2, 3]. Supplied value: {wrt}")

    if len(indices) != 2:
        raise ValueError("indices must be a 2-tuple containing indices for the metric.")

    dual_coords = [
        DualNumber(coords[0], 0.0),
        DualNumber(coords[1], 0.0),
        DualNumber(coords[2], 0.0),
        DualNumber(coords[3], 0.0),
    ]

    # Coordinate, against which, derivative will be propagated
    dual_coords[wrt].deriv = 1.0

    return _deriv(lambda q: g(dual_coords, *g_prms)[indices], coords[wrt])


def _jacobian_g(g, g_prms, coords, wrt):
    """
    Part of Jacobian of Metric

    Parameters
    ----------
    g : callable
        Metric (Contravariant) Function
    g_prms : array_like
        Tuple of parameters to pass to the metric
        E.g., ``(a,)`` for Kerr
    coords : array_like
        4-Position
    wrt : int
        Coordinate, with respect to which, the derivative
        will be calculated
        Takes values from ``[0, 1, 2, 3]``

    Returns
    -------
    numpy.ndarray
        Value of derivative of metric elements,
        w.r.t a particular coordinate, at ``coords``

    """
    J = np.zeros((4, 4))

    for i in range(4):
        for j in range(4):
            if i <= j:
                J[i, j] = _diff_g(g, g_prms, coords, (i, j), wrt)

    J = J + J.T - np.diag(np.diag(J))

    return J
