import numpy as np


class Dual:
    """
    Numbers of the form, :math:`a + b\\epsilon`, where
    :math:`\\epsilon^2 = 0` and :math:`\\epsilon \\ne 0`.
    Their addition and multiplication properties make them
    suitable for Automatic Differentiation (AD).
    EinsteinPy uses AD for solving Geodesics in arbitrary spacetimes.

    """

    def __init__(self, a, b):
        """
        Constructor

        Parameters
        ----------
        a : float
            First component
        b : float
            Second component

        """
        self.a = float(a)
        self.b = float(b)

    def __str__(self):
        return f"Dual({self.a}, {self.b})"

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        if isinstance(other, Dual):
            return Dual(self.a + other.a, self.b + other.b)

        return Dual(self.a + other, self.b)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, Dual):
            return Dual(self.a - other.a, self.b - other.b)

        return Dual(self.a - other, self.b)

    def __rsub__(self, other):
        if isinstance(other, Dual):
            return Dual(other.a - self.a, other.b - self.b)

        return Dual(other, 0) - self

    def __mul__(self, other):
        if isinstance(other, Dual):
            return Dual(self.a * other.a, self.b * other.a + self.a * other.b)

        return Dual(self.a * other, self.b * other)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, Dual):
            if self.a == 0 and other.a == 0:
                return Dual(self.b / other.b, 0.0)

            return Dual(
                self.a / other.a, (self.b * other.a - self.a * other.b) / (other.a ** 2)
            )

        return Dual(self.a / other, self.b / other)

    def __rtruediv__(self, other):
        if isinstance(other, Dual):
            if self.a == 0 and other.a == 0:
                return Dual(other.b / self.b, 0.0)

            return Dual(
                other.a / self.a, (other.b * self.a - other.a * self.b) / (self.a ** 2)
            )

        return Dual(other, 0).__truediv__(self)

    def __eq__(self, other):
        return (self.a == other.a) and (self.b == other.b)

    def __ne__(self, other):
        return not (self == other)

    def __neg__(self):
        return Dual(-self.a, -self.b)

    def __pow__(self, power):
        return Dual(self.a ** power, self.b * power * self.a ** (power - 1))

    def sin(self):
        return Dual(np.sin(self.a), self.b * np.cos(self.a))

    def cos(self):
        return Dual(np.cos(self.a), -self.b * np.sin(self.a))

    def tan(self):
        return np.sin(self) / np.cos(self)

    def log(self):
        return Dual(np.log(self.a), self.b / self.a)

    def exp(self):
        return Dual(np.exp(self.a), self.b * np.exp(self.a))


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
    funcdual = func(Dual(x, 1.0))

    if isinstance(funcdual, Dual):
        return funcdual.b

    return 0.0


def _diff_g(g, g_prms, coords, indices, wrt):
    """
    Computes derivative of metric elements

    Parameters
    ----------
    g : callable
        Metric (Contravariant) Function
    coords : array_like
        4-Position
    indices : array_like
        2-tuple, containing indices, indexing a metric
        element, whose derivative will be calculated
    wrt : int
        coordsinate, with respect to which, the derivative
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
        Dual(coords[0], 0.0),
        Dual(coords[1], 0.0),
        Dual(coords[2], 0.0),
        Dual(coords[3], 0.0),
    ]

    dual_coords[wrt].b = 1.0

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
        coordsinate, with respect to which, the derivative
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
