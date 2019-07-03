from collections import namedtuple

from sympy import Array, ones, tensorproduct, zeros
from sympy.tensor.tensor import TensorIndexType

from .partial import PartialDerivative
from .tensor import AbstractTensor, Tensor, expand_tensor, indices


class Metric(AbstractTensor, TensorIndexType):
    """
    Class representing a tensor that raises and lowers indices.
    """

    # This object allows for having the metric be represented by the
    # same symbol as the tensor it associates with.
    _MetricId = namedtuple("MetricId", ["name", "antisym"])
    is_Metric = True
    _christoffel = None
    _riemann = None
    _ricci_tensor = None
    _ricci_scalar = None
    _weyl = None
    _einstein = None

    def __new__(cls, symbol, coords, matrix, **kwargs):
        """
        Create a new Metric object.

        Parameters
        ----------
        symbol : str
            Name of the tensor and the symbol to denote it by when printed.
        coords : iterable
            List of ~sympy.Symbol objects to denote the coordinates by which
            derivatives are taken with respect to.
        matrix : (list, tuple, ~sympy.Matrix, ~sympy.Array)
            Matrix representation of the tensor to be used in substitution.
            Can be of any type that is acceptable by ~sympy.Array.

        Examples
        --------
        >>> from sympy import diag, symbols
        >>> from einsteinpy.symbolic.tensor import indices, expand_tensor
        >>> from einsteinpy.symbolic.metric import Metric
        >>> t, x, y, z = symbols('t x y z')
        >>> eta = Metric('eta', [t, x, y, z], diag(1, -1, -1, -1))
        >>> mu, nu = indices('mu nu', eta)
        >>> expr = eta(mu, nu) * eta(-mu, -nu)
        >>> expand_tensor(expr)
        4

        """
        array = Array(matrix)
        if array.rank() != 2 or array.shape[0] != array.shape[1]:
            raise ValueError(
                "matrix must be square, received matrix of shape {}".format(array.shape)
            )
        obj = TensorIndexType.__new__(
            cls,
            symbol,
            metric=cls._MetricId(symbol, False),
            dim=array.shape[0],
            dummy_fmt=symbol,
            **kwargs,
        )
        obj = AbstractTensor.__new__(cls, obj, array)
        obj.metric = Tensor(obj.name, array, obj, comm="metric", covar=(-1, -1))
        obj.coords = tuple(coords)
        return obj

    def __getattr__(self, attr):
        if hasattr(self.metric, attr):
            return getattr(self.metric, attr)
        return TensorIndexType.__getattribute__(self, attr)

    def __call__(self, *args):
        return self.metric(*args)

    @property
    def partial(self):
        return PartialDerivative(self)

    @property
    def christoffel(self):
        r"""
        Returns the Christoffel symbols using the formula:
        \Gamma^\rho_{\mu\nu} =
            \frac{1}{2} g^{\sigma\rho} (\partial_\mu g_{\nu\rho} + \partial_\nu g_{\rho\mu} - \partial_\rho g_{\mu\nu})

        """
        if self._christoffel is None:
            mu, nu, si, rh = indices("mu nu sigma rho", self)
            d = self.partial
            g = self.metric
            gamma = (
                (1 / 2)
                * g(si, rh)
                * (d(-mu) * g(-nu, -rh) + d(-nu) * g(-rh, -mu) - d(-rh) * g(-mu, -nu))
            )
            syms = expand_tensor(gamma)
            self._christoffel = Tensor("Gamma", syms, self, covar=(1, -1, -1))
        return self._christoffel

    @property
    def riemann(self):
        r"""
        Returns the Riemann curvature tensor using the formula:
        R^\rho_{\sigma\mu\nu} =
            \partial_\mu \Gamma^\rho_{\nu\sigma} - \partial_\nu \Gamma^\rho_{\mu\sigma}
          + \Gamma^\rho_{\mu\lambda} \Gamma^\lambda_{\nu\sigma} - \Gamma^\rho_{\nu\lambda} \Gamma^\lambda_{\mu\sigma}

        """
        if self._riemann is None:
            mu, nu, si, rh, la = indices("mu nu sigma rho lambda", self)
            d = self.partial
            g = self.metric
            gamma = self.christoffel
            R = (
                d(-mu) * gamma(rh, -nu, -si)
                - d(-nu) * gamma(rh, -mu, -si)
                + gamma(rh, -mu, -la) * gamma(la, -nu, -si)
                - gamma(rh, -nu, -la) * gamma(la, -mu, -si)
            )
            res = expand_tensor(R)
            self._riemann = Tensor(
                "R", res, self, symmetry=[[2, 2]], covar=(1, -1, -1, -1)
            )
        return self._riemann

    @property
    def ricci_tensor(self):
        r"""
        Returns the Ricci tensor using the formula:
        R_{\mu\nu} = R^\sigma_{\mu\sigma\nu}

        """
        if self._ricci_tensor is None:
            mu, nu, si = indices("mu nu sigma", self)
            R = self.riemann
            res = expand_tensor(R(si, -mu, -si, -nu))
            self._ricci_tensor = Tensor("R", res, self, covar=(-1, -1))
        return self._ricci_tensor

    @property
    def ricci_scalar(self):
        r"""
        Returns the Ricci scalar using the formula:
        R = R^\mu_\mu

        """
        if self._ricci_scalar is None:
            mu, nu = indices("mu nu", self)
            RR = self.ricci_tensor
            res = expand_tensor(RR(mu, -mu))
            self._ricci_scalar = res
        return self._ricci_scalar

    @property
    def weyl(self):
        r"""
        Returns the Weyl conformal tensor using the formula:
        C_{\rho\sigma\mu\nu} =
            R_{\rho\sigma\mu\nu} - \frac{2}{(n - 2)} (g_{\rho[\mu} R_{\nu]\sigma} - g_{\sigma[\mu} R_{\nu]\rho})
          + \frac{2}{(n - 1)(n - 2)} g_{\rho[\mu} g_{\nu]\sigma} R

        """
        if self._weyl is None:
            n = self.dim
            if n < 3:
                raise ValueError(
                    "the Weyl tensor is only defined in dimensions of 3 or more. {} is of dimension {}".format(
                        self, n
                    )
                )
            elif n == 3:
                res = tensorproduct(zeros(3, 3), zeros(3, 3))
                self._weyl = Tensor(
                    "C", res, self, symmetry=[[2, 2]], covar=(1, -1, -1, -1)
                )
                return self._weyl
            c1 = 1 / (n - 2)
            c2 = 1 / (n - 2) / (n - 1)
            mu, nu, si, rh = indices("mu nu sigma rho", self)
            R = self.riemann
            RR = self.ricci_tensor
            RRR = self.ricci_scalar
            g = self.metric
            C = (
                R(rh, -si, -mu, -nu)
                - c1
                * (
                    g(rh, -mu) * RR(-nu, -si)
                    - g(rh, -nu) * RR(-mu, -si)
                    + g(-si, -nu) * RR(-mu, rh)
                    - g(-si, -mu) * RR(-nu, rh)
                )
                + c2 * (g(rh, -mu) * g(-nu, -si) - g(rh, -nu) * g(-mu, -si)) * RRR
            )
            res = expand_tensor(C)
            self._weyl = Tensor(
                "C", res, self, symmetry=[[2, 2]], covar=(1, -1, -1, -1)
            )
        return self._weyl

    @property
    def einstein(self):
        r"""
        Returns the Einstein tensor using the formula:
        G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu}

        """
        if self._einstein is None:
            mu, nu = indices("mu nu", self)
            g = self.metric
            R = self.ricci_tensor
            RR = self.ricci_scalar
            res = expand_tensor(R(-mu, -nu) - (1 / 2) * RR * g(-mu, -nu))
            self._einstein = Tensor("G", res, self, covar=(-1, -1))
        return self._einstein


class SpacetimeMetric(Metric):
    """
    Class representing psuedo-Riemannian metrics.
    """

    is_Spacetime = True

    def __new__(cls, symbol, coords, matrix, timelike=True, **kwargs):
        obj = super().__new__(cls, symbol, coords, matrix, **kwargs)
        if obj.dim > 4:
            raise ValueError("metrics on spacetime must be at most 4-dimensional")
        obj.is_timelike = timelike
        obj.is_spacelike = not timelike
        return obj

    def reverse_signature(self):
        self._array *= -1
        self._replacement_dict = {self: self._array}
        self.is_timelike = not self.is_timelike
        self.is_spacelike = not self.is_spacelike
        return self.signature

    @property
    def signature(self):
        sign = -1 if self.is_timelike else 1
        sig = sign * ones(1, self.dim)
        sig[0] *= -1
        return tuple(sig)
