import numpy as np
import sympy
from sympy import simplify, tensorcontraction, tensorproduct

from einsteinpy.symbolic.helpers import _change_name
from einsteinpy.symbolic.ricci import RicciScalar, RicciTensor
from einsteinpy.symbolic.riemann import RiemannCurvatureTensor
from einsteinpy.symbolic.tensor import BaseRelativityTensor, _change_config, tensor_product


class WeylTensor(BaseRelativityTensor):
    """
    Class for defining Weyl Tensor
    """

    def __init__(self, arr, syms, config="ulll", parent_metric=None, parent_spacetime=None, name="WeylTensor", simplify=True):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'ulll'.
        parent_metric : ~einsteinpy.symbolic.metric.WeylTensor
            Corresponding Metric for the Weyl Tensor. Defaults to None.
        name : str
            Name of the Tensor. Defaults to "WeylTensor"

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 4 indices

        """
        super(WeylTensor, self).__init__(
            arr=arr, syms=syms, config=config, parent_metric=parent_metric, parent_spacetime=parent_spacetime, name=name, simplify=simplify
        )
        self._order = 4
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))
        
        self._dual = None

    @classmethod
    def from_metric(cls, metric):
        """
        Get Weyl tensor calculated from a metric tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Space-time Metric from which Christoffel Symbols are to be calculated

        Raises
        ------
        ValueError
            Raised when the dimension of the tensor is less than 3

        """
        if metric.dims > 3:
            metric_cov = metric.lower_config()
            t_riemann = RiemannCurvatureTensor.from_metric(metric)
            # Riemann Tensor with covariant indices is needed
            t_riemann_cov = t_riemann.change_config("llll", metric=None)
            t_ricci = RicciTensor.from_riemann(t_riemann, parent_metric=None)
            r_scalar = RicciScalar.from_riccitensor(t_ricci, parent_metric=None)
            g = metric_cov
            dims = g.dims
            # Indexing for resultant Weyl Tensor is iklm
            C = np.zeros(shape=(dims, dims, dims, dims), dtype=int).tolist()
            for t in range(dims**4):
                i, k, l, m = (
                    t % dims,
                    (int(t / dims)) % (dims),
                    (int(t / (dims**2))) % (dims),
                    (int(t / (dims**3))) % (dims),
                )
                C[i][k][l][m] = t_riemann_cov[i, k, l, m] + (
                    (
                        (
                            t_ricci[i, m] * g[k, l]
                            - t_ricci[i, l] * g[k, m]
                            + t_ricci[k, l] * g[i, m]
                            - t_ricci[k, m] * g[i, l]
                        )
                        / (dims - 2)
                    )
                    + (
                        r_scalar.expr
                        * (g[i, l] * g[k, m] - g[i, m] * g[k, l])
                        / ((dims - 1) * (dims - 2))
                    )
                )
            C = sympy.simplify(sympy.Array(C))
            return cls(C, metric.syms, config="llll", parent_metric=metric)
        if metric.dims == 3:
            return cls(
                sympy.Array(np.zeros((3, 3, 3, 3), dtype=int)),
                metric.syms,
                config="llll",
                parent_metric=metric,
            )
        raise ValueError("Dimension of the space/space-time should be 3 or more")

    @classmethod
    def from_tensors(cls, metric, riemann, ricci, ricci_scalar):
        """
        Get Weyl tensor calculated from the metric, riemann, ricci tensor and ricci scalar

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Space-time Metric from which Christoffel Symbols are to be calculated

        Raises
        ------
        ValueError
            Raised when the dimension of the tensor is less than 3

        """
        if metric.dims > 3:
            metric_cov = metric.lower_config()
            t_riemann  = riemann
            t_riemann_cov = t_riemann.change_config("llll", metric=None)
            t_ricci = ricci
            r_scalar = ricci_scalar
            g = metric_cov
            dims = g.dims
            # Indexing for resultant Weyl Tensor is iklm
            C = np.zeros(shape=(dims, dims, dims, dims), dtype=int).tolist()
            for t in range(dims**4):
                i, k, l, m = (
                    t % dims,
                    (int(t / dims)) % (dims),
                    (int(t / (dims**2))) % (dims),
                    (int(t / (dims**3))) % (dims),
                )
                C[i][k][l][m] = t_riemann_cov[i, k, l, m] + (
                    (
                        (
                            t_ricci[i, m] * g[k, l]
                            - t_ricci[i, l] * g[k, m]
                            + t_ricci[k, l] * g[i, m]
                            - t_ricci[k, m] * g[i, l]
                        )
                        / (dims - 2)
                    )
                    + (
                        r_scalar.expr
                        * (g[i, l] * g[k, m] - g[i, m] * g[k, l])
                        / ((dims - 1) * (dims - 2))
                    )
                )
            C = sympy.simplify(sympy.Array(C))
            return cls(C, metric.syms, config="llll", parent_metric=metric)
        if metric.dims == 3:
            return cls(
                sympy.Array(np.zeros((3, 3, 3, 3), dtype=int)),
                metric.syms,
                config="llll",
                parent_metric=metric,
            )
        raise ValueError("Dimension of the space/space-time should be 3 or more")



    def lorentz_transform(self, transformation_matrix):
        """
        Performs a Lorentz transform on the tensor.

        Parameters
        ----------
            transformation_matrix : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
                Sympy Array or multi-dimensional list containing Sympy Expressions

        Returns
        -------
            ~einsteinpy.symbolic.weyl.WeylTensor
                lorentz transformed tensor(or vector)

        """
        t = super(WeylTensor, self).lorentz_transform(transformation_matrix)
        return WeylTensor(
            t.tensor(),
            syms=self.syms,
            config=self._config,
            parent_metric=None,
            name=_change_name(self.name, context="__lt"),
        )


class BelRobinsonTensor(BaseRelativityTensor):
    """
    Class for defining Bel-Robinson Tensor
    """

    def __init__(self, arr, syms, config="ulll", parent_metric=None, parent_spacetime=None, simplify=True, name="BelRobinsonTensor"):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'ulll'.
        parent_metric : ~einsteinpy.symbolic.metric.WeylTensor
            Corresponding Metric for the Weyl Tensor. Defaults to None.
        name : str
            Name of the Tensor. Defaults to "WeylTensor"

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 4 indices

        """
        super(BelRobinsonTensor, self).__init__(
            arr=arr, syms=syms, config=config, parent_metric=parent_metric, parent_spacetime=parent_spacetime, simplify=simplify, name=name
        )
        self._order = 4
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))

    @classmethod
    def from_weyl(cls, C):
        """
        Get Bel-Robinson tensor calculated from the Weyl tensor and the alternating Levi-Civita Tensor

        Parameters
        ----------
        C : ~einsteinpy.symbolic.WeylTensor
            WeylTensor

        Raises
        ------
        ValueError
            Raised when the dimension of the tensor is less than 3

        """
        C_dual = C.DualTensor

        l = tensorcontraction(tensor_product(C.change_config("llll"), C.change_config("ullu"), 0, 0).tensor(), (2, 5))
        r = tensorcontraction(tensor_product(C_dual.change_config("llll"), C_dual.change_config("ullu"), 0, 0).tensor(), (2, 5))

        return cls( (l + r) / 4, syms=C.parent_metric.syms, config="llll", parent_metric=C.parent_metric, simplify=False)
        


    def lorentz_transform(self, transformation_matrix):
        """
        Performs a Lorentz transform on the tensor.

        Parameters
        ----------
            transformation_matrix : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
                Sympy Array or multi-dimensional list containing Sympy Expressions

        Returns
        -------
            ~einsteinpy.symbolic.weyl.BelRobinsonTensor
                lorentz transformed tensor(or vector)

        """
        t = super(BelRobinsonTensor, self).lorentz_transform(transformation_matrix)
        return BelRobinsonTensor(
            t.tensor(),
            syms=self.syms,
            config=self._config,
            parent_metric=None,
            name=_change_name(self.name, context="__lt"),
        )
