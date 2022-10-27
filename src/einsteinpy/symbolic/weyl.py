import numpy as np
import sympy

from einsteinpy.symbolic.helpers import _change_name
from einsteinpy.symbolic.ricci import RicciScalar, RicciTensor
from einsteinpy.symbolic.riemann import RiemannCurvatureTensor
from einsteinpy.symbolic.tensor import BaseRelativityTensor, _change_config


class WeylTensor(BaseRelativityTensor):
    """
    Class for defining Weyl Tensor
    """

    def __init__(self, arr, syms, config="ulll", parent_metric=None, name="WeylTensor"):
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
            arr=arr, syms=syms, config=config, parent_metric=parent_metric, name=name
        )
        self._order = 4
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))

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

    def change_config(self, newconfig="llll", metric=None):
        """
        Changes the index configuration(contravariant/covariant)

        Parameters
        ----------
        newconfig : str
            Specify the new configuration. Defaults to 'llll'
        metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Parent metric tensor for changing indices.
            Already assumes the value of the metric tensor from which it was initialized if passed with None.
            Compulsory if not initialized with 'from_metric'. Defaults to None.

        Returns
        -------
        ~einsteinpy.symbolic.weyl.WeylTensor
            New tensor with new configuration. Configuration defaults to 'llll'

        Raises
        ------
        Exception
            Raised when a parent metric could not be found.

        """
        if metric is None:
            metric = self._parent_metric
        if metric is None:
            raise Exception("Parent Metric not found, can't do configuration change")
        new_tensor = _change_config(self, metric, newconfig)
        new_obj = WeylTensor(
            new_tensor,
            self.syms,
            config=newconfig,
            parent_metric=metric,
            name=_change_name(self.name, context="__" + newconfig),
        )
        return new_obj

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
