import sympy

from einsteinpy.symbolic.ricci import RicciScalar, RicciTensor
from einsteinpy.symbolic.tensor import BaseRelativityTensor, _change_config


class EinsteinTensor(BaseRelativityTensor):
    """
    Inherits from ~einsteinpy.symbolic.tensor.BaseRelativityTensor .
    Class for defining Einstein Tensor
    """

    def __init__(self, arr, syms, config="ll", parent_metric=None):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'll'.
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Einstein Tensor.
            Defaults to None.

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 2 indices

        """
        super(EinsteinTensor, self).__init__(
            arr=arr, syms=syms, config=config, parent_metric=parent_metric
        )
        self._order = 2
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))

    @classmethod
    def from_metric(cls, metric):
        t_ricci = RicciTensor.from_metric(metric)
        r_scalar = RicciScalar.from_riccitensor(t_ricci, t_ricci.parent_metric)
        einstein_tensor = (
            t_ricci.tensor() - (1 / 2) * metric.lower_config().tensor() * r_scalar.expr
        )
        return cls(
            einstein_tensor,
            metric.syms,
            config="ll",
            parent_metric=t_ricci.parent_metric,
        )

    def change_config(self, newconfig="ul", metric=None):
        """
        Changes the index configuration(contravariant/covariant)

        Parameters
        ----------
        newconfig : str
            Specify the new configuration. Defaults to 'ul'
        metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Parent metric tensor for changing indices.
            Already assumes the value of the metric tensor from which it was initialized if passed with None.
            Compulsory if somehow does not have a parent metric. Defaults to None.

        Returns
        -------
        ~einsteinpy.symbolic.einstein.EinsteinTensor
            New tensor with new configuration. Defaults to 'ul'

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
        new_obj = EinsteinTensor(
            new_tensor, self.syms, config=newconfig, parent_metric=metric
        )
        return new_obj
