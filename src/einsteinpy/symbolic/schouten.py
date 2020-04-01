from einsteinpy.symbolic.helpers import _change_name
from einsteinpy.symbolic.ricci import RicciScalar, RicciTensor
from einsteinpy.symbolic.riemann import RiemannCurvatureTensor
from einsteinpy.symbolic.tensor import BaseRelativityTensor, _change_config


class SchoutenTensor(BaseRelativityTensor):
    """
    Class for defining Schouten Tensor

    """

    def __init__(
        self, arr, syms, config="ll", parent_metric=None, name="SchoutenTensor"
    ):
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
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor
            Corresponding Metric for the Schouten Tensor. Defaults to None.
        name : str
            Name of the Tensor. Defaults to "SchoutenTensor".

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 2 indices

        """
        super(SchoutenTensor, self).__init__(
            arr=arr, syms=syms, config=config, parent_metric=parent_metric, name=name
        )
        self._order = 2
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))

    @classmethod
    def from_metric(cls, metric):
        """
        Get Schouten tensor calculated from a metric tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Space-time Metric from which Christoffel Symbols are to be calculated

        Raises
        ------
        ValueError
            Raised when the dimension of the tensor is less than 3

        """
        if metric.dims >= 3:
            t_ricci = RicciTensor.from_metric(metric)
            r_scalar = RicciScalar.from_riccitensor(t_ricci, parent_metric=None)
            dims = metric.dims
            t_schouten = (
                t_ricci.tensor()
                - (r_scalar.expr * metric.lower_config().tensor() / (2 * (dims - 1)))
            ) / (dims - 2)
            return cls(t_schouten, metric.syms, config="ll", parent_metric=metric)
        raise ValueError("Dimension of the space/space-time should be 3 or more")

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
            Compulsory if not initialized with 'from_metric'. Defaults to None.

        Returns
        -------
        ~einsteinpy.symbolic.schouten.SchoutenTensor
            New tensor with new configuration. Configuration defaults to 'ul'

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
        new_obj = SchoutenTensor(
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
            ~einsteinpy.symbolic.schouten.SchoutenTensor
                lorentz transformed tensor

        """
        t = super(SchoutenTensor, self).lorentz_transform(transformation_matrix)
        return SchoutenTensor(
            t.tensor(),
            syms=self.syms,
            config=self._config,
            parent_metric=None,
            name=_change_name(self.name, context="__lt"),
        )
