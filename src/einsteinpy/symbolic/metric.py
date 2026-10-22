import sympy

from einsteinpy.symbolic.helpers import _change_name
from einsteinpy.symbolic.tensor import BaseRelativityTensor


class MetricTensor(BaseRelativityTensor):
    """
    Class to define a metric tensor for a space-time
    """

    def __init__(self, arr, syms, config="ll", name="GenericMetricTensor"):
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
        name : str
            Name of the Metric Tensor. Defaults to "GenericMetricTensor".

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 2 indices

        """
        super(MetricTensor, self).__init__(
            arr=arr, syms=syms, config=config, parent_metric=self, name=name
        )
        self._order = 2
        self._invmetric = None
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))

    def change_config(self, newconfig="uu"):
        """
        Changes the index configuration(contravariant/covariant)

        Parameters
        ----------
        newconfig : str
            Specify the new configuration. Defaults to 'uu'

        Returns
        -------
        ~einsteinpy.symbolic.metric.MetricTensor
            New Metric with new configuration.

        Raises
        ------
        ValueError
            Raised when the length of new configuration is not 2.

        """
        if newconfig == self.config:
            return self
            
        if len(newconfig) != 2:
            raise ValueError("Configuration should have exactly 2 indices for MetricTensor")
            
        # For 'uu' or 'll' configurations, we can use the inverse matrix
        if (newconfig == "uu" and self.config == "ll") or (newconfig == "ll" and self.config == "uu"):
            inv_met = MetricTensor(
                sympy.simplify(sympy.Matrix(self.arr.tolist()).inv()).tolist(),
                self.syms,
                config=newconfig,
                name=_change_name(self.name, context="__" + newconfig),
            )
            inv_met._invmetric = self
            return inv_met
            
        # For mixed indices (like 'ul' or 'lu'), we handle differently
        # For a metric, this is equivalent to the Kronecker delta or identity matrix
        if (newconfig == "ul" or newconfig == "lu"):
            # Create an identity matrix of the same dimension as the metric
            dim = len(self.syms)
            identity_list = [[0 for _ in range(dim)] for _ in range(dim)]
            for i in range(dim):
                identity_list[i][i] = 1
                
            mixed_met = MetricTensor(
                identity_list,
                self.syms,
                config=newconfig,
                name=_change_name(self.name, context="__" + newconfig),
            )
            return mixed_met
            
        # For other configurations (which should not occur for a rank-2 tensor)
        raise ValueError(f"Invalid configuration '{newconfig}' for MetricTensor")

    def inv(self):
        """
        Returns the inverse of the Metric.
        Returns contravariant Metric if it is originally covariant or vice-versa.

        Returns
        -------
        ~einsteinpy.symbolic.metric.MetricTensor
            New Metric which is the inverse of original Metric.

        """
        if self._invmetric is None:
            if self.config == "ll":
                self._invmetric = self.change_config("uu")
            else:
                self._invmetric = self.change_config("ll")
        return self._invmetric

    def lower_config(self):
        """
        Returns a covariant instance of the given metric tensor.

        Returns
        -------
        ~einsteinpy.symbolic.metric.MetricTensor
            same instance if the configuration is already lower or
            inverse of given metric if configuration is upper

        """
        if self.config == "ll":
            return self
        return self.inv()

    def lorentz_transform(self, transformation_matrix):
        """
        Performs a Lorentz transform on the tensor.

        Parameters
        ----------
            transformation_matrix : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
                Sympy Array or multi-dimensional list containing Sympy Expressions

        Returns
        -------
            ~einsteinpy.symbolic.metric.MetricTensor
                lorentz transformed tensor

        """
        t = super(MetricTensor, self).lorentz_transform(transformation_matrix)
        return MetricTensor(
            t.tensor(),
            syms=self.syms,
            config=self._config,
            name=_change_name(self.name, context="__lt"),
        )
