import sympy

from einsteinpy.symbolic.tensor import BaseRelativityTensor


class MetricTensor(BaseRelativityTensor):
    """
    Inherits from ~einsteinpy.symbolic.tensor.BaseRelativityTensor .
    Class to define a metric tensor for a space-time
    """

    def __init__(self, arr, syms, config="ll"):
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
            arr=arr, syms=syms, config=config, parent_metric=self
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
            New Metric with new configuration. Defaults to 'uu'

        Raises
        ------
        ValueError
            Raised when new configuration is not 'll' or 'uu'. 
            This constraint is in place because we are dealing with Metric Tensor.

        """
        if newconfig == self.config:
            return self
        elif newconfig == "uu" or newconfig == "ll":
            inv_met = MetricTensor(
                sympy.simplify(sympy.Matrix(self.arr.tolist()).inv()).tolist(),
                self.syms,
                config=newconfig,
            )
            inv_met._invmetric = self
            return inv_met
        else:
            raise ValueError(
                "Configuration can't have one upper and one lower index in Metric Tensor."
            )

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
