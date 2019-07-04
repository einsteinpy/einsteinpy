import sympy

from einsteinpy.symbolic.tensor import Tensor


class MetricTensor(Tensor):
    """
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
            Tuple of crucial symbols dentoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
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
        super(MetricTensor, self).__init__(arr, config=config)
        self._order = 2
        self._invmetric = None
        if isinstance(syms, (list, tuple)):
            self.syms = syms
            self.dims = len(self.syms)
        else:
            raise TypeError("syms should be a list or tuple")
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

    def symbols(self):
        """
        Returns the symbols used for defining the time & spacial axis

        Returns
        -------
        tuple
            tuple containing (t,x1,x2,x3)
        
        """
        return self.syms
