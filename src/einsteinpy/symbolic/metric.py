from einsteinpy.symbolic.tensor import Tensor


class MetricTensor(Tensor):
    """
    Class to define a metric tensor for a space-time
    """

    def __init__(self, arr, syms):
        """
        Constructor and Initializer
        
        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        syms : tuple or list
            Tuple of crucial symbols dentoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        
        """
        super(MetricTensor, self).__init__(arr)
        if isinstance(syms, (list, tuple)):
            self.syms = syms
            self.dims = len(self.syms)
        else:
            raise TypeError("syms should be a list or tuple")

    def symbols(self):
        """
        Returns the symbols used for defining the time & spacial axis

        Returns
        -------
        tuple
            tuple containing (t,x1,x2,x3)
        
        """
        return self.syms
