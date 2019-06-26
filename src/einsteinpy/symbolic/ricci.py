import numpy as np
import sympy

from .metric import MetricTensor
from .tensor import Tensor


class RicciTensor(Tensor):
    """
    Class for defining Ricci Tensor
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
        super(RicciTensor, self).__init__(arr)
        if isinstance(syms, (list, tuple)):
            self.syms = syms
            self.dims = len(self.syms)
        else:
            raise TypeError("syms should be a list or tuple")

    @classmethod
    def from_christoffels(cls, chris):
        """
        Get Ricci Tensor calculated from a Christoffel Symbols

        Parameters
        ----------
        chris : ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            Christoffel Symbols from which Ricci Tensor to be calculated

        """
        arr, syms = chris.tensor(), chris.symbols()
        dims = len(syms)
        ricci_list = (np.zeros(shape=(dims, dims), dtype=int)).tolist()
        for i in range(dims ** 2):
            # r,t each goes from 0 to (dims-1)
            # hack for codeclimate. Could be done with 2 nested for loops
            r = (int(i / dims)) % (dims)
            t = (int(i / (dims ** 1))) % (dims)

            temp = sympy.diff(arr[t], syms[r])
            for p in range(dims):
                temp += arr[p, r, p, t]
            ricci_list[t][r] = sympy.simplify(temp)
        return cls(ricci_list, syms)

    @classmethod
    def from_metric(cls, metric):
        """
        Get Ricci Tensor calculated from a Metric Tensor

        Parameters
        ----------
        metric : ~~einsteinpy.symbolic.metric.MetricTensor
            Metric Tensor

        """
        rc = RicciTensor.from_metric(metric)
        return cls.from_christoffels(rc)

    @classmethod
    def from_riemann(cls, riemann):
        """
        Get Ricci Tensor calculated from a Metric Tensor

        Parameters
        ----------
        metric : ~~einsteinpy.symbolic.metric.MetricTensor
           Metric Tensor

        """
        return cls(sympy.tensorcontraction(riemann.tensor(), (0, 2)), riemann.syms)

    def symbols(self):
        """
        Returns the symbols used for defining the time & spacial axis

        Returns
        -------
        tuple
            tuple containing (x1,x2)

        """
        return self.syms
