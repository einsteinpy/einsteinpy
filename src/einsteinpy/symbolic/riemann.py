import numpy as np
import sympy

from .christoffel import ChristoffelSymbols
from .metric import MetricTensor
from .tensor import Tensor


class RiemannCurvatureTensor(Tensor):
    """
    Class for defining Riemann Curvature Tensor
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
        super(RiemannCurvatureTensor, self).__init__(arr)
        if isinstance(syms, (list, tuple)):
            self.syms = syms
            self.dims = len(self.syms)
        else:
            raise TypeError("syms should be a list or tuple")

    @classmethod
    def from_christoffels(cls, chris):
        """
        Get Riemann Tensor calculated from a Christoffel Symbols

        Parameters
        ----------
        chris : ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            Christoffel Symbols from which Riemann Curvature Tensor to be calculated
        
        """
        arr, syms = chris.tensor(), chris.symbols()
        dims = len(syms)
        riemann_list = (np.zeros(shape=(dims, dims, dims, dims), dtype=int)).tolist()
        for i in range(dims ** 4):
            # t,s,r,n each goes from 0 to (dims-1)
            # hack for codeclimate. Could be done with 4 nested for loops
            n = i % dims
            r = (int(i / dims)) % (dims)
            s = (int(i / (dims ** 2))) % (dims)
            t = (int(i / (dims ** 3))) % (dims)
            temp = sympy.diff(arr[t, s, n], syms[r]) - sympy.diff(arr[t, r, n], syms[s])
            for p in range(dims):
                temp += arr[p, s, n] * arr[t, p, r] - arr[p, r, n] * arr[t, p, s]
            riemann_list[t][s][r][n] = sympy.simplify(temp)
        return cls(riemann_list, syms)

    @classmethod
    def from_metric(cls, metric):
        """
        Get Riemann Tensor calculated from a Metric Tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Metric Tensor from which Riemann Curvature Tensor to be calculated
        
        """
        ch = ChristoffelSymbols.from_metric(metric)
        return cls.from_christoffels(ch)

    def symbols(self):
        """
        Returns the symbols used for defining the time & spacial axis

        Returns
        -------
        tuple
            tuple containing (t,x1,x2,x3)
        
        """
        return self.syms
