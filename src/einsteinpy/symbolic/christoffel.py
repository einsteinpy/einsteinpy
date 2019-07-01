import numpy as np
import sympy

from einsteinpy.symbolic.tensor import Tensor


class ChristoffelSymbols(Tensor):
    """
    Class for defining christoffel symbols
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
        super(ChristoffelSymbols, self).__init__(arr)
        if isinstance(syms, (list, tuple)):
            self.syms = syms
            self.dims = len(self.syms)
        else:
            raise TypeError("syms should be a list or tuple")

    @classmethod
    def from_metric(cls, metric):
        """
        Get Christoffel symbols calculated from a metric tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Space-time Metric from which Christoffel Symbols are to be calculated
        
        """
        dims = metric.dims
        tmplist = np.zeros((dims, dims, dims), dtype=int).tolist()
        mat, syms = metric.tensor(), metric.symbols()
        matinv = sympy.Matrix(mat.tolist()).inv()
        for t in range(dims ** 3):
            # i,j,k each goes from 0 to (dims-1)
            # hack for codeclimate. Could be done with 3 nested for loops
            k = t % dims
            j = (int(t / dims)) % (dims)
            i = (int(t / (dims ** 2))) % (dims)
            tmpvar = 0
            for n in range(dims):
                tmpvar += (matinv[i, n] / 2) * (
                    sympy.diff(mat[n, j], syms[k])
                    + sympy.diff(mat[n, k], syms[j])
                    - sympy.diff(mat[j, k], syms[n])
                )
            tmplist[i][j][k] = tmpvar
        return cls(tmplist, syms)

    def symbols(self):
        """
        Returns the symbols used for defining the time & spacial axis

        Returns
        -------
        tuple
            tuple containing (t,x1,x2,x3)
        
        """
        return self.syms
