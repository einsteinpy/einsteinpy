import numpy as np
import sympy

from einsteinpy.symbolic.tensor import Tensor, _change_config


class ChristoffelSymbols(Tensor):
    """
    Class for defining christoffel symbols
    """

    def __init__(self, arr, syms, config="ull", parent_metric=None):
        """
        Constructor and Initializer
        
        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        syms : tuple or list
            Tuple of crucial symbols dentoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'ull'.
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor
            Metric Tensor from which Christoffel symbol is calculated. Defaults to None.

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 3 indices
        
        """
        super(ChristoffelSymbols, self).__init__(arr, config=config)
        self._order = 3
        self._parent_metric = parent_metric
        if isinstance(syms, (list, tuple)):
            self.syms = syms
            self.dims = len(self.syms)
        else:
            raise TypeError("syms should be a list or tuple")
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))

    @property
    def parent_metric(self):
        """
        Returns the Metric from which Christoffel Symbol was derived, if available.
        """
        return self._parent_metric

    @classmethod
    def from_metric(cls, metric):
        """
        Get Christoffel symbols calculated from a metric tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Space-time Metric from which Christoffel Symbols are to be calculated
        
        """
        if metric.config == "uu":
            metric_cov = metric.inv()
        else:
            metric_cov = metric
        dims = metric.dims
        tmplist = np.zeros((dims, dims, dims), dtype=int).tolist()
        mat, syms = metric_cov.tensor(), metric_cov.symbols()
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
        return cls(tmplist, syms, config="ull", parent_metric=metric)

    def change_config(self, newconfig="lll", metric=None):
        """
        Changes the index configuration(contravariant/covariant)

        Parameters
        ----------
        newconfig : str
            Specify the new configuration. Defaults to 'lll'
        metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Parent metric tensor for changing indices.
            Already assumes the value of the metric tensor from which it was initialized if passed with None. 
            Compulsory if not initialized with 'from_metric'. Defaults to None.

        Returns
        -------
        ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            New tensor with new configuration. Defaults to 'lll'

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
        new_obj = ChristoffelSymbols(
            new_tensor, self.syms, config=newconfig, parent_metric=metric
        )
        return new_obj

    def symbols(self):
        """
        Returns the symbols used for defining the time & spacial axis

        Returns
        -------
        tuple
            tuple containing (t,x1,x2,x3)
        
        """
        return self.syms
