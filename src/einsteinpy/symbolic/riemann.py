import numpy as np
import sympy

from einsteinpy.symbolic.christoffel import ChristoffelSymbols
from einsteinpy.symbolic.tensor import Tensor, _change_config


class RiemannCurvatureTensor(Tensor):
    """
    Class for defining Riemann Curvature Tensor
    """

    def __init__(self, arr, syms, config="ulll", parent_metric=None):
        """
        Constructor and Initializer
        
        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        syms : tuple or list
            Tuple of crucial symbols dentoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'ulll'.

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 4 indices
        
        """
        super(RiemannCurvatureTensor, self).__init__(arr, config=config)
        self._order = 4
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
        Returns the Parent Metric, if available.
        """
        return self._parent_metric

    @classmethod
    def from_christoffels(cls, chris, parent_metric=None):
        """
        Get Riemann Tensor calculated from a Christoffel Symbols

        Parameters
        ----------
        chris : ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            Christoffel Symbols from which Riemann Curvature Tensor to be calculated
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corrosponding Metric for the Riemann Tensor. 
            None if it should inherit the Parent Metric of Christoffel Symbols.
            Defaults to None.
        
        """
        if not chris.config == "ull":
            chris = chris.change_config(newconfig="ull", metric=parent_metric)
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
        if parent_metric is None:
            parent_metric = chris.parent_metric
        return cls(riemann_list, syms, config="ulll", parent_metric=parent_metric)

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
        return cls.from_christoffels(ch, parent_metric=None)

    def change_config(self, newconfig="llll", metric=None):
        """
        Changes the index configuration(contravariant/covariant)

        Parameters
        ----------
        newconfig : str
            Specify the new configuration. Defaults to 'llll'
        metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Parent metric tensor for changing indices.
            Already assumes the value of the metric tensor from which it was initialized if passed with None. 
            Compulsory if not initialized with 'from_metric'. Defaults to None.

        Returns
        -------
        ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            New tensor with new configuration. Defaults to 'llll'

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
        new_obj = RiemannCurvatureTensor(
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
