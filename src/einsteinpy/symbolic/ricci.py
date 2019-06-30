import numpy as np
import sympy

from .christoffel import ChristoffelSymbols
from .metric import MetricTensor
from .riemann import RiemannCurvatureTensor
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
    def from_riemann(cls, riemann):
        """
        Get Ricci Tensor calculated from Riemann Tensor

        Parameters
        ----------
        riemann : ~einsteinpy.symbolic.riemann.RiemannCurvatureTensor
           Riemann Tensor

        """
        return cls(sympy.tensorcontraction(riemann.tensor(), (0, 2)), riemann.syms)

    @classmethod
    def from_christoffels(cls, chris):
        """
        Get Ricci Tensor calculated from Christoffel Tensor

        Parameters
        ----------
        chris : ~einsteinpy.symbolic.christoffel.ChristoffelSymbols
            Christoffel Tensor

        """
        rt = RiemannCurvatureTensor.from_christoffels(chris)
        return cls.from_riemann(rt)

    @classmethod
    def from_metric(cls, metric):
        """
        Get Ricci Tensor calculated from Metric Tensor

        Parameters
        ----------
        metric : ~einsteinpy.symbolic.metric.MetricTensor
            Metric Tensor

        """
        rc = ChristoffelSymbols.from_metric(metric)
        return cls.from_christoffels(rc)

    def symbols(self):
        """
        Returns the symbols used for defining the time & spacial axis

        Returns
        -------
        tuple
            tuple containing (x1,x2)

        """
        return self.syms
