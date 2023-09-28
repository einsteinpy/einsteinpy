import numpy as np
import sympy
from sympy.functions.special.tensor_functions import LeviCivita
#from sympy import simplify, tensorcontraction, tensorproduct

from .helpers import _change_name
from .tensor import BaseRelativityTensor, _change_config, tensor_product


class LeviCivitaAlternatingTensor(BaseRelativityTensor):

    def __init__(self, arr, syms, config="ll", parent_metric=None, parent_spacetime=None, name="LeviCivitaAlternatingTensor"):
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
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Corresponding Metric for the Ricci Tensor.
            Defaults to None.
        parent_spacetime : ~einsteinpy.symbolic.spacetime.GenericSpacetime or None
            Spacetime object associated with this Tensor.
        name : str
            Name of the Tensor. Defaults to "RicciTensor".

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 2 indices

        """
        super(LeviCivitaAlternatingTensor, self).__init__(
            arr=arr, syms=syms, config=config, parent_metric=parent_metric, parent_spacetime=parent_spacetime, name=name,
        )
        self._order = 4
        if not len(config) == self._order:
            raise ValueError("config should be of length {}".format(self._order))


    @classmethod
    def from_metric(cls, metric):
        eps = sympy.MutableDenseNDimArray(np.zeros((4,)*4))
        for mu in range(4):
            for nu in range(4):
                for alpha in range(4):
                    for beta in range(4):
                        eps[mu,nu,alpha,beta] = LeviCivita(mu, nu, alpha, beta)
        eps = sympy.sqrt(-metric.determinant()) * eps
        return cls(sympy.Array(eps), syms=metric.syms, config="llll", parent_metric=metric)


    def GetDualTensor(self, T):
        """
        Calculates the dual tensor of a rank 4 tensor

        Parameters
        ----
        T: ~einsteinpy.symbolic.tensor.BaseRelativityTensor
            The tensor of rank 4

        Returns
        ----
        ~einsteinpy.symbolic.BaseRelativityTensor
            The dual tensor

        """
        T = T.change_config("llll")
        eps = self.change_config( "uull")

        dual = 1./2. * sympy.tensorproduct(T.arr, eps.arr)
        dual = sympy.tensorcontraction(sympy.tensorcontraction(dual, (2, 4)), (2, 4))
        return BaseRelativityTensor(dual, syms=T.syms, config="llll", parent_metric=self.parent_metric)
