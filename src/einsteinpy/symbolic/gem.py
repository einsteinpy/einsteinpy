import numpy as np
import sympy

from .helpers import simplify_sympy_array
from .tensor import BaseRelativityTensor, Tensor, tensor_product
from .vector import GenericVector
from .spacetime import LeviCivitaAlternatingTensor
from .optdecomposition import OPTDecompositionTensor





class GravitoElectricTensor(OPTDecompositionTensor):
    """
    Class for defining the gravitoelectric tensor as part of a 1+3 decomposition.

    Attributes
    ----------
    arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
        Raw Tensor in sympy array
    normal_vector : GenericVector
        The normal unit timelike vector used in the 1+3 decomposition
    syms : list or tuple
        List of symbols denoting space and time axis
    dims : int
        dimension of the space-time.
    variables : list
        free variables in the tensor expression other than the variables describing space-time axis.
    functions : list
        Undefined functions in the tensor expression.
    name : str or None
        Name of the tensor. Defaults to "GenericTensor".

    """

    def __init__(
        self,
        arr,
        nvec,
        syms,
        config="ll",
        parent_metric=None,
        variables=list(),
        functions=list(),
        name="GravitoElectricTensor",
    ):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        nvec : GenericVector
            The normal unit timelike vector used in the 1+3 decomposition
        syms : tuple or list
            List of crucial symbols dentoting time-axis and/or spacial axis.
            For example, in case of 4D space-time, the arrangement would look like [t, x1, x2, x3].
        config : str
            Configuration of contravariant and covariant indices in tensor.
            'u' for upper and 'l' for lower indices. Defaults to 'll'.
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Metric Tensor for some particular space-time which is associated with this Tensor.
        variables : tuple or list or set
            List of symbols used in expressing the tensor,
            other than symbols associated with denoting the space-time axis.
            Calculates in real-time if left blank.
        functions : tuple or list or set
            List of symbolic functions used in epressing the tensor.
            Calculates in real-time if left blank.
        name : str or None
            Name of the Tensor. Defaults to "GenericTensor".

        Raises
        ------
        TypeError
            Raised when arr is not a list, sympy array or numpy array.
        TypeError
            Raised when config is not of type str or contains characters other than 'l' or 'u'
        TypeError
            Raised when arguments syms, variables, functions have data type other than list, tuple or set.
        TypeError
            Raised when argument parent_metric does not belong to MetricTensor class and isn't None.
        ValueError
            Raised when argument ``syms`` does not agree with shape of argument ``arr``

        """
        super(GravitoElectricTensor, self).__init__(arr=arr, nvec=nvec, syms=syms, config=config, parent_metric=parent_metric, variables=variables, functions=functions, name=name)



    @classmethod
    def from_weyl(cls, weyl, nvec, metric=None):
        """
        Get gravitoelectric tensor from the weyl tensor and a normal timelike unit vector

        Parameters
        ----------
        weyl  :  ~einsteinpy.symbolic.weyl.WeylTensor
            WeylTensor from which to calculate the gravitoelectric tensor
        nvec  :  ~einsteinpy.symbolic.vector.GenericVector
            The normal timelike unit vector as part of the 1+3 decomposition

        """
        if metric is None:
            metric = weyl.parent_metric

        C = weyl.change_config(newconfig="llll", metric=metric)
        u = nvec.change_config("u", metric=metric)

        E = tensor_product(C, u, 1, 0)
        E = tensor_product(E, u, 2, 0)

        return cls(simplify_sympy_array(E.tensor()), nvec, weyl.syms, config="ll", parent_metric=metric)


class GravitoMagneticTensor(OPTDecompositionTensor):
    """
    Class for defining the gravitomagnetic tensor as part of a 1+3 decomposition.

    Attributes
    ----------
    arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
        Raw Tensor in sympy array
    normal_vector : GenericVector
        The normal unit timelike vector used in the 1+3 decomposition
    syms : list or tuple
        List of symbols denoting space and time axis
    dims : int
        dimension of the space-time.
    variables : list
        free variables in the tensor expression other than the variables describing space-time axis.
    functions : list
        Undefined functions in the tensor expression.
    name : str or None
        Name of the tensor. Defaults to "GenericTensor".

    """

    def __init__(
        self,
        arr,
        nvec,
        syms,
        config="ll",
        parent_metric=None,
        variables=list(),
        functions=list(),
        name="GravitoMagneticTensor",
    ):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        nvec : GenericVector
            The normal unit timelike vector used in the 1+3 decomposition
        syms : tuple or list
            List of crucial symbols dentoting time-axis and/or spacial axis.
            For example, in case of 4D space-time, the arrangement would look like [t, x1, x2, x3].
        config : str
            Configuration of contravariant and covariant indices in tensor.
            'u' for upper and 'l' for lower indices. Defaults to 'll'.
        parent_metric : ~einsteinpy.symbolic.metric.MetricTensor or None
            Metric Tensor for some particular space-time which is associated with this Tensor.
        variables : tuple or list or set
            List of symbols used in expressing the tensor,
            other than symbols associated with denoting the space-time axis.
            Calculates in real-time if left blank.
        functions : tuple or list or set
            List of symbolic functions used in epressing the tensor.
            Calculates in real-time if left blank.
        name : str or None
            Name of the Tensor. Defaults to "GenericTensor".

        Raises
        ------
        TypeError
            Raised when arr is not a list, sympy array or numpy array.
        TypeError
            Raised when config is not of type str or contains characters other than 'l' or 'u'
        TypeError
            Raised when arguments syms, variables, functions have data type other than list, tuple or set.
        TypeError
            Raised when argument parent_metric does not belong to MetricTensor class and isn't None.
        ValueError
            Raised when argument ``syms`` does not agree with shape of argument ``arr``

        """
        super(GravitoMagneticTensor, self).__init__(arr=arr,
                                                    nvec=nvec,
                                                    syms=syms,
                                                    config=config,
                                                    parent_metric=parent_metric,
                                                    variables=variables,
                                                    functions=functions,
                                                    name=name)



    @classmethod
    def from_weyl(cls, weyl, nvec, metric=None, levi_civita=None):
        """
        Get gravitomagnetic tensor from the weyl tensor and a normal timelike unit vector

        Parameters
        ----------
        weyl  :  ~einsteinpy.symbolic.weyl.WeylTensor
            WeylTensor from which to calculate the gravitoelectric tensor
        nvec  :  ~einsteinpy.symbolic.vector.GenericVector
            The normal timelike unit vector as part of the 1+3 decomposition

        """
        if metric is None:
            metric = weyl.parent_metric
        if levi_civita is None:
            levi_civita = LeviCivitaAlternatingTensor.from_metric(metric)

        C = weyl.change_config(newconfig="llll", metric=metric)

        C_s = levi_civita.GetDualTensor(C)
        u = nvec.change_config("u", metric=metric)

        H = tensor_product(C_s, u, 1, 0)
        H = tensor_product(H, u, 2, 0)

        return cls(simplify_sympy_array(H.tensor()), nvec, weyl.syms, config="ll", parent_metric=metric)
