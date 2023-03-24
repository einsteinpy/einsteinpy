import numpy as np
import sympy
from sympy import simplify, tensorcontraction, tensorproduct

from .helpers import expand_sympy_array
from .tensor import BaseRelativityTensor, Tensor, tensor_product, _change_name
from .vector import GenericVector


class OPTDecompositionTensor(BaseRelativityTensor):
    """
    Generic class for defining 1+3 decomposition tensors in General Relativity.

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
    simplify : Bool
            Whether to call the simplify routine on initiation

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
        name="GenericOPTDecompositionTensor",
        simplify=True,
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
        super(OPTDecompositionTensor, self).__init__(arr=arr, syms=syms, config=config, parent_metric=parent_metric, variables=variables, functions=functions, name=name, simplify=simplify)

        # Make sure we have a unit vector ?
        self.normal_vector = nvec


    def change_config(self, config, metric=None):
        cls = self.__class__
        self.__class__ = BaseRelativityTensor
        t = self.change_config(config, metric)
        self.__class__ = cls
        return self.__class__(t.arr, nvec=self.normal_vector, syms=t.syms, config=t.config, parent_metric=t.parent_metric)

    def symmetric_part(self, i=0, j=1):
        """

        """
        t = super(OPTDecompositionTensor, self).symmetric_part(i,j)
        return OPTDecompositionTensor(t.arr, nvec=self.normal_vector, syms=t.syms, config=t.config, parent_metric=t.parent_metric)

    def antisymmetric_part(self, i=0, j=1):
        """

        """
        t = super(OPTDecompositionTensor, self).antisymmetric_part(i,j)
        return OPTDecompositionTensor(t.arr, nvec=self.normal_vector, syms=t.syms, config=t.config, parent_metric=t.parent_metric)


    def subs(self, *args):
        """
        Substitute the variables/expressions in a Tensor with other sympy variables/expressions.

        Parameters
        ----------
        args : one argument or two argument
            - two arguments, e.g foo.subs(old, new)
            - one iterable argument, e.g foo.subs([(old1, new1), (old2, new2)]) for multiple substitutions at once.

        Returns
        -------
        ~Instance of self:
            Tensor with substituted values

        """
        return self.__class__(expand_sympy_array(self.tensor()).subs(*args), nvec=self.normal_vector, syms=self.syms, config=self.config, parent_metric=self._parent_metric, name=self.name)


class OPTMetric(OPTDecompositionTensor):
    """
    Class to define a metric in a 1+3 decomposition
    """

    def __init__(self, arr, nvec, syms, config="ll", name="GenericMetricTensor"):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        nvec : GenericVector
            The normal unit timelike vector used in the 1+3 decomposition
        syms : tuple or list
            Tuple of crucial symbols denoting time-axis, 1st, 2nd, and 3rd axis (t,x1,x2,x3)
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'll'.
        name : str
            Name of the Metric Tensor. Defaults to "GenericMetricTensor".

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy Array
        TypeError
            syms is not a list or tuple
        ValueError
            config has more or less than 2 indices

        """
        super(OPTMetric, self).__init__(
            arr=arr, nvec=nvec, syms=syms, config=config, parent_metric=self, name=name
        )
        self.normal_vector._parent_metric = self
        self._order = 2
        self._invmetric = None
        self._proj_tensor = None
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
            return self.inv()

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
            newconfig = "ll" if self.config == "uu" else "uu"
            inv_met = OPTMetric(
                sympy.simplify(sympy.Matrix(self.arr.tolist()).inv()).tolist(),
                nvec=self.normal_vector,
                syms=self.syms,
                config=newconfig,
                name=_change_name(self.name, context="__" + newconfig),
            )
            inv_met._invmetric = self
            self._invmetric =  inv_met
        return self._invmetric

    def lower_config(self):
        """
        Returns a covariant instance of the given metric tensor.

        Returns
        -------
        ~einsteinpy.symbolic.metric.MetricTensor
            same instance if the configuration is already lower or
            inverse of given metric if configuration is upper

        """
        if self.config == "ll":
            return self
        return self.inv()


    @property
    def projector_tensor(self):
        """

        """
        if self._proj_tensor is None:
            u = self.normal_vector.change_config("l")
            if self.config == "uu":
                g = self.change_config("ll").tensor()
            else:
                g = self.tensor()
            uu = tensorproduct(u.tensor(), u.tensor())
            self._proj_tensor = OPTDecompositionTensor(g + uu, nvec=self.normal_vector, syms=self.syms, config="ll", parent_metric=self)
        return self._proj_tensor

    @projector_tensor.setter
    def projector_tensor(self, value):
        self._proj_tensor = value

    def determinant(self):
        """
        Returns the determinant of the given metric tensor

        Returns
        ---
            ~sympy.core.mul.Mul
                Sympy multiplication object
        """
        return sympy.Matrix(self.lower_config().tensor()).det()




