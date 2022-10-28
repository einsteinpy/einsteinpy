import numpy as np
import sympy
from sympy import simplify, tensorcontraction, tensorproduct
from sympy.core.expr import Expr
from sympy.core.function import AppliedUndef, UndefinedFunction

from einsteinpy.symbolic.helpers import (
    _change_name,
    simplify_sympy_array,
    sympy_to_np_array,
)


def _config_checker(config):
    # check if the string for config contains 'u' and 'l' only
    if not isinstance(config, str):
        return False
    for ch in config:
        if (not ch == "l") and (not ch == "u"):
            return False
    return True


def _difference_list(newconfig, oldconfig):
    # defines a list of actions to be taken on a tensor
    difflist = list()
    for n_ch, o_ch in zip(newconfig, oldconfig):
        if n_ch == o_ch:
            difflist.append(0)
        elif n_ch == "u":
            difflist.append(1)
        else:
            difflist.append(-1)
    return difflist


def _change_config(tensor, metric, newconfig):
    # check length and validity of new configuration
    if not (len(newconfig) == len(tensor.config) and _config_checker(newconfig)):
        raise ValueError

    # seperate the contravariant & covariant metric tensors
    met_dict = {
        -1: metric.lower_config().tensor(),
        1: metric.lower_config().inv().tensor(),
    }

    # main code
    def chain_config_change():
        t = sympy.Array(tensor.tensor())
        difflist = _difference_list(newconfig, tensor.config)
        for i, action in enumerate(difflist):
            if action == 0:
                continue
            else:
                t = simplify(
                    tensorcontraction(tensorproduct(met_dict[action], t), (1, 2 + i))
                )
                # reshuffle the indices
                dest = list(range(len(t.shape)))
                dest.remove(0)
                dest.insert(i, 0)
                t = sympy.permutedims(t, dest)
        return t

    return chain_config_change()


def tensor_product(tensor1, tensor2, i=None, j=None):
    """Tensor Product of ``tensor1`` and ``tensor2``

    Parameters
    ----------
    tensor1 : ~einsteinpy.symbolic.BaseRelativityTensor
    tensor2 : ~einsteinpy.symbolic.BaseRelativityTensor
    i : int, optional
        contract ``i``th index of ``tensor1``
    j : int, optional
        contract ``j``th index of ``tensor2``


    Returns
    -------
    ~einsteinpy.symbolic.BaseRelativityTensor
        tensor of appropriate rank

    Raises
    ------
    ValueError
        Raised when ``i`` and ``j`` both indicate 'u' or 'l' indices
    """
    product = tensorproduct(tensor1.arr, tensor2.arr)

    if (i or j) is None:
        newconfig = tensor1.config + tensor2.config
    else:
        if tensor1.config[i] == tensor2.config[j]:
            raise ValueError(
                "Index summation not allowed between %s and %s indices"
                % (tensor1.config[i], tensor2.config[j])
            )

        product = simplify(tensorcontraction(product, (i, len(tensor1.config) + j)))

        con = tensor1.config[:i] + tensor1.config[i + 1 :]
        fig = tensor2.config[:j] + tensor2.config[j + 1 :]
        newconfig = con + fig

    return BaseRelativityTensor(
        product,
        syms=tensor1.syms,
        config=newconfig,
        parent_metric=tensor1.parent_metric,
        variables=tensor1.variables,
        functions=tensor1.functions,
    )


class Tensor:
    """
    Base Class for Tensor manipulation
    """

    def __init__(self, arr, config="ll", name=None):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array, multi-dimensional list containing Sympy Expressions, or Sympy Expressions or int or float scalar
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'll'.
        name : str or None
            Name of the tensor.

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy array
        TypeError
            Raised when config is not of type str or contains characters other than 'l' or 'u'
        ValueError
            Raised when ``config`` implies order of Tensor different than that indicated by shape of ``arr``

        """

        if isinstance(arr, (list, tuple, np.ndarray, int, float, np.number, Expr)):
            self.arr = sympy.Array(arr)
        elif isinstance(arr, sympy.Array):
            self.arr = arr
        else:
            raise TypeError("Only multi-dimensional list or Sympy Array is expected")
        if _config_checker(config):
            self._config = config
            self._order = len(config)
        else:
            raise TypeError(
                "config is either not of type 'str' or does contain characters other than 'l' or 'u'"
            )
        if len(self.arr.shape) != len(config):
            raise ValueError(
                "invalid shape of array for tensor of order implied by config: '{}'".format(
                    config
                )
            )
        self.name = name

    @property
    def order(self):
        """
        Returns the order of the Tensor

        """
        return self._order

    @property
    def config(self):
        """
        Returns the configuration of covariant and contravariant indices

        """
        return self._config

    def __getitem__(self, index):
        return self.arr[index]

    def __str__(self):
        """
        Returns a String with a readable representation of the object of class Tensor

        """
        representation = "Tensor"
        if self.name is not None:
            representation = " ".join((representation, self.name))
        representation += "\n"
        representation += self.arr.__str__()
        return representation

    def __repr__(self):
        """
        Returns a String with a representation of the state of the object of class Tensor

        """
        interpretable_representation = self.__class__.__name__
        interpretable_representation += self.arr.__repr__()
        return interpretable_representation

    def tensor(self):
        """
        Returns the sympy Array

        Returns
        -------
        ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
            Sympy Array object

        """
        return self.arr

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
        ~einsteinpy.symbolic.tensor.Tensor:
            Tensor with substituted values

        """
        return Tensor(self.tensor().subs(*args))

    def simplify(self, set_self=True):
        """
        Returns a simplified Tensor

        Parameters
        ----------
        set_self : bool
            Replaces the tensor contained the class with its simplified version, if ``True``.
            Defaults to ``True``.

        Returns
        -------
        ~einsteinpy.symbolic.tensor.Tensor
            Simplified Tensor

        """
        if set_self:
            self.arr = simplify_sympy_array(self.tensor())
            return self.tensor()
        # return sympy.simplify(self.tensor())  # this used to work with older sympy versions
        return simplify_sympy_array(self.tensor())


class BaseRelativityTensor(Tensor):
    """
    Generic class for defining tensors in General Relativity.
    This would act as a base class for other Tensorial quantities in GR.

    Attributes
    ----------
    arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
        Raw Tensor in sympy array
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
        syms,
        config="ll",
        parent_metric=None,
        variables=list(),
        functions=list(),
        name="GenericTensor",
    ):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
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
        super(BaseRelativityTensor, self).__init__(arr=arr, config=config, name=name)

        if len(self.arr.shape) != 0 and self.arr.shape[0] != len(syms):
            raise ValueError("invalid shape of argument arr for syms: {}".format(syms))

        # Cannot implement the check that parent metric belongs to the class MetricTensor
        # Due to the issue of cyclic imports, would find a workaround
        self._parent_metric = parent_metric
        if isinstance(syms, (list, tuple)):
            self.syms = syms
            self.dims = len(self.syms)
        else:
            raise TypeError("syms should be a list or tuple")

        if isinstance(variables, (list, tuple, set)) and isinstance(
            functions, (list, tuple, set)
        ):
            # compute free variables and functions if list if empty
            if not variables:
                self.variables = [
                    v for v in self.arr.free_symbols if v not in self.syms
                ]
                self.variables.sort(key=(lambda var: var.name))
            else:
                self.variables = list(variables)
            if not functions:
                self.functions = [
                    f
                    for f in self.arr.atoms(AppliedUndef).union(
                        self.arr.atoms(UndefinedFunction)
                    )
                ]
            else:
                self.functions = list(functions)

        else:
            raise TypeError(
                "arguments variables and functions should be a list, tuple or set"
            )

    @property
    def parent_metric(self):
        """
        Returns the Metric from which Tensor was derived/associated, if available.
        """
        return self._parent_metric

    def symbols(self):
        """
        Returns the symbols used for defining the time & spacial axis

        Returns
        -------
        tuple
            tuple containing (t,x1,x2,x3) in case of 4D space-time

        """
        return self.syms

    def tensor_lambdify(self, *args):
        """
        Returns lambdified function of symbolic tensors.
        This means that the returned functions can accept numerical values and return numerical quantities.

        Parameters
        ----------
            *args
                The variable number of arguments accept sympy symbols.
                The returned function accepts arguments in same order as initially defined in ``*args``.
                Uses sympy symbols from class attributes ``syms`` and ``variables`` (in the same order) if no ``*args`` is passed
                Leaving ``*args`` empty is recommended.

        Returns
        -------
            tuple
                arguments to be passed in the returned function in exact order.
            function
                Lambdified function which accepts and returns numerical quantities.

        """

        if len(args) == 0:
            numeric_arr = sympy.lambdify(
                [*self.syms, *self.variables], self.arr, modules="numpy"
            )
            arg_list = (*self.syms, *self.variables)
        else:
            numeric_arr = sympy.lambdify(args, self.arr, modules="numpy")
            arg_list = tuple(args)
        return arg_list, numeric_arr

    def lorentz_transform(self, transformation_matrix):
        """
        Performs a Lorentz transform on the tensor.

        Parameters
        ----------
            transformation_matrix : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
                Sympy Array or multi-dimensional list containing Sympy Expressions

        Returns
        -------
            ~einsteinpy.symbolic.tensor.BaseRelativityTensor
                lorentz transformed tensor(or vector)

        """
        tm = sympy.Array(transformation_matrix)
        t = self.tensor()
        for i in range(self.order):
            if self.config[i] == "u":
                t = simplify(tensorcontraction(tensorproduct(tm, t), (1, 2 + i)))
            else:
                t = simplify(tensorcontraction(tensorproduct(tm, t), (0, 2 + i)))
            dest = list(range(len(t.shape)))
            dest.remove(0)
            dest.insert(i, 0)
            t = sympy.permutedims(t, dest)

        return BaseRelativityTensor(
            t,
            syms=self.syms,
            config=self.config,
            parent_metric=None,
            variables=self.variables,
            functions=self.functions,
            name=_change_name(self.name, context="__lt"),
        )
