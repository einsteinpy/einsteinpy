import numpy as np
import sympy
from sympy import simplify, tensorcontraction, tensorproduct


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
    if metric.config == "ll":
        met_dict = {-1: metric.tensor(), 1: metric.inv().tensor()}
    else:
        met_dict = {-1: metric.inv().tensor(), 1: metric.tensor()}

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
                tmp = np.array(t).reshape(t.shape)
                source, dest = (
                    [p for p in range(len(t.shape))],
                    [p for p in range(len(t.shape))],
                )
                dest.pop(i)
                dest.insert(0, i)
                tmp = np.moveaxis(tmp, source, dest)
                t = sympy.Array(tmp)
        return t

    return chain_config_change()


class Tensor:
    """
    Base Class for Tensor manipulation
    """

    def __init__(self, arr, config="ll"):
        """
        Constructor and Initializer
        
        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions
        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'll'.

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy array
        TypeError
            Raised when config is not of type str or contains characters other than 'l' or 'u'
        

        """
        if isinstance(arr, (list, tuple)):
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

    def simplify(self):
        """
        Returns a simplified Tensor

        Returns
        -------
        ~einsteinpy.symbolic.tensor.Tensor
            Simplified Tensor

        """
        return sympy.simplify(self.tensor())
