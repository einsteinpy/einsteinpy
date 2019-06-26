import numpy as np
import sympy


class Tensor:
    """
    Base Class
    """

    def __getitem__(self, index):
        return self.arr[index]

    def __setitem__(self, index, value):
        self.arr[index] = value

    def __init__(self, arr):
        """
        Constructor and Initializer
        
        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray or list
            Sympy Array or multi-dimensional list containing Sympy Expressions

        Raises
        ------
        TypeError
            Raised when arr is not a list or sympy array
        

        """
        if isinstance(arr, (list, tuple)):
            self.arr = sympy.Array(arr)
        elif isinstance(arr, sympy.Array):
            self.arr = arr
        else:
            raise TypeError("Only multi-dimensional list or Sympy Array is expected")

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
