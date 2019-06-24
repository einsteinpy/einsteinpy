import numpy as np
import sympy


class Tensor:
    """
    Base Class
    """

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

    def simplify(self):
        """
        Returns a simplified Tensor

        Returns
        -------
        ~einsteinpy.symbolic.tensor.Tensor
            Simplified Tensor

        """
        return sympy.simplify(self.tensor())
