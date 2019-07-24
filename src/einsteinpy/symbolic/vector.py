from einsteinpy.symbolic.tensor import Tensor


class GenericVector(Tensor):
    """
    Class to represent a vector in arbitrary space-time symbolically
    """

    def __init__(self, arr, config="l"):
        """
        Constructor and Initializer

        Parameters
        ----------
        arr : ~sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
            Sympy Array containing Sympy Expressions

        config : str
            Configuration of contravariant and covariant indices in tensor. 'u' for upper and 'l' for lower indices. Defaults to 'l'.

        Raises
        ------
        TypeError
            Dimension should be equal to 1

        """
        super(GenericVector, self).__init__(arr, config=config)
        self.arr = arr
        if self.arr.rank() == 1:
            self._order = 1
        else:
            raise TypeError("Dimension should be equal to 1")
