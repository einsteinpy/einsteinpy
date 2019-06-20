import sympy
import tensor_symbolic as Tensor


class ChristoffelSymbols(Tensor):
    def __init__(self, list2d, syms):
        """
        Constructor and Initializer
        :param list2d:
        :param syms:
        """
        super.__init__(syms)
        self.list2d = list2d
        dims = self.dims
        self.mat = sympy.Matrix(self.list2d)
        self.mat_inv = self.mat.inv()
        self.christlist = self.create_christlist(dims, self.generic_list)

    def create_christlist(self, dims, generic_list):

        """
         Method to calculate christoffel symbols of a given metric

         Parameters
         ----------
         list2d : list
             2d list (Matrix) representing metric, containing ~sympy expressions
         syms : list
             1d list containing representaion of [x0,x1,x2...] in ~sympy expressions

         Returns
         -------
         list
             3d list of ~sympy expressions containing christoffel symbols
         """

        _counterlist = [i for i in range(dims ** 3)]
        for t in _counterlist:
            # i,j,k each goes from 0 to (dims-1)
            # hack for codeclimate. Could be done with 3 nested for loops
            k = t % dims
            j = (int(t / dims)) % (dims)
            i = (int(t / (dims ** 2))) % (dims)
            temp = 0
            for n in range(dims):
                temp += (self.mat_inv[i, n] / 2) * (
                    sympy.diff(self.list2d[n][j], self.syms[k])
                    + sympy.diff(self.list2d[n][k], self.syms[j])
                    - sympy.diff(self.list2d[j][k], self.syms[n])
                )
            generic_list[i][j][k] = temp
            return generic_list

    def schwarzschild_christoffels(self, symbolstr="t r theta phi"):
        """
        Returns the 3d list of christoffel symbols of Schwarzschild Metric.

        Parameters
        ----------
        symbolstr : string
            symbols to be used to define schwarzschild space, defaults to 't r theta phi'

        Returns
        -------
        list
            3d list of christoffel symbols for schwarzschild metric

        """
        list2d = [[0 for i in range(4)] for i in range(4)]
        syms = sympy.symbols(symbolstr)
        c, a = sympy.symbols("c a")
        list2d[0][0] = 1 - (a / syms[1])
        list2d[1][1] = -1 / ((1 - (a / syms[1])) * (c ** 2))
        list2d[2][2] = -1 * (syms[1] ** 2) / (c ** 2)
        list2d[3][3] = -1 * (syms[1] ** 2) * (sympy.sin(syms[2]) ** 2) / (c ** 2)
        christoffels = ChristoffelSymbols(list2d, syms)
        return christoffels

    def kerr_christoffels(symbolstr="t r theta phi"):
        """
        Returns the 3d list of christoffel symbols of Kerr metric(BL coordinates) in Plank units : G=1, c=1.

        Parameters
        ----------
        symbolstr : string
            symbols to be used to define kerr space in BL coordinates, defaults to 't r theta phi'

        Returns
        -------
        list
            3d list of christoffel symbols for kerr metric

        """
        list2d = [[0 for i in range(4)] for i in range(4)]
        syms = sympy.symbols(symbolstr)
        a, R = sympy.symbols("a R")
        A = syms[1] ** 2 - R * syms[1] + a ** 2
        sigma = syms[1] ** 2 + (a ** 2) * (sympy.cos(syms[2]) ** 2)
        list2d[0][0] = (R * syms[1] / sigma) - 1
        list2d[1][1] = sigma / A
        list2d[2][2] = sigma
        list2d[3][3] = (
            (sympy.sin(syms[2]) ** 2)
            * (
                (a ** 2 + syms[1] ** 2) ** 2
                - (a ** 2) * (A * (sympy.sin(syms[2]) ** 2))
            )
        ) / sigma
        list2d[3][0] = -1 * (R * a * (syms[1])) * (sympy.sin(syms[2]) ** 2) / sigma
        list2d[0][3] = list2d[3][0]
        christoffels = ChristoffelSymbols(list2d, syms)
        return christoffels
