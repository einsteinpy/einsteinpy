from sympy.core.symbol import Symbol


class SymbolicConstant(Symbol):
    """
    This class inherits from ~sympy.core.symbol.Symbol 

    Parameters
    ----------
    name : str
        Short, commonly accepted name of the constant. 
        For example, 'c' for Speed of light.
    descriptive_name : str
        The extended name of the constant. 
        For example, 'Speed of Light' for 'c'. 
        Defaults to None.

    """

    def __new__(cls, name, descriptive_name=None, **assumptions):
        instance = super(SymbolicConstant, cls).__new__(cls, name, **assumptions)
        instance._descriptive_name = descriptive_name
        return instance

    @property
    def descriptive_name(self):
        """
        Returns the extended name of the constant

        """
        return self._descriptive_name


c = SymbolicConstant("c", "Speed Of Light")
G = SymbolicConstant("G", "Gravitational Constant")
Cosmo_Const = SymbolicConstant("Lambda", "Cosmological Constant")


def get_constant(name):
    """
    Returns a symbolic instance of the constant

    Parameters
    ----------
    name : str
        Name of the constant. 
        Currently available names are 'c', 'G', 'Cosmo_Const'.
    
    Returns
    -------
    ~einsteinpy.symbolic.constants.SymbolicConstant
        An instance of the required constant

    """
    const_dict = {"c": c, "G": G, "Cosmo_Const": Cosmo_Const}
    return const_dict[name]
