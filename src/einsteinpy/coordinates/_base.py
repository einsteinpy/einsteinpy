import astropy.units as u


class _BaseCoordinates:
    """
    Parent Class for defining coordinate systems
    """

    def __init__(self, param1, param2, param3):

        self.param1 = param1
        self.param2 = param2
        self.param3 = param3

    def norm(self):
        raise NotImplementedError
