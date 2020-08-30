from einsteinpy.geodesic import Geodesic


class Timelike(Geodesic):
    """
    Class for defining Time-like Geodesics

    """

    def __init__(
        self,
        position,
        momentum,
        a=0.0,
        end_lambda=50.0,
        step_size=0.0005,
        return_cartesian=True,
        julia=True,
    ):
        """
        Constructor

        Parameters
        ----------
        position : array_like
            Length-3 Array, containing the initial 3-Position
        momentum : array_like
            Length-3 Array, containing the initial 3-Momentum
        a : float, optional
            Dimensionless Spin Parameter of the Black Hole
            ``0 <= a <= 1``
            Defaults to ``0.`` (Schwarzschild Black Hole)
        end_lambda : float, optional
            Affine Parameter value, where integration will end
            Equivalent to Proper Time for Timelike Geodesics
            Defaults to ``50.``
        step_size : float, optional
            Size of each geodesic integration step
            A fixed-step, symplectic VerletLeapfrog integrator is used
            Defaults to ``0.0005``
        return_cartesian : bool, optional
            Whether to return calculated positions in Cartesian Coordinates
            This only affects the coordinates. The momenta dimensionless quantities,
            and are returned in Spherical Polar Coordinates.
            Defaults to ``True``
        julia : bool, optional
            Whether to use the julia backend
            Defaults to ``True``

        """
        super().__init__(
            position=position,
            momentum=momentum,
            a=a,
            end_lambda=end_lambda,
            step_size=step_size,
            time_like=True,
            return_cartesian=return_cartesian,
            julia=julia,
        )
