"""
In differential geometry, a geodesic is a curve representing in some sense the
shortest path between two points in a surface.
"""
import astropy.units as u

from einsteinpy.metric import Schwarzschild


class Geodesic:
    """
    Class for defining Timelike Geodesics in different Geometries

    """

    def __init__(
        self, body, end_lambda, step_size=1e-3, time=0 * u.s, metric=Schwarzschild
    ):
        """
        Parameters
        ----------
        body : ~einsteinpy.bodies.Body
            Test particle for which Timelike Geodesic is to be calculated.
        end_lambda : float
            Lambda(proper time in seconds) where iterations will stop
        step_size : float, optional
            Size of each increment in proper time.
            Defaults to ``1e-3``.
        time : ~astropy.units.s, optional
            Time of start, Defaults to 0 seconds.
        metric : ~einsteinpy.metric.schwarzschild.Schwarzschild or ~einsteinpy.metric.kerr.Kerr or ~einsteinpy.metric.kerrnewman.KerrNewman, optional
            Class of the spacetime metric in which geodesics are to be calculated.
            Defaults to ``Schwarzschild``.
        """
        self.body = body
        self.attractor = body.parent
        self.metric = metric.from_coords(
            coords=self.body.coordinates,
            M=self.attractor.mass,
            q=self.body.q,
            Q=self.attractor.q,
            time=time,
            a=self.attractor.a,
        )
        self._trajectory = self.calculate_trajectory(  # - ?????
            end_lambda=end_lambda,
            OdeMethodKwargs={"stepsize": step_size},
            return_cartesian=True,
        )[1]

    def __repr__(self):
        return "Geodesic object body name= ({0}) , metric=({1}) , parent name=({2}) , parent mass=({3})".format(
            self.body.name,
            self.metric.name,
            self.body.parent.name,
            self.body.parent.mass,
        )

    def __str__(self):
        return "Geodesic object body name= ({0}) , metric=({1}) , parent name=({2}) , parent mass=({3})".format(
            self.body.name,
            self.metric.name,
            self.body.parent.name,
            self.body.parent.mass,
        )

    # DRAFT CHANGES/ADDITIONS - ????
    def calculate_trajectory(
        self,
        start_lambda=0.0,
        end_lambda=10.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        coords="BL",  # ????
        return_cartesian=False,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation
        """
        raise NotImplementedError

    def calculate_trajectory_iterator(
        self,
        start_lambda=0.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        coords="BL",  # ????
        return_cartesian=False,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation.
        Yields an iterator.
        """
        raise NotImplementedError

    # DRAFT CHANGES/ADDITIONS - ????

    @property
    def trajectory(self):
        return self._trajectory


# DRAFT CHANGES/ADDITIONS - ????
class NullGeodesic:
    """
    Class for defining Null Geodesics in different Geometries
    """

    def __init__(self, init_coords, end_lambda, step_size=1e-3, metric=Schwarzschild):
        """
        Parameters
        ----------
        init_coords : einsteinpy.coordinates.* - ?????
            Initial 4-Position and 4-Momentum of the Photon
        end_lambda : float
            Affine Parameter value, where iterations will end
        step_size : float, optional
            Size of each increment in proper time
            Defaults to ``1e-3``
        metric : einsteinpy.metric.* - ?????
            Spacetime Geometry
            Defaults to ``Schwarzschild``.
        """
        raise NotImplementedError

    def __repr__(self):
        raise NotImplementedError

    def __str__(self):
        raise NotImplementedError

    def calculate_trajectory(
        self, end_lambda=10.0, OdeMethodKwargs={"stepsize": 1e-3}, coords="BL",  # ????
    ):
        """
        Calculates photon trajectory, by solving Geodesic Equation

        Allow users to convert coordinates after calculations are done
        So as to not slow down the function exit - ?????
        """
        raise NotImplementedError

    def calculate_trajectory_iterator(
        self, OdeMethodKwargs={"stepsize": 1e-3}, coords="BL",  # ????
    ):
        """
        Calculates photon trajectory, by solving Geodesic Equation
        Yields an iterator
        """
        raise NotImplementedError

    # @property
    # def trajectory(self):
    #     return self._trajectory


# DRAFT CHANGES/ADDITIONS - ????
