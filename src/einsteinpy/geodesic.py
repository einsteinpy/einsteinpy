"""
In differential geometry, a geodesic is a curve representing in some sense the
shortest path between two points in a surface.
"""
import astropy.units as u

from einsteinpy.metric import Schwarzschild


class TimelikeGeodesic:
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
        self._trajectory = self.calculate_trajectory( # - ????
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
        coords="BL", # ????
        return_cartesian=False,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation
        """
        pass

    def calculate_trajectory_iterator(
        self,
        start_lambda=0.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        coords="BL", # ????
        return_cartesian=False,
    ):
        """
        Calculate trajectory in manifold according to geodesic equation.
        Yields an iterator.
        """
        pass
    # DRAFT CHANGES/ADDITIONS - ????

    @property
    def trajectory(self):
        return self._trajectory


# DRAFT CHANGES/ADDITIONS - ????
class NullGeodesic(Geodesic):
    """
    Class for defining Null Geodesics in different Geometries
    """

    def __init__(
        self, coords, end_lambda, start_lambda=0, step_size=1e-3, metric=Schwarzschild
    ):
        """
        Parameters
        ----------
        coords : einsteinpy.coordinates.*
            Initial Coordinates of the Photon
        end_lambda : float
            Affine Parameter value, where iterations will end
        start_lambda : float, optional
            Affine Parameter value, where iterations start
            Defaults to 0
        step_size : float, optional
            Size of each increment in proper time
            Defaults to ``1e-3``
        metric : einsteinpy.metric.*
            Geometry to work in
            Defaults to ``Schwarzschild``.
        """
        self.metric = metric.from_coords(
            coords=self.body.coordinates,
            M=self.attractor.mass,
            q=self.body.q,
            Q=self.attractor.q,
            time=time,
            a=self.attractor.a,
        )
        self._trajectory = self.calculate_trajectory( # - ????
            end_lambda=end_lambda,
            OdeMethodKwargs={"stepsize": step_size},
            return_cartesian=True,
        )[1]

    def __repr__(self):
        pass

    def __str__(self):
        pass


    def calculate_trajectory(
        self,
        start_lambda=0.0,
        end_lambda=10.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        coords="BL", # ????
        return_cartesian=False,
    ):
        """
        Calculates photon trajectory, by solving Geodesic Equation
        """
        pass

    def calculate_trajectory_iterator(
        self,
        start_lambda=0.0,
        stop_on_singularity=True,
        OdeMethodKwargs={"stepsize": 1e-3},
        coords="BL", # ????
        return_cartesian=False,
    ):
        """
        Calculates photon trajectory, by solving Geodesic Equation
        Yields an iterator
        """
        pass

    @property
    def trajectory(self):
        return self._trajectory
# DRAFT CHANGES/ADDITIONS - ????