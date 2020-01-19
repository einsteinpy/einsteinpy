import astropy.units as u

from einsteinpy.metric import Schwarzschild


class Geodesic:
    """
    Class for defining geodesics of different geometries.
    """

    def __init__(
        self, body, end_lambda, step_size=1e-3, time=0 * u.s, metric=Schwarzschild
    ):
        """
        Parameters
        ----------

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
        self._trajectory = self.metric.calculate_trajectory(
            end_lambda=end_lambda,
            OdeMethodKwargs={"stepsize": step_size},
            return_cartesian=True,
        )[1]

    def __repr__(self):
        return "body name= ({0}) , metric=({1}) , parent name=({2}) , parent mass=({3})".format(
            self.body.name,
            self.metric.name,
            self.body.parent.name,
            self.body.parent.mass,
        )

    def __str__(self):
        return "body name= ({0}) , metric=({1}) , parent name=({2}) , parent mass=({3})".format(
            self.body.name,
            self.metric.name,
            self.body.parent.name,
            self.body.parent.mass,
        )

    @property
    def trajectory(self):
        return self._trajectory
