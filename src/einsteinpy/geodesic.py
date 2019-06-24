from einsteinpy.metric import Schwarzschild


class Geodesic:
    """
    Class for defining geodesics of different geometries.
    """

    def __init__(
        self, body, attractor, time, end_lambda, step_size=1e-3, metric=Schwarzschild
    ):
        """
        Parameters
        ----------

        """
        self.body = body
        self.attractor = attractor
        self.metric = metric.from_coords(
            coords=self.body.coordinates,
            M=self.attractor.mass,
            q=self.body.q,
            Q=self.attractor.q,
            time=time,
            a=self.attractor.a,
        )
        self._trajectory = self.metric.calculate_trajectory(
            end_lambda=end_lambda, OdeMethodKwargs={"stepsize": step_size}
        )[1]

    @property
    def trajectory(self):
        return self._trajectory
