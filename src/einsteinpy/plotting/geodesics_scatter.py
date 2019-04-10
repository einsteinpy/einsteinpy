import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from einsteinpy.metric import Schwarzschild


class ScatterGeodesicPlotter:
    """
    Class for plotting static matplotlib plots.
    """

    def __init__(self, mass, time=0 * u.s):
        self.mass = mass
        self.time = time
        self._attractor_present = False

    def _plot_attractor(self):
        self._attractor_present = True
        plt.scatter(0, 0, color="black")

    def plot(self, pos_vec, vel_vec, end_lambda=10, step_size=1e-3):
        swc = Schwarzschild.from_spherical(pos_vec, vel_vec, self.time, self.mass)

        vals = swc.calculate_trajectory(
            end_lambda=end_lambda, OdeMethodKwargs={"stepsize": step_size}
        )[1]

        time = vals[:, 0]
        r = vals[:, 1]
        # Currently not being used (might be useful in future)
        # theta = vals[:, 2]
        phi = vals[:, 3]

        pos_x = r * np.cos(phi)
        pos_y = r * np.sin(phi)

        plt.scatter(pos_x, pos_y, s=1, c=time, cmap="Oranges")

        if not self._attractor_present:
            self._plot_attractor()

    def show(self):
        plt.show()

    def save(self, name="scatter_geodesic.png"):
        plt.savefig(name)
