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
        """

        Parameters
        ----------
        pos_vec : list
            list of r, theta & phi components along with ~astropy.units.
        vel_vec : list
            list of velocities of r, theta & phi components along with ~astropy.units.
        end_lambda : float, optional
            Lambda where iteartions will stop.
        step_size : float, optional
            Step size for the ODE.

        """

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

    def animate(
        self, pos_vec, vel_vec, end_lambda=10, step_size=1e-3, interval=50
    ):
        """
        Function to generate animated plots of geodesics.

        Parameters
        ----------
        pos_vec : list
            list of r, theta & phi components along with ~astropy.units.
        vel_vec : list
            list of velocities of r, theta & phi components along with ~astropy.units.
        end_lambda : float, optional
            Lambda where iteartions will stop.
        step_size : float, optional
            Step size for the ODE.
        interval : int, optional
            Control the time between frames. Add time in milliseconds.

        """

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
        frames = pos_x.shape[0]
        x_max, x_min = max(pos_x), min(pos_x)
        y_max, y_min = max(pos_y), min(pos_y)
        margin_x = (x_max - x_min) * 0.1
        margin_y = (y_max - y_min) * 0.1

        fig = plt.figure()

        plt.xlim(x_min - margin_x, x_max + margin_x)
        plt.ylim(y_min - margin_y, y_max + margin_y)
        pic = plt.scatter([], [], s=1)
        plt.scatter(0, 0, color="black")

        def _update(frame):
            pic.set_offsets(np.vstack((pos_x[: frame + 1], pos_y[: frame + 1])).T)
            pic.set_array(time[: frame + 1])
            return (pic,)

        self.animated = FuncAnimation(fig, _update, frames=frames, interval=interval)

    def show(self):
        plt.show()

    def save(self, name="scatter_geodesic.png"):
        plt.savefig(name)
