import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from einsteinpy.bodies import Body
from einsteinpy.geodesic import Geodesic
from einsteinpy.metric import Schwarzschild


class ScatterGeodesicPlotter:
    """
    Class for plotting static matplotlib plots.
    """

    def __init__(self, time=0 * u.s, attractor_color="black", cmap_color="Oranges"):
        """
        Constructor.

        Parameters
        ----------
        time : ~astropy.units.quantity.Quantity
            Time of start, defaults to 0 seconds.
        attractor_color : string, optional
            Color which is used to denote the attractor. Defaults to black.
        cmap_color : string, optional
            Color used in function plot.

        """
        self.time = time
        self._attractor_present = False
        self.attractor_color = attractor_color
        self.cmap_color = cmap_color

    def _plot_attractor(self):
        self._attractor_present = True
        plt.scatter(0, 0, color=self.attractor_color)

    def plot(self, geodesic):
        """

        Parameters
        ----------
        coords : ~einsteinpy.coordinates.velocity.SphericalDifferential
            Position and velocity components of particle in Spherical Coordinates.
        end_lambda : float, optional
            Lambda where iteartions will stop.
        step_size : float, optional
            Step size for the ODE.

        """
        self.mass = geodesic.attractor.mass

        vals = geodesic.trajectory

        time = vals[:, 0]
        r = vals[:, 1]
        # Currently not being used (might be useful in future)
        # theta = vals[:, 2]
        phi = vals[:, 3]

        pos_x = r * np.cos(phi)
        pos_y = r * np.sin(phi)

        plt.scatter(pos_x, pos_y, s=1, c=time, cmap=self.cmap_color)

        if not self._attractor_present:
            self._plot_attractor()

    def animate(self, geodesic, interval=50):
        """
        Function to generate animated plots of geodesics.

        Parameters
        ----------
        geodesic : ~einsteinpy.geodesic.Geodesic
            Geodesic of the body
        interval : int, optional
            Control the time between frames. Add time in milliseconds.

        """

        vals = geodesic.trajectory

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
        pic = plt.scatter([], [], s=1, c=[])
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
