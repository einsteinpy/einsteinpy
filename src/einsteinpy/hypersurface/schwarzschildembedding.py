import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from mpl_toolkits import mplot3d


class SchwarzschildEmbedding:
    """
    Class for Utility functions for Schwarzschild Embedding surface to
    implement gravitational lensing

    Attributes
    ----------
    input_units : list
        list of input units of M
    units_list : list
        customized units to handle values of M and render plots
        within grid range
    r_init : ~astropy.units.m

    """

    def __init__(self, M):

        """
        Constructor
        Initialize mass and embedding initial radial coordinate in appropiate units
        in order to render the plots of the surface in finite grid. The initial r
        is taken to be just greater than schwarzschild radius but it is important
        to note that the embedding breaks at r < 9m/4.

        Parameters
        ----------
        M : ~astropy.units.kg
            Mass of the body

        """
        self.input_units = [M.unit]
        self.units_list = [u.kg * 10e22, u.m / M.to(u.kg * 10e22).value]
        M = M.to(self.units_list[0])
        self.M = M
        self.r_init = (((3 * self.M.value + 0.0001) / self.M.value) * u.m).to(
            self.units_list[1]
        )

    def gradient(self, r):

        """
        Calculate gradient of Z coordinate w.r.t r to update the value of r and
        thereby get value of spherical radial coordinate R.

        Parameters
        ----------
        r : float
            schwarzschild coordinate at which gradient is supposed to be obtained

        Returns
        -------
        float
            gradient of Z w.r.t r at the point r (passed as argument)

        """
        R = r / np.sqrt(1 - (2 * self.M.value / r))
        num_one = 1 - (3 * self.M.value / r)
        num_two = np.sqrt(
            ((4 * self.M.value * r - 9 * self.M.value * self.M.value) * R)
            / (r - 3 * self.M.value) ** 2
        )
        deno = np.sqrt(1 - (2 * self.M.value / r)) ** 3

        return num_one * num_two / deno

    def radial_coord(self, r):

        """
        Returns spherical radial coordinate (of the embedding) from given schwarzschild
        coordinate.

        Parameters
        ----------
        r : float

        Returns
        -------
        float
            spherical radial coordinate of the 3d embedding

        """
        return r / np.sqrt(1 - (2 * self.M.value / r))

    def get_values(self, alpha):

        """
        Obtain the Z coordinate values and corrosponding R values for range of
        r as 9m/4 < r < 9m.

        Parameters
        ----------
        alpha : float
            scaling factor to obtain the step size for incrementing r

        Returns
        -------
        tuple
            (list, list) : values of R (x_axis) and Z (y_axis)

        """
        x_axis = []
        y_axis = []
        r_initial = self.r_init.value
        r_step = self.M.value / alpha

        z = 0
        r = r_initial
        while r < 9 * self.M.value:
            x_axis.append(self.radial_coord(r))
            y_axis.append(z)
            z = z + self.gradient(r) * r_step
            r = r + r_step

        z = 0
        r = r_initial
        while r > (9 * self.M.value / 4):
            x_axis.append(self.radial_coord(r))
            y_axis.append(z)
            z = z + self.gradient(r) * r_step
            r = r - r_step

        return x_axis, y_axis

    def get_values_surface(self, alpha):

        """
        Obtain the same values as of the get_values function but reshapes them to obtain
        values for all points on the solid of revolution about Z axis (as the
        embedding is symmetric in angular coordinates).

        Parameters
        ----------
        alpha : float
            scaling factor to obtain the step size for incrementing r
        
        Returns
        -------
        tuple
            (~numpy.array of X, ~numpy.array of Y, ~numpy.array of Z) values in cartesian coordinates
            obtained after applying solid of revolution
        
        """
        r_initial = self.r_init.value
        r_step = self.M.value / alpha
        phi_values = np.linspace(0, 2 * np.pi, 60)
        R_values = []
        z_values = []

        z = 0
        r = r_initial
        while r < 20 * self.M.value:
            R_values.append(self.radial_coord(r))
            z_values.append(z)
            z = z + self.gradient(r) * r_step
            r = r + r_step

        R_values = np.array(R_values)
        R_values, phi_values = np.meshgrid(R_values, phi_values)

        X = R_values * np.cos(phi_values)
        Y = R_values * np.sin(phi_values)
        x_len = X.shape[0]
        Z = np.array(z_values)
        z_values = np.array(z_values)
        for i in range(0, x_len - 1):
            Z = np.concatenate((Z, z_values), axis=0)

        Z.reshape((X.shape[0], X.shape[1]))
        return X, Y, Z

    def plot_hypersurface(self, plot_type="wireframe", alpha=100):

        """
        Plots the surface thus obtained for the embedding.

        Parameters
        ----------
        plot_type : str
            type of texture for the plots - wireframe / surface, defaults to 'wireframe'
        alpha : float
            scaling factor to obtain the step size for incrementing r, defaults to 100

        """
        fig = plt.figure()
        ax = plt.axes(projection="3d")
        X, Y, Z = self.get_values_surface(alpha)
        shape_tuple = X.shape
        Z = Z.reshape((shape_tuple[0], shape_tuple[1]))
        if plot_type == "wireframe":
            ax.plot_wireframe(X, Y, Z, color="black")
        elif plot_type == "surface":
            ax.plot_surface(
                X, Y, Z, rstride=1, cstride=1, cmap="cubehelix", edgecolor="none"
            )

    def show(self):
        """
        Show the plot made by plot_hypersurface()
        """
        plt.show()
