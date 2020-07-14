from matplotlib import pyplot as plt


class HypersurfacePlotter:
    """
    Class for plotting and visualising hypersurfaces
    """

    def __init__(self, embedding, plot_type="wireframe", alpha=100):
        """
        Constructor for plotter.

        Parameters
        ----------
        embedding : ~einsteinpy.hypersurface.SchwarzschildEmbedding
            The embedding of the hypersurface.
        plot_type : str, optional
            type of texture for the plots - wireframe / surface, defaults to 'wireframe'
        alpha : float, optional
            scaling factor to obtain the step size for incrementing r, defaults to 100
        """
        self.embedding = embedding
        self.plot_type = plot_type
        self.alpha = alpha

    def plot(self):
        """
        Plots the surface thus obtained for the embedding.
        """
        fig = plt.figure()
        ax = plt.axes(projection="3d")
        X, Y, Z = self.embedding.get_values_surface(self.alpha)
        shape_tuple = X.shape
        Z = Z.reshape((shape_tuple[0], shape_tuple[1]))
        if self.plot_type == "wireframe":
            ax.plot_wireframe(X, Y, Z, color="black")
        elif self.plot_type == "surface":
            ax.plot_surface(
                X, Y, Z, rstride=1, cstride=1, cmap="cubehelix", edgecolor="none"
            )

    def show(self):
        """
        Shows the plot.
        """
        plt.show()
