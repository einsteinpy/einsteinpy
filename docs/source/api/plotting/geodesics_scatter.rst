Static 2D Plotting module
=========================

This module contains the basic classes for static plottings in
2-dimensions for scatter and line:

.. automodule:: einsteinpy.plotting.geodesics_scatter
    :members:

Color
=====
* Attractor : 
    User can give the color to attractor of his/her choice.
    It can be passed while making the object of geodesics_static class.
    Default color of attractor is "black".

    .. code-block:: python
        
        self.attractor_color = attractor_color
        plt.scatter(0, 0, color=self.attractor_color)

* Geodesic : 
    User can give the color to the orbit of the particle moving around the attractor of his/her choice.
    It can be passed while making the object of geodesics_scatter class.
    Default color is "Oranges".
    
    .. code-block:: python
        
        self.cmap_color = cmap_color
        plt.scatter(pos_x, pos_y, s=1, c=time, cmap=self.cmap_color)
        