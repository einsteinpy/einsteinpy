Static Geodesic Plotting
========================

This module contains the methods for static geodesic plotting.

.. automodule:: einsteinpy.plotting.senile.geodesics_static
    :members:

Auto and Manual scaling
=======================
EinsteinPy supports Automatic and Manual scaling of the attractor to make plots look better since radius of 
attractor can be really small and not visible.

* Manual_Scaling : 
    If the user provides the attractor_radius_scale, then the autoscaling will not work.
    This is checked by initialising the  attractor_radius_scale by -1 and if the user enters the value then 
    it will be >0 so the value won't remain -1 which is easily checked.
    
    The radius is multiplied to the value given in attractor_radius_scale
    
    .. code-block:: python

        radius = radius * self.attractor_radius_scale

* Auto Scaling : 
    If the user does not provide the attractor_radius_scale, the value will be initialised to -1 and 
    then we will call the auto scaling function. In autoscaling, the attractor radius is first 
    initialised to the minimum distance between the attractor and the object moving around it. 
    Now, if this radius is greater than the 1/12th of minimum of range of X and Y coordinates then,
    the radius is initialised to this minimum. This is done so that the plots are easy to look at.

    minrad_nooverlap : Stores the minimum distance between the particle and attractor
    .. code-block:: python

        for i in range(0, len(self.xarr)):
            minrad_nooverlap = min( minrad_nooverlap, self.mindist(self.xarr[i], self.yarr[i]))

    minlen_plot : Stores the minimum of range of X and Y axis
    
    .. code-block:: python
            
            xlen = max(self.xarr) - min(self.xarr)
            ylen = max(self.yarr) - min(self.yarr)
            minlen_plot = min(xlen, ylen)
    
    multiplier : Stores the value which is multiplied to the radius to make it 1/12th of the minlenplot
    
    .. code-block:: python
            
            mulitplier = minlen_plot / (12 * radius)
            min_radius = radius * mulitplier


Attractor Color
===============
* Color Options : 
    User can give the color to attractor of his/her choice.
    It can be passed while calling the geodesics_static class.
    Default color of attractor is "black".

    .. code-block:: python
        
        self.attractor_color = attractor_color
        mpl.patches.Circle( (0, 0), radius.value, lw=0, color=self.attractor_color)
