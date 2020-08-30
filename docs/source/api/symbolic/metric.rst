Metric Tensor Module
====================

This module contains the class for defining a Metric of an arbitrary spacetime, symbolically. Note that, only the coordinate symbols are to 
be supplied to the ``syms`` parameter, while ``arr`` takes the metric (as a SymPy array), which may contain several constants. The symbols in 
``syms`` define the basis to perform several operations on the metric, such as symbolic differentiation. Symbols, for the constants in the 
metric, should be defined independently and used directly in the specification of ``arr``. Please check the metric definitions in 
``einsteinpy.symbolic.predefined`` for examples of doing this.

.. automodule:: einsteinpy.symbolic.metric
    :members:
    :show-inheritance:
