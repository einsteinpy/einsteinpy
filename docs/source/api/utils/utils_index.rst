Utils
=====

This module provides common utility classes and functions used throughout the EinsteinPy 
project. These utilities support the core functionality of the library and provide 
specialized mathematical tools for relativistic calculations.

Key utilities include:

- **DualNumber**: A class implementing dual numbers for Forward Mode Automatic Differentiation. 
  This is used for computing derivatives efficiently, which is particularly useful in 
  numerical relativity for gradient-based computations.
- **ScalarFactor**: Utilities for working with scale factors in cosmology
- **Exceptions**: Custom exception classes for better error handling and debugging

These utilities form the foundation for many higher-level operations in EinsteinPy.

.. toctree::
    :maxdepth: 2

    dual
    scalar_factor
    exceptions
