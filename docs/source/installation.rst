Installation
============

Stable version
~~~~~~~~~~~~~~

You can install the most recent ``stable`` version of EinsteinPy in the 
following ways. At the moment, using ``pip`` is recommended, since 
the package on ``conda-forge`` is currently a version behind ``PyPI``. 
We are working on updating it. 

- Using ``pip``:

  .. code-block:: sh

       $ pip install einsteinpy

- Using ``conda``:

  .. code-block:: sh

       $ conda install -c conda-forge einsteinpy

Latest version
~~~~~~~~~~~~~~

For installing the ``latest`` version, you can use any of the following methods:

- Installation from clone:

  .. code-block:: sh

       $ git clone https://github.com/einsteinpy/einsteinpy.git
       $ cd einsteinpy/
       $ pip install .

- Install using pip:

  .. code-block:: sh

       $ pip install git+https://github.com/einsteinpy/einsteinpy.git

**Note that the latest version is not guaranteed to be stable.**

Development version
~~~~~~~~~~~~~~~~~~~

If you want to contribute to EinsteinPy, please see the `Developer Guide`_ 
for instructions on how to set up your development environment and start 
contributing.

.. _Developer Guide: https://docs.einsteinpy.org/en/latest/dev_guide.html

----

First run
~~~~~~~~~

Now that you have EinsteinPy installed, you can try out some of the 
examples listed in the `Examples`_ section. Below, we demonstrate a simple 
example of plotting a precessing timelike geodesic in Schwarzschild spacetime.

.. _`Examples` : https://docs.einsteinpy.org/en/latest/jupyter.html

.. code-block:: python

   from einsteinpy.plotting import StaticGeodesicPlotter
   from einsteinpy.examples import precession

   gpl = StaticGeodesicPlotter()
   gpl.plot2D(precession())
   gpl.show()

.. figure:: _static/precess.png
   :align: center
   :figwidth: 450px
   :alt: Precession in Schwarzschild spacetime
