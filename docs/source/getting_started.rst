Getting started
===============

Overview
--------

EinsteinPy is a easy-to-use python library which provides a
user-friendly interface for supporting numerical relativity and
relativistic astrophysics research. The library is an attempt to provide
programming and numerical environment for a lot of numerical relativity
problems like geodesics plotter, gravitational lensing and ray tracing,
solving and simulating relativistic hydrodynamical equations, plotting
of black hole event horizons, solving Einstein’s field equations and
simulating various dynamical systems like binary merger etc.


Who can use?
------------

Most of the numerical relativity platforms currently available in the
gravitational physics research community demands a heavy programming
experience in languages like C, C++ or their wrappers on some other non
popular platforms. Many of the people working in the field of
gravitational physics have theoretical background and does not have any
or have little programming experience and they find using these
libraries mind-boggling. EinsteinPy is motivated by this problem and
provide a high level of abstraction that shed away from user all the
programming and algorithmic view of the implemented numerical methods
and enables anyone to simulate complicated system like binary merger
with just 20-25 lines of python code.

Even people who does not know any python programming can also follow up
with the help of tutorials and documentation given for the library. We
aim to provide all steps, from setting up your library environment to
running your first geodesic plotter with example jupyter notebooks.

So now you are motivated enough so let’s first start with installing the
library.


Installation
------------

It’s as easy as running one command!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stable Versions:
~~~~~~~~~~~~~~~~

For installation of the latest ``stable`` version of EinsteinPy:

- Using pip:

  .. code-block:: sh

       $ pip install einsteinpy

- Using conda:

  .. code-block:: sh

       $ conda install -c conda-forge einsteinpy

Latest Versions
~~~~~~~~~~~~~~~

For installing the development version, you can do two things:

- Installation from clone:

  .. code-block:: sh

       $ git clone https://github.com/einsteinpy/einsteinpy.git
       $ cd einsteinpy/
       $ python setup.py install

- Install using pip:

  .. code-block:: sh

       $ pip install git+https://github.com/einsteinpy/einsteinpy.git

Development Version
~~~~~~~~~~~~~~~~~~~

  .. code-block:: sh

       $ git clone your_account/einsteinpy.git
       $ pip install --editable /path/to/einsteinpy[dev]

Please open an issue `here`_ if you feel any
difficulty in installation!

.. _`here` : https://github.com/einsteinpy/einsteinpy/issues


Running your first code using the library
-----------------------------------------

Various examples can be found in the `examples`_ folder.

.. _`examples` : https://einsteinpy-project.readthedocs.io/en/latest/jupyter.html

Contribute
----------

EinsteinPy is an open source library which is under heavy development.
To contribute kindly do visit :

https://github.com/einsteinpy/einsteinpy/

and also check out current posted issues and help us expand this
awesome library.
