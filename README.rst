.. einsteinpy

.. image:: https://blog.einsteinpy.org/img/logo.png
   :target: https://einsteinpy.org/
   :alt: EinsteinPy logo
   :width: 675px
   :align: center

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat-square
   :target: http://www.astropy.org/

.. |mailing| image:: https://img.shields.io/badge/mailing%20list-groups.io-8cbcd1.svg?style=flat-square
   :target: https://groups.io/g/einsteinpy-dev

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2582387.svg
   :target: https://doi.org/10.5281/zenodo.2582387

.. |gitter| image:: https://img.shields.io/gitter/room/EinsteinPy-Project/EinsteinPy.svg?logo=gitter&style=flat-square
   :alt: Join the chat at https://gitter.im/EinsteinPy-Project/EinsteinPy
   :target: https://gitter.im/EinsteinPy-Project/EinsteinPy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

.. |riotchat| image:: https://img.shields.io/matrix/einsteinpy:matrix.org.svg?logo=riot&style=flat-square
   :target: https://app.element.io/#/room/#einsteinpy:matrix.org

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
   :target: https://github.com/einsteinpy/einsteinpy/raw/main/COPYING

.. |docslatest| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat-square
   :target: https://docs.einsteinpy.org/en/latest/?badge=latest

.. |docsstable| image:: https://img.shields.io/badge/docs-stable-brightgreen.svg?style=flat-square
   :target: https://docs.einsteinpy.org/en/stable/?badge=stable

.. |circleci| image:: https://img.shields.io/circleci/project/github/einsteinpy/einsteinpy/main.svg?style=flat-square&logo=circleci
   :target: https://circleci.com/gh/einsteinpy/einsteinpy

.. |ghactions| image:: https://img.shields.io/github/workflow/status/einsteinpy/einsteinpy/Test_MacOS?logo=github&style=flat-square
   :target: https://github.com/einsteinpy/einsteinpy/actions?query=branch%3Amain

.. |codecov| image:: https://img.shields.io/codecov/c/github/einsteinpy/einsteinpy.svg?style=flat-square
   :target: https://codecov.io/github/einsteinpy/einsteinpy?branch=main

.. |appveyor| image:: https://img.shields.io/appveyor/ci/shreyasbapat/einsteinpy/main?logo=appveyor&style=flat-square
   :target: https://ci.appveyor.com/project/shreyasbapat/einsteinpy/branch/main

.. |orcid-shreyas| image:: https://img.shields.io/badge/id-0000--0002--0870--4665-a6ce39.svg
   :target: https://orcid.org/0000-0002-0870-4665

:Name: EinsteinPy
:Website: https://einsteinpy.org/
:Version: 0.5.dev0

|astropy| |docslatest| |gitter| |riotchat| |mailing| |license| 

|circleci| |ghactions| |appveyor| |codecov|

**EinsteinPy** is an open-source pure Python package dedicated to the study of problems arising 
in General Relativity and gravitational physics. Using EinsteinPy, it is possible to approach 
problems symbolically as well as numerically. On the symbolic side, EinsteinPy provides a robust 
API, which allows users to access several predefined metrics or to define custom metrics and perform 
symbolic calculations on them. Computation of quantities, such as terms of metric tensors, 
Christoffel symbols, Riemann Curvature tensor, Ricci tensor, stress-energy tensor and more are 
all supported and extensible, since the symbolic modules are built on top of SymPy. On 
the numerical side, EinsteinPy provides tools to calculate and visualize geodesics, metric 
singularities and hypersurface embeddings in certain spacetimes. We hope to extend the package to 
include more features in the future. EinsteinPy is released under the MIT license.

Documentation
=============

|docslatest| |docsstable|

Complete documentation, including a user guide and an API reference, can be perused on
the wonderful `Read the Docs`_.

.. _`Read the Docs`: https://docs.einsteinpy.org/en/latest/

Examples
========

.. |mybinder| image:: https://img.shields.io/badge/launch-binder-e66581.svg?style=flat-square
   :target: https://mybinder.org/v2/gh/einsteinpy/einsteinpy/main?filepath=index.ipynb

|mybinder|

Several tutorial Jupyter notebooks on specific applications of EinsteinPy can be found 
in the `examples`_ directory. You can launch a Jupyter notebook instance in the cloud 
using `binder`_ to run and edit these notebooks without installing anything. Try it out!

.. _examples: https://github.com/einsteinpy/einsteinpy/tree/main/docs/source/examples
.. _binder: https://mybinder.org/v2/gh/einsteinpy/einsteinpy/main?filepath=index.ipynb

Requirements
============

EinsteinPy requires the following Python packages:

* ``NumPy``, for basic numerical routines
* ``SciPy``, for solving ordinary differential equations
* ``SymPy``, for symbolic calculations
* ``Astropy``, for handling conversion between physical units
* ``Matplotlib``, for producing static visualizations
* ``Plotly``, for producing interactive visualizations
* ``Numba``, for accelerating the code

EinsteinPy is currently tested on Linux, Windows and macOS on Python 3.7 and 3.8, against the latest ``NumPy``.

==============  ===============  ===================
Platform        Site             Status
==============  ===============  ===================
Linux           CircleCI         |circleci|
macOS           Github Actions   |ghactions|
Windows x64     Appveyor         |appveyor|
==============  ===============  ===================

Installation
============

Currently, the recommended way to install EinsteinPy is using ``pip``
from `PyPI <https://pypi.org/project/einsteinpy/>`_::

  $ pip install einsteinpy

Or, you can install the package using `conda <https://anaconda.org/conda-forge/einsteinpy>`_::

  $ conda install einsteinpy --channel conda-forge

Note that the package on ``conda-forge`` is currently a version behind ``PyPI``. We are working on updating it. 

For Debian/Ubuntu/Mint users, the package is installable via `apt <https://packages.debian.org/sid/python3-einsteinpy>`_ (Ubuntu 19.04 onwards)::

  $ sudo apt install python3-einsteinpy

If you prefer to install from source to stay on the latest but likely unstable version, 
you can do so using the method described `here <https://docs.einsteinpy.org/en/latest/getting_started.html#installation>`_.


Problems
========

If the installation fails or you find something that doesn't work as expected,
please open an issue in the `issue tracker`_.

.. _`issue tracker`: https://github.com/einsteinpy/einsteinpy/issues

Contributing
============

EinsteinPy is a community project. Hence, all contributions are more than
welcome! For more information, head to `CONTRIBUTING`_ or see the `developer guide`_.

.. _`CONTRIBUTING`: https://github.com/einsteinpy/einsteinpy/blob/main/CONTRIBUTING.rst
.. _`developer guide`: https://docs.einsteinpy.org/en/latest/dev_guide.html

Support
=======

|gitter| |riotchat| |mailing|

Please join our `[matrix]`_ channel or `Gitter`_ chat room for general discussions and further queries.

.. _`[matrix]`: https://matrix.to/#/#einsteinpy:matrix.org
.. _`Gitter`: https://gitter.im/EinsteinPy-Project/EinsteinPy

Release announcements take place on our `mailing list`_. Feel free to join!

.. _`mailing list`: https://groups.io/g/einsteinpy-dev

If you still have a doubt, write to us directly at `all@einsteinpy.org <mailto:all@einsteinpy.org>`_.

Citing
======

If you use EinsteinPy in your project, please `drop us a line <mailto:all@einsteinpy.org>`_. 
You can also use the DOI to cite it in your publications. This is the latest one:

|doi|

And this is an example citation format:

 Shreyas Bapat et al (2021). EinsteinPy 0.4.0 (v0.4.0). Zenodo. https://doi.org/10.5281/zenodo.2582387

License
=======

|license|

EinsteinPy is released under the MIT license, thereby allowing commercial
use of the library. Please refer to `COPYING`_ for more details.

.. _`COPYING`: https://github.com/einsteinpy/einsteinpy/blob/main/COPYING

FAQ
===

Why "EinsteinPy"?
-----------------

EinsteinPy borrows the name of the famous physicist, Nobel laureate and revolutionary 
human, Dr. Albert Einstein. This is a small tribute on our part for the amazing work 
he did for humanity!


Can I do <`insert nerdy thing`> with EinsteinPy?
------------------------------------------------

EinsteinPy is focused on general relativity. One can always discuss probable features in discussion 
forums and the mailing list and also work with the maintainers to try to implement them. 
We welcome every contribution to EinsteinPy. Please see `CONTRIBUTING`_ for more details.

What's the future of the project?
---------------------------------

EinsteinPy is actively maintained and we hope to receive an influx of new contributors.
The best way to get an idea about the roadmap is to view the `milestones`_ of
the project.

.. _`Milestones`: https://github.com/einsteinpy/einsteinpy/milestones

Inspiration
-----------

The documentation and code structure is shamelessly inspired by `poliastro`_. We wholeheartedly thank the ``poliastro``
developers that made this possible. EinsteinPy is nothing without its supporters and community.

.. _`poliastro`: https://docs.poliastro.space/
