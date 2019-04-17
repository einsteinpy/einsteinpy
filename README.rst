.. einsteinpy

.. image:: https://einsteinpy.org/img/logo.png
   :target: https://einsteinpy.org/
   :alt: EinsteinPy logo
   :width: 675px
   :align: center

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat-square
   :target: http://www.astropy.org/

.. |mailing| image:: https://img.shields.io/badge/mailing%20list-groups.io-8cbcd1.svg?style=flat-square
   :target: https://groups.io/g/einsteinpy-dev

.. |doi| image:: https://zenodo.org/badge/168302584.svg?style=flat-square
   :target: https://zenodo.org/badge/latestdoi/168302584

.. |gitter| image:: https://img.shields.io/gitter/room/EinsteinPy-Project/EinsteinPy.svg?logo=gitter&style=flat-square
   :alt: Join the chat at https://gitter.im/EinsteinPy-Project/EinsteinPy
   :target: https://gitter.im/EinsteinPy-Project/EinsteinPy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

.. |riotchat| image:: https://img.shields.io/matrix/einsteinpy:matrix.org.svg?logo=riot&style=flat-square
   :target: https://riot.im/app/#/room/#einsteinpy:matrix.org

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
   :target: https://github.com/einsteinpy/einsteinpy/raw/master/COPYING

.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat-square
   :target: https://docs.einsteinpy.org/en/latest/?badge=latest

.. |travisci| image:: https://img.shields.io/travis/einsteinpy/einsteinpy/master.svg?style=flat-square&logo=travis
   :target: https://travis-ci.org/einsteinpy/einsteinpy

.. |codeclimate| image:: https://img.shields.io/codeclimate/maintainability/einsteinpy/einsteinpy.svg?logo=code-climate&style=flat-square
   :target: https://codeclimate.com/github/einsteinpy/einsteinpy/maintainability
   :alt: Maintainability

.. |circleci| image:: https://img.shields.io/circleci/project/github/einsteinpy/einsteinpy/master.svg?style=flat-square&logo=circleci
   :target: https://circleci.com/gh/einsteinpy/einsteinpy

.. |codecov| image:: https://img.shields.io/codecov/c/github/einsteinpy/einsteinpy.svg?style=flat-square
   :target: https://codecov.io/github/einsteinpy/einsteinpy?branch=master

.. |appveyor| image:: https://img.shields.io/appveyor/ci/shreyasbapat/einsteinpy.svg?logo=appveyor&style=flat-square
   :target: https://ci.appveyor.com/project/shreyasbapat/einsteinpy

.. |orcid-shreyas| image:: https://img.shields.io/badge/id-0000--0002--0870--4665-a6ce39.svg
   :target: https://orcid.org/0000-0002-0870-4665

:Name: EinsteinPy
:Website: https://einsteinpy.org/
:Version: 0.2.dev0

|astropy| |mailing| |gitter| |riotchat| |license| |docs|

|circleci| |travisci| |appveyor| |codecov| |codeclimate|

EinsteinPy is an open source pure Python package dedicated to problems arising in 
General Relativity and relativistic physics, such as goedesics calculation for vacuum 
solutions for Einstein's field equations, calculation of various quantities in these 
geometries like Schwarzschild Radius and event horizon. The library also has functions 
for Symbolic calculations in GR like Christoffel Symbols and much more is planned. 
The library aims to solve Einstein's field equations for arbitarily complicated 
matter distribution as one of the main goals. 
It is released under the MIT license.

Documentation
=============

|docs|

Complete documentation, including a user guide and an API reference, can be read on
the wonderful `Read the Docs`_.

https://docs.einsteinpy.org/

.. _`Read the Docs`: https://readthedocs.org/

Examples
========

.. |mybinder| image:: https://img.shields.io/badge/launch-binder-e66581.svg?style=flat-square
   :target: https://beta.mybinder.org/v2/gh/einsteinpy/einsteinpy/master?filepath=index.ipynb

|mybinder|

In the examples directory you can find several Jupyter notebooks with specific
applications of einsteinpy. You can consider theses Jupyter Notebooks as tutorials for einsteinpy.
You can launch a cloud Jupyter server using `binder`_ to edit
the notebooks without installing anything. Try it out!

https://beta.mybinder.org/v2/gh/einsteinpy/einsteinpy/master?filepath=index.ipynb

.. _binder: https://beta.mybinder.org/

Requirements
============

EinsteinPy requires the following Python packages:

* NumPy, for basic numerical routines
* Astropy, for physical units and time handling
* numba (optional), for accelerating the code
* matplotlib, for geodesics plotting and visualisations.
* SciPy, for solving ordinary differential equations.
* SymPy (optional), for symbolic calculations related to GR.

EinstienPy is usually tested on Linux, Windows and OS X on Python
3.5, 3.6 and 3.7 against latest NumPy.

==============  ============  ===================
Platform        Site          Status
==============  ============  ===================
Linux           CircleCI      |circleci|
OS X            Travis CI     |travisci|
Windows x64     Appveyor      |appveyor|
==============  ============  ===================

Installation
============

The easiest and fastest way to get the package up and running is to
install EinsteinPy using `conda <http://conda.io>`_::

  $ conda install einsteinpy --channel conda-forge

Please check out the `guide for alternative installation methods`_.

.. _`guide for alternative installation methods`: https://einsteinpy.github.io/installation/

Testing
=======

|codecov|

If installed correctly, the tests can be run using pytest::

  $ python -c "import einsteinpy.testing; einsteinpy.testing.test()"
  ============================= test session starts ==============================
  platform linux -- Python 3.7.1, pytest-4.3.1, py-1.8.0, pluggy-0.9.0
  rootdir: /home/shreyas/Local Forks/einsteinpy, inifile: setup.cfg
  plugins: remotedata-0.3.1, openfiles-0.3.1, doctestplus-0.3.0, cov-2.5.1, arraydiff-0.3
  collected 56 items                                                             
  [...]
  ==================== 56 passed, 1 warnings in 28.19 seconds ====================
  $

Problems
========

If the installation fails or you find something that doesn't work as expected,
please open an issue in the `issue tracker`_.

.. _`issue tracker`: https://github.com/einsteinpy/einsteinpy/issues

Contributing
============

.. image:: https://img.shields.io/waffle/label/einsteinpy/einsteinpy/1%20-%20Ready.svg?style=flat-square
   :target: https://waffle.io/einsteinpy/einsteinpy
   :alt: 'Stories in Ready'

EinsteinPy is a community project, hence all contributions are more than
welcome! For more information, head to `CONTRIBUTING.rst`_.

.. _`CONTRIBUTING.rst`: https://github.com/einsteinpy/einsteinpy/blob/master/CONTRIBUTING.rst

Support
=======

|mailing|

Release announcements and general discussion take place on our `mailing list`_.
Feel free to join!

.. _`mailing list`: https://groups.io/g/einsteinpy-dev

https://groups.io/g/einsteinpy-dev

Please join our `[matrix]`_ channel or `gitter`_ chat room for further queries.

.. _`[matrix]`: https://matrix.to/#/#einsteinpy:matrix.org

.. _`gitter`: https://gitter.im/EinsteinPy-Project/EinsteinPy

Citing
======

If you use EinsteinPy on your project, please
`drop us a line <mailto:einsteinpy.project@gmail.com>`_.

You can also use the DOI to cite it in your publications. This is the latest
one:

|doi|

And this is an example citation format::

 Shreyas Bapat et al.. (2019). EinsteinPy: einsteinpy 0.1.0. Zenodo. 10.5281/zenodo.2582388


License
=======

|license|

EinsteinPy is released under the MIT license, hence allowing commercial
use of the library. Please refer to `COPYING`_.

.. _`COPYING`: https://github.com/einsteinpy/einsteinpy/blob/master/COPYING

FAQ
===

What's up with the name?
------------------------

EinsteinPy comes from the name of the famous physicist, nobel laureate, revolutionary person, Prof. Albert Einstein.
This is a small tribute from our part for the amazing work he did for the science.

Can I do <insert awesome thing> with EinsteinPy?
------------------------------------------------

EinsteinPy is focused on general relativity.  One can always discuss probable features on the mailing list and try to implement it.
We welcome every contribution and will be happy to include it in EinsteinPy.

What's the future of the project?
---------------------------------

EinsteinPy is actively maintained and we hope to receive an influx of new contributors.
The best way to get an idea of the roadmap is to see the `Milestones`_ of
the project.

.. _`Milestones`: https://github.com/einsteinpy/einsteinpy/milestones

Inspiration
-----------

The whole documentation, and code structure is shamelessly inspired by `poliastro`_ . We really thank the developers to
help us acheive this.

.. _`poliastro`: https://docs.poliastro.space/
