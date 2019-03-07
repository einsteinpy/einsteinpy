.. einsteinpy

.. image:: http://einsteinpy.github.io/img/logo.png
   :target: http://einsteinpy.github.io/
   :alt: EinsteinPy logo
   :width: 675px
   :align: center

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat-square
   :target: http://www.astropy.org/

.. |mailing| image:: https://img.shields.io/badge/mailing%20list-groups.io-8cbcd1.svg?style=flat-square
   :target: https://groups.io/g/einsteinpy-dev

.. |doi| image:: https://zenodo.org/badge/168302584.svg?style=flat-square
   :target: https://zenodo.org/badge/latestdoi/168302584

.. |gitter| image:: https://badges.gitter.im/EinsteinPy-Project/EinsteinPy.svg
   :alt: Join the chat at https://gitter.im/EinsteinPy-Project/EinsteinPy
   :target: https://gitter.im/EinsteinPy-Project/EinsteinPy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

.. |riotchat| image:: https://img.shields.io/matrix/einsteinpy:matrix.org.svg?style=flat-square
   :target: https://riot.im/app/#/room/#einsteinpy:matrix.org

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
   :target: https://github.com/einsteinpy/einsteinpy/raw/0.1.0/COPYING

.. |docs| image:: https://img.shields.io/badge/docs-v0.1.0-brightgreen.svg?style=flat-square
   :target: https://einsteinpy-project.readthedocs.io/en/v0.1.0/?badge=v0.1.0

.. |travisci| image:: https://img.shields.io/travis/einsteinpy/einsteinpy/0.1.0.svg?style=flat-square&logo=travis
   :target: https://travis-ci.org/einsteinpy/einsteinpy

.. |codeclimate| image:: https://api.codeclimate.com/v1/badges/6efb3f754d20777d8b8d/maintainability
   :target: https://codeclimate.com/github/einsteinpy/einsteinpy/maintainability
   :alt: Maintainability

.. |circleci| image:: https://img.shields.io/circleci/project/github/einsteinpy/einsteinpy/0.1.0.svg?style=flat-square&logo=circleci
   :target: https://circleci.com/gh/einsteinpy/einsteinpy

.. |codecov| image:: https://img.shields.io/codecov/c/github/einsteinpy/einsteinpy.svg?style=flat-square
   :target: https://codecov.io/github/einsteinpy/einsteinpy?branch=0.1.0

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/b95ml54ebspx6sm9?svg=true
   :target: https://ci.appveyor.com/project/shreyasbapat/einsteinpy

.. |orcid-shreyas| image:: https://img.shields.io/badge/id-0000--0002--0870--4665-a6ce39.svg
   :target: https://orcid.org/0000-0002-0870-4665

:Name: EinsteinPy
:Website: https://einsteinpy.github.io/
:Version: 0.1.0

|astropy| |mailing| |gitter| |riotchat| |license| |docs|

|circleci| |travisci| |appveyor| |codecov| |codeclimate|

EinsteinPy is an open source pure Python package dedicated to problems arising in General Relativity and relativistic physics, such as goedesics plotting for schwartzschild space-time model, calculation of schwartzschild radius for any mass given, symbolic calculation of various functions related to GR such as christoffel symbols. Features like visualisation of geodesics of curved black holes and 3D visualisations are some of the features which are planned.
It is released under the MIT license.

Documentation
=============

|docs|

Complete documentation, including a user guide and an API reference, can be read on
the wonderful `Read the Docs`_.

https://einsteinpy-project.readthedocs.io/

.. _`Read the Docs`: https://readthedocs.org/

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

.. _`CONTRIBUTING.rst`: https://github.com/einsteinpy/einsteinpy/blob/0.1.0/CONTRIBUTING.rst

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

.. _`COPYING`: https://github.com/einsteinpy/einsteinpy/blob/0.1.0/COPYING

FAQ
===

What's up with the name?
------------------------

EinsteinPy comes from the name of the famous physicist, nobel laureate, revolutionary person, Prof. Albert Einstein.
This is a small tribute from our part for the amazing work he did for the science.

Can I do <insert awesome thing> with EinsteinPy?
------------------------------------------------

EinsteinPy is focused on general relativity.  One can always discuss probable features on the mailing list and try to implement it.
We welcome every contribution and will be happy to include it in einteinpy.

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
