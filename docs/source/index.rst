EinsteinPy - Making Einstein possible in Python
===============================================

.. image:: https://einsteinpy.org/img/wordmark.png
   :target: https://einsteinpy.org/
   :alt: EinsteinPy logo
   :width: 675px
   :align: center

**EinsteinPy** is an open-source pure Python package dedicated to the study of problems arising 
in General Relativity and gravitational physics. Using EinsteinPy, it is possible to approach 
problems symbolically as well as numerically.

On the symbolic side, EinsteinPy provides a robust API, which allows users to access several 
predefined metrics or to define custom metrics and perform symbolic calculations on them. 
Computation of quantities, such as terms of metric tensors, Christoffel symbols, Riemann 
Curvature tensor, Ricci tensor, stress-energy tensor and more are all supported and extensible, 
since the symbolic modules are built on top of SymPy.

On the numerical side, EinsteinPy provides tools to calculate and visualize geodesics, metric 
singularities and hypersurface embeddings in certain spacetimes. We hope to extend the package to 
include more features in the future.

EinsteinPy works on recent versions of Python and is released under the MIT license, 
hence allowing commercial use of the library.

Motivation
~~~~~~~~~~

Most of the general relativity packages currently available in the
gravitational physics research community demand knowledge of 
paid-for software / computer algebra systems or programming experience 
in languages like C/C++. Many of the people working in the field of 
gravitational physics come from a theoretical background and have 
little programming experience. As such, they find using these libraries 
difficult. EinsteinPy is motivated by these problems and aims to provide 
a high level of abstraction, that enables anyone to do general relativity 
calculations without being bogged down by the implementation details or 
obscure errors. Since the library is written in Python, it is easy to 
use and understand, which also makes it a perfect fit for teaching purposes.

Moreover, EinsteinPy being `FOSS <https://en.wikipedia.org/wiki/Free_and_open-source_software>`_ 
also helps interested users understand the code and modify it to suit 
their needs, which in turn grows the project and positively adds to the 
research community.

----

.. _`GitHub`: https://github.com/einsteinpy/einsteinpy
.. _`report errors or suggest improvements`: https://github.com/einsteinpy/einsteinpy/issues
.. _`developer guide`: https://docs.einsteinpy.org/en/latest/dev_guide.html
.. _`Gitter`: https://gitter.im/EinsteinPy-Project/EinsteinPy
.. _`Element`: https://app.element.io/#/room/#einsteinpy:matrix.org
.. _`mailing list`: https://groups.io/g/einsteinpy-dev

.. |gitter| image:: https://img.shields.io/gitter/room/EinsteinPy-Project/EinsteinPy.svg?logo=gitter&style=flat-square
   :alt: Join the chat at https://gitter.im/EinsteinPy-Project/EinsteinPy
   :target: https://gitter.im/EinsteinPy-Project/EinsteinPy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

.. |element| image:: https://img.shields.io/matrix/einsteinpy:matrix.org.svg?logo=riot&style=flat-square
   :target: https://app.element.io/#/room/#einsteinpy:matrix.org

EinsteinPy is developed by an open community. Hence, all contributions are more than
welcome! Visit our `GitHub`_ to view EinsteinPy's source code and to `report errors or suggest 
improvements`_. For information on contributing code to the project, head to the `developer guide`_.

|element| |gitter|

Please join our `Element`_ channel or `Gitter`_ chat room for general discussions 
and further queries. Release announcements take place on our `mailing list`_.

.. include:: form.rst

If you still have a doubt, write to us directly at `all@einsteinpy.org <mailto:all@einsteinpy.org>`_.

Please use the links below or on the navigation bar to explore some examples or the API 
documentation. If you want to view the documentation for a specific version, please use the
pop-up menu on the bottom-right corner of the page.

----

Contents
--------

.. toctree::
    :maxdepth: 2

    getting_started
    user_guide
    metric
    jupyter
    api/index
    dev_guide
    changelog
    codeofconduct
