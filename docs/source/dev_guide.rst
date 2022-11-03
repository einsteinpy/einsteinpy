.. Keep CONTRIBUTING.rst and ./docs/sources/dev_guide.rst synchronized

Developer Guide
===============

As with any open-source project, we rely on the community for implementing 
new ideas and fixing old bugs. You can contribute to EinsteinPy in many ways, 
for example, by reporting errors in code or documentation, writing new code, 
improving documentation, or even by writing tutorials or blog posts. If you 
are interested in contributing to EinsteinPy, please go through the following 
sections.


Contributing
------------

For starters, we recommend checking out the `"good first issue"`_ tag on 
the issue tracker. Those issues should be relatively easy to fix 
and would require minimal knowledge about the library. However, you 
will need to possess some knowledge about how ``git`` works. ``git`` 
is a decentralized version control system that preserves codebase history, 
keeps track of changes and allows for collaboration. If you are new to 
``git`` and version control, try following the `GitHub Quickstart`_ tutorial.

If you already know all this and would like to contribute, then 
that's awesome! But before coding away, please make sure your 
intended addition or change is in line with the project's scope and goals. 
You can do this by chatting with us on `Element`_ or `Gitter`_, or `mailing`_ us.

All code changes & additions should be accompanied by altered or new 
unit-tests, so that the code coverage increases or stays the same. Automated 
services will ensure that your code works across the various platforms that 
EinsteinPy supports.

.. _`"good first issue"`: https://github.com/einsteinpy/einsteinpy/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22
.. _`GitHub Quickstart`: https://docs.github.com/en/get-started/quickstart/hello-world
.. _`mailing`: all@einsteinpy.org
.. _`Element`: https://app.element.io/#/room/#einsteinpy:matrix.org
.. _`Gitter`: https://gitter.im/EinsteinPy-Project/EinsteinPy


Development Pipeline
--------------------
Note that EinsteinPy is a Python 3-only library. So, you need to have Python 3 
`installed <https://www.python.org/downloads/>`_ in your system and you need 
to be familiar with it. Then, you can follow these steps to set up a 
platform-agnostic development environment and start contributing to EinsteinPy:

1. Install `git <https://git-scm.com/>`_.
2. Create an account on `GitHub <https://github.com/signup>`_, if you don't already have one.
3. `Fork`_ EinsteinPy's repository (``einsteinpy/einsteinpy``).
4. `Clone`_ your forked repository.
5. Create a virtual environment and activate it using the following commands:

.. code:: bash

   $ python -m venv <your-venv-name>
   $ source <your-venv-name>/bin/activate # Linux/macOS
   $ ./<your-venv-name>/Scripts/activate # Windows

Learn more about Python virtual environments `here <https://docs.python.org/3/tutorial/venv.html>`_.

6. After activating the environment, install EinsteinPy in editable or development mode like so: 

.. code:: bash

   $ pip install -e ./einsteinpy/[dev]
 
``[dev]`` ensures that all the dependencies required for development and 
testing are installed, while the ``-e`` or ``--editable`` flag ensures that any changes you make to the code 
take effect immediately (after you save the changes).

7. Create and switch to a new branch like so:

.. code:: bash

   $ git checkout -b <your-branch-name>

Usually, the branch name should be kept similar to the issue/topic you are working on.

8. Make changes to the code and add unit-tests (if applicable). Ensure that the code is formatted properly. See the `Code Linting & Testing`_ section below.
9. After you are done, commit your changes and `push`_ them to your forked repository.
10. Finally, open a `Pull Request`_ (PR) to the `main` branch of the ``einsteinpy/einsteinpy`` repository. You can also `request a review`_ from specific maintainers of the project. The maintainers will review your PR and might ask you to make some changes. If there are any changes required, you can make them in the same branch and push them to your forked repository. The PR will automatically update with the new changes.
11. When your PR is approved and ready to be merged, we will ask you to update the project `CHANGELOG`_ and if needed, we might also ask you to `squash your commits`_.

If you are facing issues with any of the above steps, please feel free to ask for help on `Element`_ or `Gitter`_.

.. _`Pull Request`: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request?tool=webui>`_ 
.. _`Fork`: https://docs.github.com/en/get-started/quickstart/fork-a-repo?tool=webui
.. _`Clone`: https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository?tool=webui>
.. _`push`: https://docs.github.com/en/get-started/using-git/pushing-commits-to-a-remote-repository
.. _`request a review`: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/requesting-a-pull-request-review
.. _`CHANGELOG`: https://github.com/einsteinpy/einsteinpy/blob/main/CHANGELOG
.. _`squash your commits`: https://github.com/einsteinpy/einsteinpy/wiki/How-to-refactor-commits


Code Linting & Testing
----------------------

After you push your code and open a PR to the `main` branch, your code will be 
tested automatically by the continuous integration (CI) tools. If they encounter any errors, 
you will be notified by CI with the details about the error, which 
you can use to make necessary changes to fix the error. These errors usually come 
in the form of code quality errors, failed tests and decreased coverage.

Code quality
~~~~~~~~~~~~

`Code quality` is an important aspect of any software project. It ensures that
the code is readable, maintainable, and bug-free. The quality of a piece of code 
depends on its formatting, style and complexity. To maintain consistent code 
formatting and style, we use `black`_, `isort`_ and `mypy`_. For convenience, 
we have set up `tox`_, which runs all these in a single, short command:

.. code:: bash

    $ cd ./einsteinpy/
    $ tox -e reformat

If your PR is failing quality checks, executing the above commands should fix it. 
Also, to ensure that the code remains understandable and maintainable over time, 
we enforce `Cyclomatic Complexity <https://docs.codeclimate.com/docs/cyclomatic-complexity>`_ (CC) checks on the codebase using CodeFactor/CodeClimate. 
CC is a measure of the complexity of a program. The lower the CC, the easier it is to 
understand and maintain the code. If your PR is running into CodeFactor or CodeClimate 
errors, you will have to refactor your code, so that its CC is below a certain threshold. 
To check the CC of your code locally, you can use ``radon``:

.. code:: bash

   $ radon cc ./einsteinpy/

Unit-tests & Code coverage
~~~~~~~~~~~~~~~~~~~~~~~~~~

Since you have made changes to the codebase, it is likely that some of the 
unit-tests that were previously passing will now fail. CI will alert you to these 
failures on your PR with the error details. You can use these details to 
fix the errors, commit & push the changes and make the tests pass. In case you want to 
debug the test errors locally, you can do so by running the following command:

.. code:: bash

   $ pytest --cov-report=term-missing --cov=einsteinpy ./einsteinpy/

Note that this can take a while. If the error is isolated to only one file or a few files,
you can choose to test only that or those file(s). For example, to execute all the 
``einsteinpy.metric``-related tests, you can run:

.. code:: bash

   $ pytest --cov-report=term-missing --cov=einsteinpy ./einsteinpy/tests/test_metric/

This command also reports the overall code coverage after all the tests are run. `Code coverage` 
is a measure of how much of the codebase is covered by the tests. It is a good 
practice to have a high code coverage, so that the tests can catch any bugs that 
might be introduced in the codebase. We use `codecov`_ to track the code coverage 
of the project. You can see the current `code coverage of the project
here <https://codecov.io/gh/einsteinpy/einsteinpy>`_. If your PR is failing coverage 
checks, you will have to add more tests to increase the coverage.

.. _`black`: https://black.readthedocs.io/en/stable/
.. _`isort`: https://pycqa.github.io/isort/
.. _`mypy`: https://mypy.readthedocs.io/en/latest/?badge=latest
.. _`tox`: https://tox.wiki/en/latest/
.. _`codecov`: https://about.codecov.io/


Documentation
-------------

After you have implemented your bugfix or your shiny new feature, you should also 
add some documentation, so that users, maintainers and future contributors can 
understand how to use or make changes to your code. All of EinsteinPy's non-API 
documentation is stored in text files under `docs/source <https://github.com/einsteinpy/einsteinpy/tree/main/docs/source>`_. 
If you think anything can be improved there, please edit the files and open a PR. 
The API docs comprising the docstrings of the Python code follow 
`numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html>`_ guidelines. 
If you come across any inconsistency or opportunity for improvement, feel free to 
edit the docstrings and submit a PR.

We use `Sphinx`_ to generate the overall API + non-API documentation, which is then 
`hosted here <https://docs.einsteinpy.org/en/latest/>`_, courtesy of `Read the Docs`_. To 
build the docs locally, you can run the following commands:

.. code:: bash

   $ pip install Sphinx nbsphinx jupyter_sphinx
   $ cd ./einsteinpy/docs/
   $ sphinx-build -b html source build
   $ cd ./build/
   $ python -m http.server 8000 --bind 127.0.0.1

This should open the built documentation website in a web browser at 
``http://127.0.0.1:8000``. If you want hot-reloading, you can use 
``sphinx-autobuild`` instead of ``sphinx-build``, after installing it using 
``pip install sphinx-autobuild``. This will automatically rebuild the docs and 
refresh the browser tab whenever you make changes to the source files.

In addition to the usual documentation, the `GitHub Wiki`_ for EinsteinPy is open to everybody. Please feel free to add
new content there.

.. _`Sphinx`: https://www.sphinx-doc.org/en/master/
.. _`Read the Docs`: https://readthedocs.org/
.. _`GitHub Wiki`: https://github.com/einsteinpy/einsteinpy/wiki


After your PR is merged
-----------------------

Great job 🎊🎊. Your PR just got approved & merged. Now how do you ensure your local 
``main`` branch is up-to-date with ``upstream/main``? Like so:

.. code:: bash

   $ git remote add upstream https://github.com/einsteinpy/einsteinpy.git # Set up upstream
   $ git checkout main
   $ git fetch upstream
   $ git merge upstream/main
   $ git branch -d <your-branch-name> # Delete your branch
   $ git push origin main

And now, you are all set to create another branch and contribute further to EinsteinPy. Have fun coding 😀!

Reporting bugs or suggestions
-----------------------------

Not only can things break at any given time, but different people also have different
use cases for the project. If you find anything that doesn't work as expected
or have suggestions, please refer to the `issue tracker`_ on GitHub.

.. _`issue tracker`: https://github.com/einsteinpy/einsteinpy/issues
