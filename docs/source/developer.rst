###############
Developer Guide
###############

This guide provides onboarding information to external contributors and serve as a reference for maintainers.

First-Time Contributors
=======================

DendroPy welcomes contributions from users to help improve the library and make it more useful for the community.
Thanks for getting involved!

General information for first-time open source contributors can be found `here <https://opensource.guide/how-to-contribute/>`_ and `here <https://github.com/firstcontributions/first-contributions>`_.

Filing a bug report or making a feature request --- covered elsewhere in the documentation --- can also be a great way to contribute to the project.

Linting
=======

We use a static linter to help detect source quality issues like syntax errors in the code base.
To perform linting,

.. code-block:: shell

   ./lint.sh

We use `Ruff <https://github.com/charliermarsh/ruff>`_ for linting.
Instructions to install it can be found `here <https://github.com/charliermarsh/ruff#getting-started>`_.

Testing
=======

To run all tests,

.. code-block:: shell

   python setup.py test

To run an individual test file,

.. code-block:: shell

   python -m unittest tests/test_example.py

Documentation
=============

We use `Netlify <https://netlify.com/>`_ to host the projects' documentation.
We keep deployed documentation up-to-date with the `main` branch through our Github Actions continuous integration.

Manual documentation builds aren't usually necessary, unless you're encountering a documentation-specific issue.
To build the documentation locally,

.. code-block:: shell

   make -C docs/ html

This requires the version of DendroPy you want to build for and requirements listed in `docs/requirements.txt` to be installed.

Docstrings
==========

DendroPy uses `NumPy docstring conventions <https://numpydoc.readthedocs.io/en/latest/format.html>`_ with `Sphinx cross-referencing syntax <https://www.sphinx-doc.org/en/master/usage/restructuredtext/domains.html#cross-referencing-syntax>`_.

Virtual Environments
====================

Setting up a virtual environment can help you manage dependencies and isolate project-specific packages from your system's global packages.
This can prevent conflicts and ensure that your project runs smoothly and consistently on different machines.

To set up and enter a virtual environment

.. code-block:: shell

   python -m venv env
   source env/bin/activate

Continuous Integration
======================

CI will run tests and linting each time the main branch or an open pull request is pushed to.
Note that tests are run across all supported Python versions.
So, even if tests pass locally CI may flag incompatibilities with other Python versions.

Continuous integration configuration is handled in `.github/workflows/`.

Version Bumping
===============

We use the `BumpVer <https://github.com/mbarkhau/bumpver>`_ tool to perform version bumping.
To update the project to a new version, run of the following

.. code-block:: shell

   bumpver update --patch
   bumpver update --minor
   bumpver update --major

This will create a tagged commit that updates the source (i.e., `__version__`) with the new version number

Only perform version bumping on the `main` branch.
Otherwise, you could create an inadvertent release.

Releasing to PyPi
=================

Continuous integration will automatically release any git tag that begins with `v`.
So, after version bumping simply

.. code-block:: shell

   git push origin --tags

to begin the release process.

Continuous integration will take it from there!
After quality control checks pass (i.e., the usual continuous integration), the new package version will be published to PyPi.
Note, though, that if quality control checks fail the release will be aborted.
