.. raw:: html

   <div align="center">

|logo|

.. |logo| image:: https://raw.githubusercontent.com/jeetsukumaran/DendroPy/DendroPy4/doc/source/_static/dendropy_logo.png
   :target: https://github.com/jeetsukumaran/DendroPy
   :alt: DendroPy wordmark

+--------------------+-+---------------------+-+---------------------+-+-------------------+
| *automated tests:* | | *package version:*  | | *documentation:*    | | *test coverage:*  |
|                    | |                     | |                     | |                   |
| |CI|               | | |PyPI|              | | |Docs|              | | |Coverage|        |
+--------------------+-+---------------------+-+---------------------+-+-------------------+

.. |CI| image:: https://github.com/jeetsukumaran/DendroPy/actions/workflows/ci.yaml/badge.svg
   :target: https://github.com/jeetsukumaran/DendroPy/actions/workflows/ci.yaml
   :alt: Continuous Integration build
.. |PyPI| image:: https://img.shields.io/pypi/v/DendroPy.svg
   :target: https://pypi.org/project/DendroPy/
   :alt: PyPI version
.. |Docs| image:: https://github.com/jeetsukumaran/DendroPy/actions/workflows/pages/pages-build-deployment/badge.svg
   :target: https://jeetsukumaran.github.io/DendroPy/
   :alt: Documentation Status
.. |Coverage| image:: https://codecov.io/gh/jeetsukumaran/DendroPy/graph/badge.svg?token=JwMfFOpBBD
   :target: https://codecov.io/gh/jeetsukumaran/DendroPy
   :alt: codecov coverage
.. |nbsp| unicode:: 0xA0
   :trim:

.. raw:: html

   </div>


DendroPy is a Python library for phylogenetic computing.
It provides classes and functions for the simulation, processing, and
manipulation of phylogenetic trees and character matrices, and supports the
reading and writing of phylogenetic data in a range of formats, such as NEXUS,
NEWICK, NeXML, Phylip, FASTA, etc.  Application scripts for performing some
useful phylogenetic operations, such as data conversion and tree posterior
distribution summarization, are also distributed and installed as part of the
libary.  DendroPy can thus function as a stand-alone library for phylogenetics,
a component of more complex multi-library phyloinformatic pipelines, or as a
scripting "glue" that assembles and drives such pipelines.

The primary home page for DendroPy, with detailed tutorials and documentation, is at:

    https://jeetsukumaran.github.io/DendroPy/

DendroPy is also hosted in the official Python repository:

    http://pypi.org/project/DendroPy/

Requirements and Installation
=============================

The current version of DendroPy requires Python 3.

You can install DendroPy by running::

    $ python -m pip install dendropy

For Conda users, DendroPy can be installed from the `conda-forge <https://conda-forge.org/>`_ channel::

    $ conda install -c conda-forge dendropy

More information is available here:

    https://jeetsukumaran.github.io/DendroPy/downloading.html

Documentation
=============

Full documentation is available here:

    https://jeetsukumaran.github.io/DendroPy/

This includes:

    -   `A comprehensive "getting started" primer <https://jeetsukumaran.github.io/DendroPy/primer/index.html>`_ .
    -   `API documentation <https://jeetsukumaran.github.io/DendroPy/library/index.html>`_ .
    -   `Descriptions of data formats supported for reading/writing <https://jeetsukumaran.github.io/DendroPy/schemas/index.html>`_ .
    -   Guidance for `reporting issues <https://jeetsukumaran.github.io/DendroPy/index.html#bug-reports-and-other-issues>`_, `submitting feature requests <https://jeetsukumaran.github.io/DendroPy/index.html#feature-requests>`_, and `contributing to DendroPy <https://jeetsukumaran.github.io/DendroPy/developer.html>`_.

and more.

Citing
======

If you use any portion of DendroPy v5 in your research, please cite it as:

      Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a mature Python library for phylogenetic computing. Journal of Open Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943

For BibTex users:

.. code-block:: bibtex

  @misc{dendropy5,
    doi = {10.21105/joss.06943},
    url = {https://doi.org/10.21105/joss.06943},
    year = {2024},
  publisher = {The Open Journal},
    volume = {9},
    number = {101},
    pages = {6943},
    author = {Matthew Andres Moreno and Mark T. Holder and Jeet Sukumaran},
    title = {DendroPy 5: a mature Python library for phylogenetic computing},
    journal = {Journal of Open Source Software}
  }

Earlier DendroPy versions can be cited as:

      Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library for phylogenetic computing. Bioinformatics 26: 1569-1571. https://doi.org/10.1093/bioinformatics/btq228

Consider also leaving a `star on GitHub <https://github.com/jeetsukumaran/DendroPy/stargazers>`_!

License and Warranty
====================

Please see the file "LICENSE.rst" for details.

Developers
==========

- `Jeet Sukumaran <https://sukumaranlab.org/people/>`_
- `Mark Holder <https://phylo.bio.ku.edu/content/mark-t-holder>`_
- `Matthew Andres Moreno <https://mmore500.com/>`_
