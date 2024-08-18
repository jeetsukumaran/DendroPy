.. raw:: html

   <div align="center">

.. image:: https://raw.githubusercontent.com/jeetsukumaran/DendroPy/DendroPy4/doc/source/_static/dendropy_logo.png
   :align: right
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

-----


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

The current version of DendroPy requires Python 3:

You can install DendroPy by running::

    $ python -m pip install dendropy

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

      Moreno, M. A., Sukumaran, J., and M. T. Holder. 2024. DendroPy 5: a mature Python library for phylogenetic computing. arXiv preprint arXiv:2405.14120. https://doi.org/10.48550/arXiv.2405.14120

For BibTex users:

.. code-block:: bibtex

      @misc{dendropy5,
         title = {DendroPy 5: a mature Python library for phylogenetic computing},
         author = {Moreno,  Matthew Andres and Sukumaran,  Jeet and Holder,  Mark T.},
         year = {2024},
         keywords = {Populations and Evolution (q-bio.PE),  FOS: Biological sciences,  FOS: Biological sciences},
         publisher = {arXiv},
         doi = {10.48550/ARXIV.2405.14120},
         url = {https://arxiv.org/abs/2405.14120},
         copyright = {arXiv.org perpetual, non-exclusive license}
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
