<div align="center">

[![DendroPy wordmark](https://raw.githubusercontent.com/jeetsukumaran/DendroPy/DendroPy4/doc/source/_static/dendropy_logo.png)](https://github.com/jeetsukumaran/DendroPy)

| **automated tests:** | **package version:** | **documentation:** | **test coverage:** |
|----------------------|----------------------|--------------------|--------------------|
| [![Continuous Integration build](https://github.com/jeetsukumaran/DendroPy/actions/workflows/ci.yaml/badge.svg)](https://github.com/jeetsukumaran/DendroPy/actions/workflows/ci.yaml) | [![PyPI version](https://img.shields.io/pypi/v/DendroPy.svg)](https://pypi.org/project/DendroPy/) | [![Documentation Status](https://github.com/jeetsukumaran/DendroPy/actions/workflows/pages/pages-build-deployment/badge.svg)](https://jeetsukumaran.github.io/DendroPy/) | [![codecov coverage](https://codecov.io/gh/jeetsukumaran/DendroPy/graph/badge.svg?token=JwMfFOpBBD)](https://codecov.io/gh/jeetsukumaran/DendroPy) |

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

<https://jeetsukumaran.github.io/DendroPy/>

DendroPy is also hosted in the official Python Packaging Index (PyPI):

<http://pypi.org/project/DendroPy/>

# Requirements and Installation

The current version of DendroPy requires Python 3.

You can install DendroPy by running::

    $ python -m pip install dendropy

For Conda users, DendroPy can be installed from the [conda-forge](https://conda-forge.org/) channel:

    $ conda install -c conda-forge dendropy

More information is available here:

<https://jeetsukumaran.github.io/DendroPy/downloading.html>

# Documentation

Full documentation is available here:

<https://jeetsukumaran.github.io/DendroPy/>

This includes:

-   [A comprehensive "getting started" primer](https://jeetsukumaran.github.io/DendroPy/primer/index.html).
-   [API documentation](https://jeetsukumaran.github.io/DendroPy/library/index.html).
-   [Descriptions of data formats supported for reading/writing](https://jeetsukumaran.github.io/DendroPy/schemas/index.html).
-   Guidance for [reporting issues](https://jeetsukumaran.github.io/DendroPy/index.html#bug-reports-and-other-issues), [submitting feature requests](https://jeetsukumaran.github.io/DendroPy/index.html#feature-requests), and [contributing to DendroPy](https://jeetsukumaran.github.io/DendroPy/developer.html).

and more.

# Citing

If you use any portion of DendroPy v5 in your research, please cite it as:
> Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a mature Python library for phylogenetic computing. Journal of Open Source Software, 9(101), 6943, [https://doi.org/10.21105/joss.06943](https://doi.org/10.21105/joss.06943)

For BibTex users:
```bibtex
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
````

Earlier DendroPy versions can be cited as:
> Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library for phylogenetic computing. Bioinformatics 26: 1569-1571. https://doi.org/10.1093/bioinformatics/btq228

Consider also leaving a [star on GitHub](https://github.com/jeetsukumaran/DendroPy/stargazers)!

# License and Warranty

Please see the file "LICENSE.rst" for details.

# Developers

- [Jeet Sukumaran](https://sukumaranlab.org/people/)
- [Mark Holder](https://phylo.bio.ku.edu/content/mark-t-holder)
- [Matthew Andres Moreno](https://mmore500.com/)

```
