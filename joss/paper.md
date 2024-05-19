---
title: 'DendroPy 5: a mature Python library for phylogenetic computing'
tags:
  - Python
  - phylogenetics
  - phylogeny
  - phylogenies
  - phylogeography
  - biology
  - evolution
  - evolutionary biology
  - systematics
  - coalescent
  - population
  - genetics
  - phyloinformatics
  - bioinformatics
authors:
  - name: Matthew Andres Moreno
    orcid: 0000-0003-4726-4479
    affiliation: "1, 2, 3"
  - name: Jeet Sukumaran
    orcid: 0000-0002-9222-9608
    affiliation: "4"
  - name: Mark T. Holder
    orcid: 0000-0001-5575-0536
    affiliation: "5, 6"
affiliations:
  - name: Department of Ecology and Evolutionary Biology, University of Michigan, Ann Arbor, MI, USA
   index: 1
  - name: Center for the Study of Complex Systems, University of Michigan, Ann Arbor, MI, USA
   index: 2
  - name: Michigan Institute for Data Science, University of Michigan, Ann Arbor, MI, USA
   index: 3
  - name: Department of Biology, San Diego State University, San Diego, CA, USA
   index: 4
  - name: Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS, USA
   index: 5
  - name: Biodiversity Institute, University of Kansas, Lawrence, KS, USA
   index: 6
date: 19 May 2024
bibliography: paper.bib
---

<!-- @MAM title could also be "mainstay" instead of "mature" -->

# Summary

Phylogenies not only tell the history of life, but are phylogenetic also predictive and can be applied to do things: conservation, public health, medicine, and evolutionary algorithms.


# Statement of Need

There are a broad range of packages, frameworks, toolkits, and software projects that are immediately relevant to this domain.
The broader evolutionary biology analysis software landscape divides into the following major categories:

1. High-performance specialized tools for inference (e.g., BEAST2, RAxML, MrBayes, PAUP, etc.) [@bouckaert2014beast2,@stamatakis2014raxml,@ronquist2012mrbayes,@wilgenbusch2003inferring];
2. Python phylogenetics libraries providing rich tree-centric data models and operations:
  - ETE [@huertacepas2016ete]
  - Scikit-bio [@rideout2024scikitbio]
  - tskit [@kelleher2018efficient]
  - TreeSwift and SuchTree, which feature lightweight,  high-performance tree representations [@moshiri2020treeswift,@neches2018suchtree]
  - hstrat and Phylotrack, which are specialized for working with simulation-generated phylogenies [@dolson2024phylotrack,@moreno2022hstrat].
3. Python phylogenetics libraries with genome/gene-centric data models and operations:
  - PyCogent (now succeeded by Cogent3)[@knight2007pycogent],
  - BioPython [@cock2009biopython]
4/ In R, there are a large number of often highly-specialized tools that generally interoperate via `ape.phylo` data structures [@paradis2019ape].

DendroPy provides a particularly broad portfolio of evolutionary models with comprehensive, systematic documentation
DendroPy can thus function as a stand-alone library for phylogenetics, a component of more complex multi-library phyloinformatics pipelines, or as a scripting "glue" that assembles and drives such pipelines [@sukumaran2010dendropy].

# Features

- rich object-oriented representations for manipulation of phylogenetic trees and character matrices,
- reading and writing of phylogenetic data in a range of formats, such as NEXUS, Newick, NeXML, Phylip, FASTA, etc. [@maddison1997nexus,@olsen1990newick,@vos2012nexml,@felsenstein1981evolutionary,@lipman1985rapid]
- simulation of phylogenetic trees under a range of models, including coalescent models, birth-death models, and population genetics simulations of gene trees.
- application scripts for performing some useful phylogenetic operations, such as data conversion and tree posterior distribution summarization.

# Maintenance

Sustainable long-term maintenance is the major focus of the version 5 release.
To this end, we have implemented continuous integration through GitHub Actions comprising:
- code linting with [Ruff](https://pypi.org/project/ruff/)
- build and deployment of documentation via GitHub pages, hosted at <jeetsukumaran.github.io/DendroPy>,
- unit tests, largely using the unittest framework,
- new smoke tests using [pytest](https://pypi.org/project/pytest/),
- code coverage reporting via the codecov service,
- automatic deployment to PyPI upon creation of a new version.

Other key activities in preparing this release included triage of bug reports, creating an issue organization system to manage open issues, and creating issue templates to increase the quality of future bug reports and feature requests.

The library has reached a major milestone in dropping support for end-of-life (EOL) Python versions, including Python 2.7.
This brings DendroPy in line with best practice, makes it easier to test and ensure support, and allows us to take advantage of new language features.


# Impact

Over its nearly 15-year history, DendroPy’s stability has led to adoption as a core dependency of other software and pipelines.
There are currently [85 projects](https://perma.cc/P865-JHXW) on PyPI listed as depending on DendroPy.

- PASTA [@mirarab2014pasta]
- Physcraper [@sanchezreyes2021physcraper]
- Espalier [@rasmussen2023espalier],
- MetaPhlAn [@blancomiguez2023extending]

# Acknowledgements

Thank you to Connor Yang for his contributions in increasing test coverage, and to the community for bug reports, feature suggestions, and patch contributions over the years.
This research is based upon work supported by the Eric and Wendy Schmidt AI in Science Postdoctoral Fellowship, a Schmidt Futures program.

# References

<div id="refs"></div>

\pagebreak
\appendix
