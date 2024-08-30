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
  - name: Mark T. Holder
    orcid: 0000-0001-5575-0536
    affiliation: "4, 5"
  - name: Jeet Sukumaran
    orcid: 0000-0002-9222-9608
    affiliation: "6"
affiliations:
  - name: Department of Ecology and Evolutionary Biology, University of Michigan, Ann Arbor, MI, USA
    index: 1
  - name: Center for the Study of Complex Systems, University of Michigan, Ann Arbor, MI, USA
    index: 2
  - name: Michigan Institute for Data and AI in Society, University of Michigan, Ann Arbor, MI, USA
    index: 3
  - name: Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS, USA
    index: 4
  - name: Biodiversity Institute, University of Kansas, Lawrence, KS, USA
    index: 5
  - name: Department of Biology, San Diego State University, San Diego, CA, USA
    index: 6
date: 19 May 2024
bibliography: paper.bib
---

# Summary

Modern bioinformatics has unlocked remarkable insight into the composition, structure, and history of the natural world around us.
Arguably, the central pillar of bioinformatics is phylogenetics --- the study of hereditary relatedness among organisms.
Insights from phylogenetic analysis have touched nearly every corner of biology.
Examples range across natural history [@title2024macroevolutionary], population genetics and phylogeography [@knowles2002phylogeography], conservation biology [@faithConservationEvaluationPhylogenetic1992], public health [@giardina2017inference; @voznica2022deep], medicine [@lewinsohnStatedependentEvolutionaryModels2023; @kim2006discovery], *in vivo* and *in silico* experimental evolution [@rozen2005longterm,@moreno2023toward; @lenski2003evolutionary], application-oriented evolutionary algorithms [@lalejini2024phylogeny; @hernandez2022can; @shahbandegan2022untangling], and beyond.

High-throughput genetic and phenotypic data has realized groundbreaking results, in large part, through conjunction with open-source software used to process and analyze it.
Indeed, the preceding decades have ushered in a flourishing ecosystem of bioinformatics software applications and libraries.
Over the course of its nearly fifteen-year history, the DendroPy library for phylogenetic computation in Python has established a generalist niche in serving the bioinformatics community [@sukumaran2010dendropy].
Here, we report on the recent major release of the library, DendroPy version 5.
The software release represents a major milestone in transitioning the library to a sustainable long-term development and maintenance trajectory.
As such, this work positions DendroPy to continue fulfilling a key supporting role in phyloinformatics infrastructure.

# Statement of Need

DendroPy operates within a rich ecosystem of packages, frameworks, toolkits, and software projects supporting bioinformatics and phylogenetics research.
The broader software landscape largely divides into the following major categories,

1. High-performance specialized tools for inference (e.g., *BEAST2*, *RAxML*, *MrBayes*, *PAUP*, etc.) [@bouckaert2014beast; @stamatakis2014raxml; @ronquist2012mrbayes; @wilgenbusch2003inferring];
2. Python phylogenetics libraries that provide rich tree-centric data models and operations, such as
    - *ETE*, known in particular for powerful phylogeny visualization capabilities [@huertacepas2016ete],
    - *Scikit-bio* and tskit [@rideout2024scikitbio; @kelleher2018efficient],
    - *TreeSwift* and *SuchTree*, which provide lightweight, high-performance tree representations [@moshiri2020treeswift; @neches2018suchtree], and
    - *hstrat* and *Phylotrack*, which specialize in collecting phylogenies from agent-based evolutionary simulation [@dolson2024phylotrack; @moreno2022hstrat];
3. Python phylogenetics libraries with genome/gene-centric data models and operations (e.g., *PyCogent*/*Cogent3*, *BioPython*, etc.) [@knight2007pycogent; @cock2009biopython]; and
4. Numerous R phylogenetics packages, which are often highly specialized but generally interoperate via `ape.phylo` data structures [@paradis2019ape].

DendroPy falls largely within the second camp above.
It is notable in providing a broad portfolio of evolutionary models, but also fields population genetics and sequence evolution utilities.
DendroPy is also notable for its comprehensive, systematic documentation and rich, user-extensible tree representation.
The library's use cases range across serving as a stand-alone library for phylogenetics, a component of more complex multi-library phyloinformatics pipelines, or as an interstitial "glue" that assembles and drives such pipelines.

# Features

Key features of DendroPy are:

- rich object-oriented representations for manipulation of phylogenetic trees and character matrices;
- efficient, bit-level representation of nodes' leaf bipartitions;
- loading and saving popular phylogenetic data formats, including NEXUS, Newick, NeXML, Phylip, and FASTA [@maddison1997nexus; @olsen1990newick; @vos2012nexml; @felsenstein1981evolutionary; @lipman1985rapid];
- simulation of phylogenetic trees under a range of models, including coalescent models, birth-death models, and population genetics simulations of gene trees; and
- application scripts for performing data conversion, collating taxon labels from multiple trees, and tree posterior distribution summarization.

Significant improvements have been made since DendroPy's original release [@sukumaran2010dendropy], including performance enhancements in saving and loading trees, support for distance-based tree construction, and addition of new phylogeny statistics and speciation models.

# Maintenance

The primary focus of DenroPy's version 5 release is to support sustainable long-term maintenance by reducing the effort needed for new releases.
We hope this will result in a regular release schedule (incorporating timely patches for reported issues), development of new features, and incorporation of user contributions.

The version 5 release reflects substantial investment in adopting modern software development best practices.
In version 5, DendroPy has officially dropped support for Python 2.7, as well as Python 3.X versions that have reached end-of-life.
Focusing support on Python 3.6 and higher simplifies cross-environment testing and allows future development to leverage new language features.
In addition, we have established comprehensive continuous integration (CI) infrastructure via GitHub Actions, comprising

- code linting with [Ruff](https://pypi.org/project/ruff/);
- deploying up-to-date documentation via GitHub pages;[^1]
- unit tests, largely organized within the `unittest` framework;
- new smoke tests using [pytest](https://pypi.org/project/pytest/);
- code coverage reporting via the Codecov service; and
- automatic deployment of tagged versions to PyPI.

[^1]: Documentation is hosted at [https://jeetsukumaran.github.io/DendroPy](https://jeetsukumaran.github.io/DendroPy).

Other behind-the-scenes activity in preparing this release includes repair of library components flagged by the new tooling, triage of user bug reports, applying issue tags to manage open tracker items, establishing a code of conduct, and creating issue templates to increase the quality of future bug reports and feature requests.
Altogether, these improvements serve as a foundation for future work maintaining and extending DendroPy in a manner that is reliable, stable, and responsive to user needs.

# Impact

Over its nearly 15-year history, DendroPy’s versatility and stability have driven adoption as a core dependency of many phylogenetics pipelines and bioinformatics software libraries.
Currently, [85 projects](https://perma.cc/P865-JHXW) on PyPI list DendroPy as a direct dependency.
Notable projects using DendroPy include:

- PASTA, which performs multiple sequence alignment [@mirarab2014pasta];
- Physcraper, which automates curation of gene trees [@sanchezreyes2021physcraper];
- Propinquity, the supertree pipeline [@redelings2017supertree] of the Open Tree of Life project;
- DELINEATE, software for analyses discerning true speciation from population lineages [@sukumaran2021incorporating];
- Archipelago, which models spatially explicit biographical phylogenesis [@sukumaran2015machine];
- Espalier, a utility for constructing maximum agreement forests [@rasmussen2023espalier]; and
- MetaPhlAn, which extracts information about microbial community composition from metagenomic shotgun sequencing data [@blancomiguez2023extending].

During this time, DendroPy has also directly helped enable numerous end-user phylogenetics projects.
Notable examples include work on the early natural history of birds [@jarvis2014wholegenome], the molecular evolution of the Zika virus [@faye2014molecular], and early human migration within the Americas [@garcaortiz2021genomic].
As of May 2024, Google Scholar counts 1,654 works referencing DendroPy [@sukumaran2010dendropy].

# Acknowledgements

Thank you to University of Michigan Undergraduate Research Opportunity Program participant Connor Yang for his contributions in increasing test coverage, and to our open-source community for bug reports, feature suggestions, and patch contributions over the years.
This research is based upon work supported by:

- the Eric and Wendy Schmidt AI in Science Postdoctoral Fellowship, a Schmidt Sciences program (author MAM);
- the National Science Foundation grant NSF-DEB 1937725 "COLLABORATIVE RESEARCH: Phylogenomics, spatial phylogenetics and conservation prioritization in trapdoor spiders (and kin) of the California Floristic Province" (author JS); and
- the National Science Foundation grant NSF-DEB 1457776 "Collaborative Research Developing novel methods for estimating coevolutionary processes
using tapeworms and their shark and ray hosts" (author MH).

# References

<div id="refs"></div>

\pagebreak
\appendix
