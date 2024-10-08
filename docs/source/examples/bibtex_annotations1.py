#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import dendropy

warnings.warn(
    "This example is known to be broken! "
    "It will be fixed or removed in the future. "
    "See https://github.com/jeetsukumaran/DendroPy/issues/160 for details. "
    "Patch contributions are welcome.",
)

citation = """\
@article{HeathHH2012,
	Author = {Tracy A. Heath and Mark T. Holder and John P. Huelsenbeck},
	Doi = {10.1093/molbev/msr255},
	Journal = {Molecular Biology and Evolution},
	Number = {3},
	Pages = {939-955},
	Title = {A {Dirichlet} Process Prior for Estimating Lineage-Specific Substitution Rates.},
	Url = {http://mbe.oxfordjournals.org/content/early/2011/11/04/molbev.msr255.abstract},
	Volume = {29},
	Year = {2012}
	}
"""


dataset = dendropy.DataSet.get(
        data="(A,(B,(C,(D,E))));",
        schema="newick")
dataset.annotations.add_citation(citation)
print(dataset.as_string(schema="nexml"))
