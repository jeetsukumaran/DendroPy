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
@article {MEN:MEN2925,
author = {BROWN, JEREMY M. and SAVIDGE, KEVIN and McTAVISH, EMILY JANE B.},
title = {DIM SUM: demography and individual migration simulated using a Markov chain},
journal = {Molecular Ecology Resources},
volume = {11},
number = {2},
publisher = {Blackwell Publishing Ltd},
issn = {1755-0998},
url = {http://dx.doi.org/10.1111/j.1755-0998.2010.02925.x},
doi = {10.1111/j.1755-0998.2010.02925.x},
pages = {358--363},
keywords = {dispersal, Markov chain, migration, phylogeography, simulation},
year = {2011},
}
"""

dataset = dendropy.DataSet.get(
        data="(A,(B,(C,(D,E))));",
        schema="newick")
dataset.annotations.add_citation(citation,
        store_as="prism")
dataset.annotations.add_citation(citation,
        store_as="dublin")
print(dataset.as_string(schema="nexml"))
