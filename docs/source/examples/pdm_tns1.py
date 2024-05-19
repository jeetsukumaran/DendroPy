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

tns = dendropy.TaxonNamespace()
mle_tree = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus",
        taxon_namespace=tns)
pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=open("pythonidae.mle.weighted.pdm.csv"),
        delimiter=",",
        taxon_namespace=tns)
