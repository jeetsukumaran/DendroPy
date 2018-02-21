#! /usr/bin/env python
# -*- coding: utf-8 -*-

import dendropy

tns = dendropy.TaxonNamespace()
mle_tree = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus",
        taxon_namespace=tns)
pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=open("pythonidae.mle.weighted.pdm.csv"),
        delimiter=",",
        taxon_namespace=tns)
