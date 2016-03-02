#! /usr/bin/env python

import dendropy

pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=open("pythonidae.mle.weighted.pdm.csv"),
        delimiter=",")
upgma_tree = pdm.upgma_tree()
print(upgma_tree.as_string("newick"))
