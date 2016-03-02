#! /usr/bin/env python

import dendropy

label_transform_fn = lambda x: x.replace("_", " ")
pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=open("pythonidae.mle.weighted.pdm.csv"),
        delimiter=",",
        label_transform_fn=label_transform_fn)
nj_tree = pdm.nj_tree()
print(nj_tree.as_string("newick"))
