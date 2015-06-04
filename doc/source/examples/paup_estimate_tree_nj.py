#! /usr/bin/env python

import dendropy
from dendropy.interop import paup

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
tree = paup.estimate_tree(data,
        tree_est_criterion='nj')
print tree.as_string(schema="newick")
