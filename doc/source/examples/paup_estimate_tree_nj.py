#! /usr/bin/env python

import dendropy
from dendropy.interop import paup

data = dendropy.DnaCharacterMatrix.get_from_path("pythonidae.nex", "nexus")
tree = paup.estimate_tree(data,
        tree_est_criterion='nj')
print tree.as_string("newick")
