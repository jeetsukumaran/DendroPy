#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import dendropy
from dendropy.interop import paup

warnings.warn("This example requires paup to be installed to run.")

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
tree = paup.estimate_tree(data,
        tree_est_criterion='nj')
print(tree.as_string(schema="newick"))
