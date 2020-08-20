#! /usr/bin/env python
# -*- coding: utf-8 -*-

import dendropy
from dendropy.interop import paup

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
tree = paup.estimate_tree(data,
        tree_est_criterion='likelihood',
        num_subst=2,
        unequal_base_freqs=True,
        gamma_rates=False,
        prop_invar=False)
print(tree.as_string(schema="newick"))
