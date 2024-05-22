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
        tree_est_criterion='likelihood',
        num_subst=2,
        unequal_base_freqs=True,
        gamma_rates=False,
        prop_invar=False)
print(tree.as_string(schema="newick"))
