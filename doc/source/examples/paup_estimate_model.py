#! /usr/bin/env python

import dendropy
from dendropy.interop import paup

data = dendropy.DnaCharacterMatrix.get_from_path("pythonidae.nex", "nexus")
tree = paup.estimate_tree(data,
        tree_est_criterion='nj')
model = paup.estimate_model(data,
        tree,
        num_states=2,
        unequal_base_freqs=True,
        gamma_rates=False,
        prop_invar=False)
for k, v in model:
    print k, v
