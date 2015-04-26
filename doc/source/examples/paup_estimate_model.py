#! /usr/bin/env python

import dendropy
from dendropy.interop import paup

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
tree = paup.estimate_tree(data,
        tree_est_criterion='nj')
est_tree, est_model = paup.estimate_model(data,
        tree,
        num_states=2,
        unequal_base_freqs=True,
        gamma_rates=False,
        prop_invar=False)
for k, v in est_model.items():
    print("{}: {}".format(k, v))
