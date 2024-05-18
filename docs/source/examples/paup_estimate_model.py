#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import dendropy
from dendropy.interop import paup

warnings.warn(
    "This example is known to be broken! "
    "It will be fixed or removed in the future. "
    "See https://github.com/jeetsukumaran/DendroPy/issues/160 for details. ",
    "Patch contributions are welcome.",
)

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
tree = paup.estimate_tree(data,
        tree_est_criterion='nj')
est_tree, est_model = paup.estimate_model(data,
        tree,
        num_subst=2,
        unequal_base_freqs=True,
        gamma_rates=False,
        prop_invar=False)
for k, v in est_model.items():
    print("{}: {}".format(k, v))
