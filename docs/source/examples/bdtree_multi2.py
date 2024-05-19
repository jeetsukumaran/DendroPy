#! /usr/bin/env python
# -*- coding: utf-8 -*-

import random
import warnings
import dendropy
from dendropy.simulate import treesim

warnings.warn(
    "This example is known to be broken! "
    "It will be fixed or removed in the future. "
    "See https://github.com/jeetsukumaran/DendroPy/issues/160 for details. "
    "Patch contributions are welcome.",
)

def generate(mean, sd, num_periods):
    tree = dendropy.Tree()
    for i in range(num_periods):
        tree = treesim.birth_death_tree(birth_rate=random.gauss(mean, sd),
                                   death_rate=random.gauss(mean, sd),
                                   max_time=random.randint(1,5),
                                   tree=tree,
                                   assign_taxa=False,
                                   repeat_until_success=True)
    tree.randomly_assign_taxa(create_required_taxa=True)
    return tree

tree = generate(0.1, 0.01, 100)
print(tree.as_string(schema='newick'))
