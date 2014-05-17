#! /usr/bin/env python

import random
import dendropy
from dendropy import treesim

def generate(mean, sd, num_periods):
    tree = dendropy.Tree()
    for i in range(num_periods):
        tree = treesim.birth_death(birth_rate=random.gauss(mean, sd),
                                   death_rate=random.gauss(mean, sd),
                                   max_time=random.randint(1,5),
                                   tree=tree,
                                   assign_taxa=False,
                                   repeat_until_success=True)
    tree.randomly_assign_taxa(create_required_taxa=True)
    return tree

tree = generate(0.1, 0.01, 100)
print(tree.as_string('newick'))
