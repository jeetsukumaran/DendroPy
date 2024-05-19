# /usr/bin/env python

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

def generate(birth_rates, death_rates):
    assert len(birth_rates) == len(death_rates)
    tree = dendropy.Tree()
    for i, br in enumerate(birth_rates):
        tree = treesim.birth_death_tree(birth_rates[i],
                                   death_rates[i],
                                   max_time=random.randint(1,8),
                                   tree=tree,
                                   assign_taxa=False,
                                   repeat_until_success=True)
        print(tree.as_string(schema='newick'))
    tree.randomly_assign_taxa(create_required_taxa=True)
    return tree

tree = generate([0.1, 0.6, 0.1], [0.1, 0.6, 0.1])
print(tree.as_string(schema='newick'))
