#! /usr/bin/env python

import dendropy
from dendropy.interop import phyml

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
tree = dendropy.Tree.get(
    path="pythonidae.mle.nex",
    schema="nexus")
result = phyml.run_phyml(
    char_matrix=data,
    data_type="nt",
    subst_model="HKY85",
    gamma_cats=1,
    starting_tree=tree,
    optimization="lr"  # optimize branch length and rate parameters only
    )
print(result.stats_text)