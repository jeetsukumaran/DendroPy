#! /usr/bin/env python

import dendropy
from dendropy.interop import phyml

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
result = phyml.run_phyml(
    phyml_path='phyml',
    char_matrix=data,
    data_type="nt",
    subst_model="HKY85",
    gamma_cats=1)
print(result.best_tree.as_string(schema="newick"))
