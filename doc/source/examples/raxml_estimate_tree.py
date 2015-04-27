#! /usr/bin/env python

import dendropy
from dendropy.interop import raxml

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
rx = raxml.RaxmlRunner()
tree = rx.estimate_tree(
        char_matrix=data,
        raxml_args=["--no-bfgs"])
print(tree.as_string(schema="newick"))

