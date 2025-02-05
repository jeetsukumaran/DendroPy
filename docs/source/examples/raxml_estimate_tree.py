#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import dendropy
from dendropy.interop import raxml

warnings.warn("This example requires Raxml to be installed to run.")


data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
rx = raxml.RaxmlRunner()
tree = rx.estimate_tree(
        char_matrix=data,
        raxml_args=["--no-bfgs"])
print(tree.as_string(schema="newick"))

