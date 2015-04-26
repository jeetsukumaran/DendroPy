#! /usr/bin/env python

import dendropy
from dendropy.interop import raxml

data = dendropy.DnaCharacterMatrix.get(
    path="pythonidae.nex",
    schema="nexus")
rx = raxml.RaxmlRunner()
tree = rx.estimate_tree(data)
