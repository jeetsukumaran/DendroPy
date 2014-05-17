#! /usr/bin/env python

import dendropy
from dendropy.interop import raxml

data = dendropy.DnaCharacterMatrix.get_from_path("pythonidae.nex", "nexus")
rx = raxml.RaxmlRunner()
tree = rx.estimate_tree(data)
