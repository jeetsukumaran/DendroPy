#! /usr/bin/env python

import dendropy
from dendropy import treesim

taxa = dendropy.TaxonNamespace(["z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8"])
tree = treesim.pure_kingman(
        taxon_namespace=taxa,
        pop_size=10000)
print(tree.as_string(schema="newick"))
print(tree.as_ascii_plot())
