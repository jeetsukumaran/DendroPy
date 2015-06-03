#! /usr/bin/env python

import dendropy
import random

categories = {
    "A" : "N/A",
    "B" : "N/A",
    "C" : "N/A",
    "D" : "N/A",
    "E" : "N/A"
}
tree = dendropy.Tree.get(
        data="(A,(B,(C,(D,E))));",
        schema="newick")
for taxon in tree.taxon_namespace:
    taxon.category = categories[taxon.label]
    taxon.annotations.add_bound_attribute("category")
for node in tree.postorder_node_iter():
    node.pop_size = None
    node.annotations.add_bound_attribute("pop_size")
for node in tree.postorder_node_iter():
    node.pop_size = random.randint(100, 10000)
    if node.taxon is not None:
        if node.pop_size >= 8000:
            node.taxon.category = "large"
        elif node.pop_size >= 6000:
            node.taxon.category = "medium"
        elif node.pop_size >= 4000:
            node.taxon.category = "small"
        elif node.pop_size >= 2000:
            node.taxon.category = "tiny"
print tree.as_string(schema="nexml")

