#! /usr/bin/env python

import dendropy
from dendropy.calculate import treecompare

tree0 = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus")

morelia_taxa = set(taxon for taxon in tree0.taxon_namespace
        if taxon.label.startswith("Morelia"))
morelia_labels = set([t.label for t in morelia_taxa])
non_morelia_taxa = set(taxon for taxon in tree0.taxon_namespace
        if not taxon.label.startswith("Morelia"))
non_morelia_labels = set([t.label for t in non_morelia_taxa])

tree1 = tree0.extract_tree_with_taxa(taxa=morelia_taxa)
tree2 = tree0.extract_tree_with_taxa_labels(labels=morelia_labels)
tree3 = tree0.extract_tree_without_taxa(taxa=non_morelia_taxa)
tree4 = tree0.extract_tree_without_taxa_labels(labels=non_morelia_labels)

print tree1.as_string("newick")
print tree2.as_string("newick")
print tree3.as_string("newick")
print tree4.as_string("newick")

assert treecompare.weighted_robinson_foulds_distance(tree1, tree2) == 0.0
assert treecompare.weighted_robinson_foulds_distance(tree2, tree3) == 0.0
assert treecompare.weighted_robinson_foulds_distance(tree3, tree4) == 0.0
