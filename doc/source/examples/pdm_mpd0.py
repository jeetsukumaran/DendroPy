#! /usr/bin/env python

import dendropy

tree = dendropy.Tree.get(
    path="pythonidae.mle.nex",
    schema="nexus")
pdm = tree.phylogenetic_distance_matrix()

# MPD of entire tree
print(pdm.mean_pairwise_distance())

# MNTD of entire tree
print(pdm.mean_nearest_taxon_distance())

# Statistics of a "community" consisting of first
# 8 taxa
community_taxa = set(tree.taxon_namespace[:8])
filter_fn = lambda taxon : taxon in community_taxa
print(pdm.mean_pairwise_distance(filter_fn=filter_fn))
print(pdm.mean_nearest_taxon_distance(filter_fn=filter_fn))

# Statistics of a "community" consisting of
# species not in the Python genus
filter_fn = lambda taxon : not taxon.label.startswith("Python")
print(pdm.mean_pairwise_distance(filter_fn=filter_fn))
print(pdm.mean_nearest_taxon_distance(filter_fn=filter_fn))
