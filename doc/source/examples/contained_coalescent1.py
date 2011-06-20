#! /usr/bin/env python

import dendropy
from dendropy import treesim


sp_tree_str = """\
[&R] (A:10,(B:6,(C:4,(D:2,E:2):2):2):4)
"""

sp_tree = dendropy.Tree.get_from_string(sp_tree_str, "newick")
gene_to_species_map = dendropy.TaxonSetMapping.create_contained_taxon_mapping(
        containing_taxon_set=sp_tree.taxon_set,
        num_contained=3)
gene_tree = treesim.contained_coalescent(containing_tree=sp_tree,
    gene_to_containing_taxon_map=gene_to_species_map)
print(gene_tree.as_string('newick'))
print(gene_tree.as_ascii_plot())

