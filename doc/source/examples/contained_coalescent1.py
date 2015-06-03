#! /usr/bin/env python

import dendropy
from dendropy.simulate import treesim


sp_tree_str = """\
[&R] (A:10,(B:6,(C:4,(D:2,E:2):2):2):4);
"""

sp_tree = dendropy.Tree.get(data=sp_tree_str, schema="newick")
gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=sp_tree.taxon_namespace,
        num_contained=3)
gene_tree = treesim.contained_coalescent_tree(containing_tree=sp_tree,
    gene_to_containing_taxon_map=gene_to_species_map)
print(gene_tree.as_string(schema='newick'))
print(gene_tree.as_ascii_plot())

