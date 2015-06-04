#! /usr/bin/env python

import dendropy
from dendropy.simulate import treesim
from dendropy.model import reconcile

# simulation parameters and output
num_reps = 10

# population tree descriptions
stepwise_tree_str = "[&R](A:120000,(B:80000,(C:40000,D:40000):40000):40000):100000;"
frag_tree_str = "[&R](A:120000,B:120000,C:120000,D:120000):100000;"

# taxa and trees
containing_taxa = dendropy.TaxonNamespace()
stepwise_tree = dendropy.Tree.get(
        data=stepwise_tree_str,
        schema="newick",
        taxon_namespace=containing_taxa)
frag_tree = dendropy.Tree.get(
        data=frag_tree_str,
        schema="newick",
        taxon_namespace=containing_taxa)

# taxon set association
genes_to_species = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=containing_taxa,
        num_contained=8)

# convert to containing tree
stepwise_tree = reconcile.ContainingTree(stepwise_tree,
            contained_taxon_namespace=genes_to_species.domain_taxon_namespace,
            contained_to_containing_taxon_map=genes_to_species)
frag_tree = reconcile.ContainingTree(frag_tree,
            contained_taxon_namespace=genes_to_species.domain_taxon_namespace,
            contained_to_containing_taxon_map=genes_to_species)

# for each rep
for rep in range(num_reps):
    gene_tree1 = treesim.contained_coalescent_tree(containing_tree=stepwise_tree,
        gene_to_containing_taxon_map=genes_to_species,
        default_pop_size=40000)
    stepwise_tree.embed_tree(gene_tree1)
    gene_tree2 = treesim.contained_coalescent_tree(containing_tree=frag_tree,
        gene_to_containing_taxon_map=genes_to_species,
        default_pop_size=40000)
    frag_tree.embed_tree(gene_tree2)

# write results

# returns dictionary with contained trees as keys
# and number of deep coalescences as values
stepwise_deep_coals = stepwise_tree.deep_coalescences()
stepwise_out = open("stepwise.txt", "w")
for tree in stepwise_deep_coals:
    stepwise_out.write("%d\n" % stepwise_deep_coals[tree])

# returns dictionary with contained trees as keys
# and number of deep coalescences as values
frag_deep_coals = frag_tree.deep_coalescences()
frag_out = open("frag.txt", "w")
for tree in frag_deep_coals:
    frag_out.write("%d\n" % frag_deep_coals[tree])








