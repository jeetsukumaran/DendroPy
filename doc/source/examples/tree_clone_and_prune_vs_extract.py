#! /usr/bin/env python

import dendropy
import timeit

tree0 = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus")
taxa_to_prune = set(taxon for taxon in tree0.taxon_namespace
        if taxon.label.startswith("Morelia"))

def f1():
    tree1 = dendropy.Tree(tree0)
    tree1.prune_taxa(taxa_to_prune)

def f2():
    tree1 = tree0.extract_tree(
            node_filter_fn=lambda taxon: taxon not in taxa_to_prune)

t1 = timeit.Timer(f1)
print(min(t1.repeat(10,10))/10)

t2 = timeit.Timer(f2)
print(min(t2.repeat(10,10))/10)
