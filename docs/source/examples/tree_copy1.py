import dendropy

# original list
s1 = "(A,(B,C));"
tree1 = dendropy.Tree.get(
        data=s1,
        schema="newick")

# taxon namespace-scoped deep copy by calling Tree.clone(1)
# I.e. Everything cloned, but with Taxon and TaxonNamespace references shared
tree2 = tree1.clone(depth=1)

# taxon namespace-scoped deep copy by copy-construction
# I.e. Everything cloned, but with Taxon and TaxonNamespace references shared
tree3 = dendropy.Tree(tree1)

# *different* tree instances, with different nodes and edges
for nd1, nd2, nd3 in zip(tree1, tree2, tree3):
    assert nd1 is not nd2
    assert nd1 is not nd3
    assert nd2 is not nd3

# Note: TaxonNamespace is still shared
# I.e. Everything cloned, but with Taxon and TaxonNamespace references shared
assert tree2.taxon_namespace is tree1.taxon_namespace
assert tree3.taxon_namespace is tree1.taxon_namespace
