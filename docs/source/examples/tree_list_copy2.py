import dendropy

# original list
s1 = "(A,(B,C));(B,(A,C));(C,(A,B));"
treelist1 = dendropy.TreeList.get(
        data=s1,
        schema="newick")

# taxon namespace-scoped deep copy by calling Tree.clone(1)
# I.e. Everything cloned, but with Taxon and TaxonNamespace references shared
treelist2 = treelist1.clone(depth=1)

# taxon namespace-scoped deep copy by copy-construction
# I.e. Everything cloned, but with Taxon and TaxonNamespace references shared
treelist3 = dendropy.TreeList(treelist1)

# *different* tree instances
for t1, t2, t3 in zip(treelist1, treelist2, treelist3):
    assert t1 is not t2
    assert t1 is not t3
    assert t2 is not t3

# Note: TaxonNamespace is still shared
# I.e. Everything cloned, but with Taxon and TaxonNamespace references shared
assert treelist2.taxon_namespace is treelist1.taxon_namespace
assert treelist3.taxon_namespace is treelist1.taxon_namespace
