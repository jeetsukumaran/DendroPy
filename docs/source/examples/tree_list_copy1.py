import dendropy

# original list
s1 = "(A,(B,C));(B,(A,C));(C,(A,B));"
treelist1 = dendropy.TreeList.get(
        data=s1,
        schema="newick")

# shallow copy by calling Tree.clone(0)
treelist2 = treelist1.clone(depth=0)

# shallow copy by slicing
treelist3 = treelist1[:]

# same tree instances are shared
for t1, t2 in zip(treelist1, treelist2):
    assert t1 is t2
for t1, t2 in zip(treelist1, treelist3):
    assert t1 is t2

# note: (necessarily) sharing same TaxonNamespace
assert treelist2.taxon_namespace is treelist1.taxon_namespace
assert treelist3.taxon_namespace is treelist1.taxon_namespace

