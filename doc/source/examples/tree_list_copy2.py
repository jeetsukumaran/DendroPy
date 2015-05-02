import copy
import dendropy

# original list
s1 = "(A,(B,C));(B,(A,C));(C,(A,B));"
treelist1 = dendropy.TreeList.get(
        data=s1,
        schema="newick")

# deep copy by copy-construction
treelist2 = dendropy.TreeList(treelist1)

# *different* tree instances
for t1, t2 in zip(treelist1, treelist2):
    assert t1 is not t2

# note: TaxonNamespace is still shared
assert treelist2.taxon_namespace is treelist1.taxon_namespace
