import copy
import dendropy

# original list
s1 = "(A,(B,C));(B,(A,C));(C,(A,B));"
treelist1 = dendropy.TreeList.get(
        data=s1,
        schema="newick")

# Full deep copy by calling copy.deepcopy()
# I.e. Everything cloned including Taxon and TaxonNamespace instances
treelist2 = copy.deepcopy(treelist)

# *different* tree instances
for t1, t2 in zip(treelist1, treelist2):
    assert t1 is not t2

# Note: TaxonNamespace is also different
assert treelist2.taxon_namespace is not treelist1.taxon_namespace
for tx1 in treelist1.taxon_namespace:
    assert tx1 not in treelist2.taxon_namespace
for tx2 in treelist2.taxon_namespace:
    assert tx2 not in treelist1.taxon_namespace
