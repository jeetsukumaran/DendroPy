import copy
import dendropy

# original list
s1 = "(A,(B,C));"
tree1 = dendropy.Tree.get(
        data=s1,
        schema="newick")

# Full deep copy by calling copy.deepcopy()
# I.e. Everything cloned including Taxon and TaxonNamespace instances
tree2 = copy.deepcopy(tree1)

# *different* tree instances
for nd1, nd2 in zip(tree1, tree2):
    assert nd1 is not nd2

# Note: TaxonNamespace is also different
assert tree2.taxon_namespace is not tree1.taxon_namespace
for tx1 in tree1.taxon_namespace:
    assert tx1 not in tree2.taxon_namespace
for tx2 in tree2.taxon_namespace:
    assert tx2 not in tree1.taxon_namespace
