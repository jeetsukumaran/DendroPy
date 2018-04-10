import dendropy
from dendropy.calculate import treecompare

trees = dendropy.TreeList.get(
        path="pythonidae.random.bd0301.tre",
        schema="nexus")

print(len(trees))

tree = dendropy.Tree.get(path="pythonidae.mle.nex", schema="nexus")

# As we did not specify a |TaxonNamespace| instance to use above, by default
# 'tree' will get its own, distinct |TaxonNamespace|
original_tree_taxon_namespace = tree.taxon_namespace
print(id(original_tree_taxon_namespace))
assert tree.taxon_namespace is not trees.taxon_namespace

# This operation adds the |Tree|, 'tree', to the |TreeList|, 'trees',
# *and* migrates the |Taxon| objects of the tree over to the |TaxonNamespace|
# of 'trees'. This will break things if the tree is contained in another
# |TreeList| with a different |TaxonNamespace|!
trees.append(tree)

# In contrast to before, the |TaxonNamespace| of 'tree' is not the same
# as the |TaxonNamespace| of 'trees. The |Taxon| objects have been imported
# and/or remapped based on their label.
assert tree.taxon_namespace is trees.taxon_namespace
print(id(original_tree_taxon_namespace))
