import dendropy
from dendropy.calculate import treecompare

tree_str1 = "((A,B),C);"
tree_str2 = "((A,B),C);"

tree1 = dendropy.Tree.get(data=tree_str1, schema="newick")
tree2 = dendropy.Tree.get(data=tree_str1, schema="newick")
for taxon in tree1.taxon_namespace:
    if taxon in tree2.taxon_namespace:
        # this branch is never visited
        print("Taxon '{}': found in both trees".format(taxon.label))

## Following will result in:
## dendropy.utility.error.TaxonNamespaceIdentityError: Non-identical taxon namespace references: ...
# print(treecompare.symmetric_difference(tree1, tree2))
