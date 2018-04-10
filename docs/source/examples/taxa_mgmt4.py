import dendropy
from dendropy.calculate import treecompare

tree_str1 = "((A,B),C);"

tree_list1 = dendropy.TreeList()
tree_list1.read(data=tree_str1, schema="newick")
tree_list2 = dendropy.TreeList()
tree_list2.read(data=tree_str1, schema="newick")

for taxon in tree_list1.taxon_namespace:
    if taxon in tree_list2.taxon_namespace:
        # this branch is never visited
        print("Taxon '{}': found in both trees".format(taxon.label))

## Following will result in:
## dendropy.utility.error.TaxonNamespaceIdentityError: Non-identical taxon namespace references: ...
# print(treecompare.symmetric_difference(tree_list1[0], tree_list2[0]))
