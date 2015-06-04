import dendropy
from dendropy.calculate import treecompare

tree_str1 = "((A,B),C);"

tree_list1 = dendropy.TreeList()
tree_list1.read(data=tree_str1, schema="newick")
tree_list2 = dendropy.TreeList(taxon_namespace=tree_list1.taxon_namespace)
tree_list2.read(data=tree_str1, schema="newick")

# Results in: 0
print(treecompare.symmetric_difference(tree_list1[0], tree_list2[0]))
