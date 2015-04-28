import dendropy
from dendropy.calculate import treecompare

tree_str1 = "((A,B),C);"
tree_str2 = "((A,B),C);"

ds = dendropy.DataSet()
ds.read(data=tree_str1, schema="newick")
ds.read(data=tree_str1, schema="newick")

print(len(ds.taxon_namespaces))
print(ds.tree_lists[0].taxon_namespace is ds.tree_lists[1].taxon_namespace)
print(ds.tree_lists[0].taxon_namespace[0] is ds.tree_lists[1].taxon_namespace[0])

# Results in:
# 2
# False
# False
