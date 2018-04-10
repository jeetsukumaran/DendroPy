import dendropy
from dendropy.calculate import treecompare

tree_str1 = "((A,B),C);"
tree_str2 = "((A,B),C);"

ds2 = dendropy.DataSet()
ds2.read(data=tree_str1, schema="newick")
ds2.read(
        data=tree_str1,
        schema="newick",
        taxon_namespace=ds2.tree_lists[0].taxon_namespace)

print(len(ds2.taxon_namespaces))
print(ds2.tree_lists[0].taxon_namespace is ds2.tree_lists[1].taxon_namespace)
print(ds2.tree_lists[0].taxon_namespace[0] is ds2.tree_lists[1].taxon_namespace[0])

# Results in:
# 1
# True
# True

ds2 = dendropy.DataSet()
tns = ds2.new_taxon_namespace()
ds2.read(
        data=tree_str1,
        schema="newick",
        taxon_namespace=tns)
ds2.read(
        data=tree_str1,
        schema="newick",
        taxon_namespace=tns)

print(len(ds2.taxon_namespaces))
print(ds2.tree_lists[0].taxon_namespace is ds2.tree_lists[1].taxon_namespace)
print(ds2.tree_lists[0].taxon_namespace[0] is ds2.tree_lists[1].taxon_namespace[0])

# Results in:
# 1
# True
# True

