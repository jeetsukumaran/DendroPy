import dendropy
from dendropy.calculate import treecompare

trees = dendropy.TreeList.get(
        path="pythonidae.random.bd0301.tre",
        schema="nexus")

for tree in trees:
    print(tree.as_string("newick"))

print(len(trees))

print(trees[4].as_string("nexus"))
print(treecompare.robinson_foulds_distance(trees[0], trees[1]))
print(treecompare.weighted_robinson_foulds_distance(trees[0], trees[1]))

first_10_trees = trees[:10]
last_10_trees = trees[-10:]

# Note that the TaxonNamespace is propogated to slices
assert first_10_trees.taxon_namespace is trees.taxon_namespace
assert first_10_trees.taxon_namespace is trees.taxon_namespace


print(id(trees[4]))
print(id(trees[5]))
trees[4] = trees[5]
print(id(trees[4]))
print(id(trees[5]))
print(trees[4] in trees)

trees.remove(trees[-1])
tx = trees.pop()
print(trees.index(trees[0]))

trees.sort(key=lambda t:t.label)
trees.reverse()
trees.clear()
