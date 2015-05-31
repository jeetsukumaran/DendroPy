import dendropy
from dendropy.calculate import treecompare

s1 = "(a,(b,(c,d)));"
s2 = "((a,b),(c,d));"

tns = dendropy.TaxonNamespace()

unrooted_tree1 = dendropy.Tree.get(
        data=s1,
        schema='newick',
        taxon_namespace=tns)
unrooted_tree2 = dendropy.Tree.get(
        data=s2,
        schema='newick',
        taxon_namespace=tns)

rooted_tree1 = dendropy.Tree.get(
        data=s1,
        schema='newick',
        rooting="force-rooted",
        taxon_namespace=tns)
rooted_tree2 = dendropy.Tree.get(
        data=s2,
        schema='newick',
        rooting="force-rooted",
        taxon_namespace=tns)

## Unweighted Robinson-Foulds distance (rooted) = 2
print(treecompare.symmetric_difference(rooted_tree1, rooted_tree2))
## Unweighted Robinson-Foulds distance (unrooted) = 0
print(treecompare.symmetric_difference(unrooted_tree1, unrooted_tree2))
## Unweighted Robinson-Foulds distance (rooted1 to unrooted2) = 3
print(treecompare.symmetric_difference(rooted_tree1, unrooted_tree2))
## Unweighted Robinson-Foulds distance (unrooted1 to rooted2) = 5
print(treecompare.symmetric_difference(unrooted_tree1, rooted_tree2))

