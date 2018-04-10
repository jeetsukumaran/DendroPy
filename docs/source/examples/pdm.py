import dendropy

tree = dendropy.Tree.get(
    path="pythonidae.mle.nex",
    schema="nexus")
pdc = tree.phylogenetic_distance_matrix()
for i, t1 in enumerate(tree.taxon_namespace[:-1]):
    for t2 in tree.taxon_namespace[i+1:]:
        print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))
