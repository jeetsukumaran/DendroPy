import dendropy
from dendropy.calculate import treemeasure

tree = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus")
pdc = tree.phylogenetic_distance_matrix()
for taxon1 in tree.taxon_namespace:
    for taxon2 in tree.taxon_namespace:
        mrca = pdc.mrca(taxon1, taxon2)
        print("MRCA of '{}' and '{}' is: {}".format(taxon1.label, taxon2.label, mrca.description()))
