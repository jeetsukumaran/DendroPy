import dendropy

tree = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus")
pdm = tree.phylogenetic_distance_matrix()
for idx1, taxon1 in enumerate(tree.taxon_namespace):
    for taxon2 in tree.taxon_namespace:
        mrca = pdm.mrca(taxon1, taxon2)
        weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)
        unweighted_patristic_distance = pdm.path_edge_count(taxon1, taxon2)
        print("'{}' vs '{}': {} (distance (weighted-edges, unweighted-edges) = {}, {})".format(
            taxon1.label,
            taxon2.label,
            mrca.bipartition.split_as_bitstring(),
            weighted_patristic_distance,
            unweighted_patristic_distance))
