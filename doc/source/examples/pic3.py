import dendropy
taxa = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(
    path="primates.cc.tre",
    schema="newick",
    taxon_namespace=taxa)
chars = dendropy.ContinuousCharacterMatrix.get(
    path="primates.cc.nex",
    schema="nexus",
    taxon_namespace=taxa)
