import dendropy
taxa = dendropy.TaxonNamespace()
d1 = dendropy.DnaCharacterMatrix.get(
        path="primates.chars.subsets-1stpos.nexus",
        schema="nexus",
        taxon_namespace=taxa)
print("d1: {} sequences, {} characters".format(len(d1), d1.max_sequence_size))
d2 = dendropy.DnaCharacterMatrix.get(
        path="primates.chars.subsets-2ndpos.nexus",
        schema="nexus",
        taxon_namespace=taxa)
print("d2: {} sequences, {} characters".format(len(d2), d2.max_sequence_size))
d3 = dendropy.DnaCharacterMatrix.get(
        path="primates.chars.subsets-3rdpos.nexus",
        schema="nexus",
        taxon_namespace=taxa)
print("d3: {} sequences, {} characters".format(len(d3), d3.max_sequence_size))
d_all = dendropy.DnaCharacterMatrix.concatenate([d1,d2,d3])
print("d_all: {} sequences, {} characters".format(len(d_all), d_all.max_sequence_size))
print("Subsets: {}".format(d_all.character_subsets))
