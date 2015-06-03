import dendropy
taxa = dendropy.TaxonNamespace()
paths = [
        "primates.chars.subsets-1stpos.nexus",
        "primates.chars.subsets-2ndpos.nexus",
        "primates.chars.subsets-3rdpos.nexus",
        ]
d_all = dendropy.DnaCharacterMatrix.concatenate_from_paths(
        paths=paths,
        schema="nexus")
print("d_all: {} sequences, {} characters".format(len(d_all), d_all.max_sequence_size))
print("Subsets: {}".format(d_all.character_subsets))
