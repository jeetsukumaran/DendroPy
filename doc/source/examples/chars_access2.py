import dendropy

dna = dendropy.DnaCharacterMatrix.get(
        path="primates.chars.nexus",
        schema="nexus")

# iterate over taxa
for taxon in dna:
    print("{}: {}".format(taxon.label, dna[taxon]))

# iterate over the sequences
for seq in dna.values():
    print(seq)

# iterate over taxon/sequence pairs
for taxon, seq in dna.items():
    print("{}: {}".format(taxon.label, seq))

