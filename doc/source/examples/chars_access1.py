import dendropy

dna = dendropy.DnaCharacterMatrix.get(
        path="primates.chars.nexus",
        schema="nexus")

# access by dereferencing taxon label
s1 = dna["Macaca sylvanus"]

# access by taxon index
s2 = dna[0]
s3 = dna[4]
s4 = dna[-2]

# access by taxon instance
t = dna.taxon_namespace.get_taxon(label="Macaca sylvanus")
s5 = dna[t]

