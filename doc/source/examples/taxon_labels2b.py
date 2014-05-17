#! /usr/bin/env python

import dendropy

nexus1 = """
#NEXUS

begin taxa;
    dimensions ntax=2;
    taxlabels Python_regius Python_sebae;
end;

begin characters;
    dimensions nchar=5;
    format datatype=dna gap=- missing=? matchchar=.;
    matrix
        Python_regius ACGTA
        Python_sebae   ACGTA
    ;
end;
"""

fasta1 = """
>Python regius
AAAA
>Python sebae
ACGT
"""

d = dendropy.DataSet()
d.attach_taxon_set()
d.read_from_string(nexus1, "nexus")
d.read_from_string(fasta1, "dnafasta")
print(d.taxon_sets[0].description(2))
