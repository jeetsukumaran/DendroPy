#! /usr/bin/env python

import dendropy

nexus1 = """
#NEXUS

begin taxa;
    dimensions ntax=2;
    taxlabels 'Python regius' 'Python sebae';
end;

begin characters;
    dimensions nchar=5;
    format datatype=dna gap=- missing=? matchchar=.;
    matrix
        'Python regius' ACGTA
        'Python sebae'   ACGTA
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
tns = d.new_taxon_namespace()
d.attach_taxon_namespace(tns)
d.read(data=nexus1, schema="nexus")
d.read(data=fasta1, schema="fasta", data_type="dna")
print(d.taxon_namespaces[0].description(2))
