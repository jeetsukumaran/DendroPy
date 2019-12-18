#! /usr/bin/env python
# -*- coding: utf-8 -*-

import dendropy

seqstr = """\
#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=13;
    TAXLABELS a1 a2 a3 b1 b2 b3 c1 c2 c3 c4 c5 d1 d2;
END;
BEGIN CHARACTERS;
    DIMENSIONS NCHAR=7;
    FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=.;
    MATRIX
        a1 ACCTTTG
        a2 ACCTTTG
        a3 ACCTTTG
        b1 ATCTTTG
        b2 ATCTTTG
        b3 ACCTTTG
        c1 ACCCTTG
        c2 ACCCTTG
        c3 ACCCTTG
        c4 ACCCTTG
        c5 ACCCTTG
        d1 ACAAAAG
        d2 ACCAAAG
    ;
END
"""
seqs = dendropy.DnaCharacterMatrix.get_from_string(seqstr, 'nexus')
taxon_namespace = seqs.taxon_namespace

tax_parts = taxon_namespace.partition(membership_func=lambda x: x.label[0])

for s in tax_parts.subsets():
    print(s.description())
