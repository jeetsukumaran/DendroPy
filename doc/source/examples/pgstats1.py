#! /usr/bin/env python

import dendropy

seqs = dendropy.DnaCharacterMatrix.get_from_path("pacfish.nex", schema="nexus")
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_set):
    if t.label.startswith('EPAC'):
        p1.append(seqs[t])
    else:
        p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)

