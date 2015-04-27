#! /usr/bin/env python

import dendropy
from dendropy.calculate import popgenstat

seqs = dendropy.DnaCharacterMatrix.get(
        path="orti1994.nex",
        schema="nexus")
p1 = []
p2 = []
for idx, t in enumerate(seqs.taxon_namespace):
    if t.label.startswith('EPAC'):
        p1.append(seqs[t])
    else:
        p2.append(seqs[t])
pp = popgenstat.PopulationPairSummaryStatistics(p1, p2)

print('Average number of pairwise differences (total): %s' \
    % pp.average_number_of_pairwise_differences)
print('Average number of pairwise differences (between populations): %s' \
    % pp.average_number_of_pairwise_differences_between)
print('Average number of pairwise differences (within populations): %s' \
    % pp.average_number_of_pairwise_differences_within)
print('Average number of pairwise differences (net): %s' \
    % pp.average_number_of_pairwise_differences_net)
print('Number of segregating sites: %s' \
    % pp.num_segregating_sites)
print("Watterson's theta: %s" \
    % pp.wattersons_theta)
print("Wakeley's Psi: %s" \
    % pp.wakeleys_psi)
print("Tajima's D: %s" \
    % pp.tajimas_d)
