# /usr/bin/env python

import dendropy
from dendropy.utility import container
from dendropy.utility.textprocessing import StringIO

phylogeny_str = """\
[&R]((((spA:0.10,spB:0.10):0.67,(((spC:0.08,spD:0.08):0.24,(spE:0.13,spF:0.13):0.19):0.40,(spG:0.12,((spH:0.04,spI:0.04):0.07,spJ:0.12):0.00):0.60):0.05):0.78,((spK:0.05,(spL:0.04,spM:0.04):0.01):0.31,spN:0.37):1.19):1.22,spO:2.79);
"""

assemblage_data_table_str = """\
.,spA,spB,spC,spD,spE,spF,spG,spH,spI,spJ,spK,spL,spM,spN,spO
C1,15,0,10,0,12,0,7,0,0,0,0,0,2,2,1
C2,0,25,0,6,4,0,2,0,0,3,23,0,0,0,0
C3,30,0,10,0,0,10,9,4,0,0,0,0,0,10,0
C4,0,0,0,0,0,0,0,10,20,1,2,25,4,0,0
C5,0,0,0,0,0,0,0,35,14,10,0,0,0,0,0
"""

# read the tree
tree = dendropy.Tree.get(
        data=phylogeny_str,
        schema="newick",
        )

# obtain the PhylogeneticDistanceMatrix corresponding to the taxon-to-taxon
# distances of the above tree
pdm = tree.phylogenetic_distance_matrix()


# read the assemblage data into a table,
# being sure to specify an appropriate data type!
assemblage_data = container.DataTable.from_csv(
        src=StringIO(assemblage_data_table_str),
        default_data_type=int)

# generate the assemblage definitions
assemblage_names = []
assemblage_memberships = []
for row_name in assemblage_data.row_name_iter():
    assemblage_names.append(row_name)
    member_labels = set([col_name for col_name in assemblage_data.column_name_iter() if assemblage_data[row_name, col_name] > 0])
    member_taxa = set([t for t in pdm.taxon_namespace if t.label in member_labels])
    assemblage_memberships.append(member_taxa)

# calculate the SES statistics for each assemblage
results_mpd = pdm.standardized_effect_size_mean_pairwise_distance(
        assemblage_memberships=assemblage_memberships)

# inspect the results
print("Phylogenetic Community Standardized Effect Size Statistics:")
assert len(results_mpd) == len(assemblage_memberships)
assert len(results_mpd) == len(assemblage_names)
for assemblage_name, assemblage_membership, result in zip(assemblage_names, assemblage_memberships, results_mpd, ):
    print("# Assemblage '{}' ({})".format(
        assemblage_name,
        sorted([t.label for t in assemblage_membership])))
    print("   -     MPD: {}".format(result.obs))
    print("   - SES MPD: {}".format(result.z))
    print("   - p-value: {}".format(result.p))
