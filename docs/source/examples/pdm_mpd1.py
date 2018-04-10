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

# print the communities
print("Assemblage Memberships:")
for row_name in assemblage_data.row_name_iter():
    members = [col_name for col_name in assemblage_data.column_name_iter() if assemblage_data[row_name, col_name] > 0]
    print("{}: {}".format(row_name, members))

# calculate the statistics for each community:
print("Phylogenetic Community Statistics:")
for row_name in assemblage_data.row_name_iter():
    members = [col_name for col_name in assemblage_data.column_name_iter() if assemblage_data[row_name, col_name] > 0]
    filter_fn = lambda taxon: taxon.label in set(members)
    mpd = pdm.mean_pairwise_distance(filter_fn=filter_fn)
    mntd = pdm.mean_nearest_taxon_distance(filter_fn=filter_fn)
    print("{}: MPD={}, MNTD={}".format(row_name, mpd, mntd))
