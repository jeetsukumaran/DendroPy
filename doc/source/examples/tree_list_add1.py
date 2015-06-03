import dendropy

trees = dendropy.TreeList()
trees.read(path="sometrees.nex", schema="nexus", tree_offset=10)
trees.read(data="(A,(B,C));((A,B),C);", schema="newick")

