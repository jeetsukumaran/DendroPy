import dendropy
dataset = dendropy.DataSet.get(
    path="primates.cc.combined.nex",
    schema="nexus")
tree = dataset.tree_lists[0][0]
chars = dataset.char_matrices[0]

