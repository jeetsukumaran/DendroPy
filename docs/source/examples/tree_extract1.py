import dendropy
from dendropy.calculate import treecompare

tree0 = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus")
for idx, nd in enumerate(tree0):
    nd.label = "hello, world{}".format(idx)
    nd.edge.label = "world, hello{}".format(idx)
    nd.annotations["color"] = "blue"
    nd.edge.annotations["taste"] = "sweet"
tree1 = tree0.extract_tree()

assert tree0.taxon_namespace is tree1.taxon_namespace
assert treecompare.weighted_robinson_foulds_distance(
        tree0, tree1) == 0.0

for nd in tree1:
    original_node = nd.extraction_source
    print("{} on extracted tree corresponds to {} on original tree".format(
        nd, original_node))
    ## basic attributes copied
    assert nd.label == original_node.label
    assert nd.edge.label == original_node.edge.label
    assert nd.edge.length == original_node.edge.length
    ## but not annotations
    assert len(nd.annotations) == 0 and len(original_node.annotations) > 0
    assert len(nd.edge.annotations) == 0 and len(original_node.edge.annotations) > 0


