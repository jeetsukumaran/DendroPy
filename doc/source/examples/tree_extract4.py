#! /usr/bin/env python

import dendropy

class PhyloNode(dendropy.Node):
    pass

class PhyloTree(dendropy.Tree):
    @classmethod
    def node_factory(cls, *args, **kwargs):
        return PhyloNode(*args, **kwargs)

tree0 = dendropy.Tree.get(
        data="(a,(b,(c,d)));",
        schema="newick",
        )

# Default: use dendropy.Tree for tree
# and dendropy.Node for node
tree2 = tree0.extract_tree(tree_factory=PhyloTree)
print(type(tree2)) # PhyloTree
for nd in tree2:
    print(type(nd)) # PhyloNode

# Node factory defaults to ``node_factory`` method
# of instance returned by ``tree_factory`` if
# ``tree_factory is specified.
tree2 = tree0.extract_tree(tree_factory=PhyloTree)
print(type(tree2)) # PhyloTree
for nd in tree2:
    print(type(nd)) # PhyloNode

# equivalent to above
tree3 = tree0.extract_tree(tree_factory=PhyloTree,
        node_factory=PhyloNode)
print(type(tree2))
for nd in tree2:
    print(type(nd))

# Use dendropy.Tree for tree but
# PhyloNode for node
tree4 = tree0.extract_tree(node_factory=PhyloNode)
print(type(tree2))
for nd in tree2:
    print(type(nd))
