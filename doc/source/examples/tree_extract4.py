#! /usr/bin/env python

import dendropy


# Our custom node class. We do not actually need to derive from dendropy.Node.
# Doing so, though, ensures that all the required methods exist instead of
# having to write them ourselves.
class PhyloNode(dendropy.Node):
    pass

# Our custom tree class. Note that the definition of node factory
# (class-)method here ensures that nodes created are of the type we we want.
# Many DendroPy methods, include 'Tree.extract_tree()' check and preferentially
# use the return value of tree's class' ``node_factory()`` method when building
# trees. Note also that we do not actually need to derive from dendropy.Tree.
# Doing so, though, ensures that all the required methods exist instead of
# having to write them ourselves.
class PhyloTree(dendropy.Tree):
    @classmethod
    def node_factory(cls, *args, **kwargs):
        return PhyloNode(*args, **kwargs)

# The original tree using dendropy.Tree for the tree and dendropy.Node for the
# nodes.
tree0 = dendropy.Tree.get(
        data="(a,(b,(c,d)));",
        schema="newick",
        )
print(type(tree0)) # dendropy.Tree
for nd in tree0:
    print(type(nd)) # dendropy.Node

# Default extraction: use dendropy.Tree for tree
# and dendropy.Node for node
tree1 = tree0.extract_tree()
print(type(tree1)) # dendropy.Tree
for nd in tree1:
    print(type(nd)) # dendropy.Node

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
print(type(tree3)) # PhyloTree
for nd in tree3:
    print(type(nd)) # PhyloNode

# Use dendropy.Tree for tree but
# PhyloNode for node
tree4 = tree0.extract_tree(node_factory=PhyloNode)
print(type(tree4)) # dendropy.Tree
for nd in tree4:
    print(type(nd)) # PhyloNode
