#! /usr/bin/env python

import dendropy

# tree assumed to be unrooted unless '[&R]' is specified
tree = dendropy.Tree.get(path="pythonidae.mle.nex", schema="nexus")

# forces tree to be rooted
tree = dendropy.Tree.get(path="pythonidae.mle.nex",
        schema="nexus",
        rooting="force-rooted")

# forces tree to be unrooted
tree = dendropy.Tree.get(path="pythonidae.mle.nex",
        schema="nexus",
        rooting="force-unrooted")

tree = dendropy.Tree()
tree.read(
        path="pythonidae.mle.nex",
        schema="nexus",
        rooting="force-rooted")


tree_list = dendropy.TreeList.get(
        path="pythonidae.mcmc.nex",
        schema="nexus",
        rooting="force-rooted")

tree_list = dendropy.TreeList()
tree_list.read(
    path="pythonidae.mcmc.nex",
    schema="nexus",
    rooting="force-rooted")
