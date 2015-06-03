#! /usr/bin/env python

import dendropy

# tree assumed to be unrooted unless '[&R]' is specified
tree = dendropy.Tree.get(path="pythonidae.mle.nex", schema="nexus")

# same as above, explicit
tree = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus",
        rooting="default-unrooted",
        )

# forces tree to be rooted
tree = dendropy.Tree.get(path="pythonidae.mle.nex",
        schema="nexus",
        rooting="force-rooted")

# forces tree to be unrooted
tree = dendropy.Tree.get(path="pythonidae.mle.nex",
        schema="nexus",
        rooting="force-unrooted")

# forces trees to be rooted
tree_list = dendropy.TreeList.get(
        path="pythonidae.mcmc.nex",
        schema="nexus",
        rooting="force-rooted")

# forces trees to default to rooted, unless '[&U]' is specified
tree_list = dendropy.TreeList()
tree_list.read(
    path="pythonidae.mcmc.nex",
    schema="nexus",
    rooting="default-rooted")
