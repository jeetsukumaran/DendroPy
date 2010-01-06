#! /usr/bin/env python

import dendropy

# tree assumed to be unrooted unless '[&R]' is specified
tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus")

# forces tree to be rooted
tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus", as_rooted=True)

# forces tree to be unrooted
tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus", as_rooted=False)

# forces tree to be unrooted
tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus", as_unrooted=True)

# forces tree to be rooted
tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus", as_unrooted=False)

# also applies to constructors ...
tree = dendropy.Tree(
        stream=open("pythonidae.mle.nex", "rU"),
        schema="nexus",
        as_rooted=True)

# and 'read_from_*' methods
tree = dendropy.Tree()
tree.read_from_path("pythonidae.mle.nex", "nexus", as_rooted=True)

# and TreeList constructor, 'get_from_*', and 'read_from_*' methods
tree_list = dendropy.TreeList(
        stream=open("pythonidae.mcmc.nex", "rU"),
        schema="nexus",
        as_rooted=True)

tree_list = dendropy.TreeList.get_from_path(
        "pythonidae.mcmc.nex",
        "nexus", as_rooted=True)

tree_list = dendropy.TreeList()
tree_list.read_from_path("pythonidae.mcmc.nex", "nexus", as_rooted=True)
