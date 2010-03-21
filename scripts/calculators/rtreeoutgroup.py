#!/usr/bin/env python
def rtree_outgroup_labels(tree):
    """Takes a tree (which will be treated as rooted), and returns a list of labels
    of the nodes that would serve as the "outgroup" if you were to define the
    largest clade in the tree to be the "ingroup".

    Adds "n_leaves_under" and "in_biggest" attributes to the nodes in the tree."""
    node = None
    # add an n_leaves_under attribute
    for node in tree.postorder_node_iter():
        e = node.edge
        p = getattr(e, "tail_node", None)
        if p:
            p.n_leaves_under = getattr(p, "n_leaves_under", 0) +  getattr(node, "n_leaves_under", 1)

    # find the child of the root with the largest number of descendants
    seed_node = tree.seed_node
    ch = seed_node.child_nodes()
    f = ch[0]
    f.in_biggest = False
    biggest_clade, bc_size = f, getattr(f, "n_leaves_under", 1)
    for nd in ch[1:]:
        nk = getattr(nd, "n_leaves_under", 1)
        if nd > bc_size:
            biggest_clade, bc_size = nd, nk
        nd.in_biggest = False
    # Mark the biggest clade, and accumulate out all unmarked leaf names
    biggest_clade.in_biggest = True
    outgroup_labels = []
    for node in tree.preorder_node_iter():
        par = node.parent_node
        if node == seed_node or par == seed_node:
            continue
        node.in_biggest = par.in_biggest
        if (not node.in_biggest) and (not node.child_nodes()):
            outgroup_labels.append(node.label)
    return outgroup_labels

if __name__ == '__main__':
    import sys
    from dendropy import dataio

    from optparse import OptionParser
    usage = """Takes a tree (which will be treated as rooted), and reports the list of taxa
    as the "outgroup" such that the largest clade in the tree is the "ingroup".
This script can be useful when you are simulating rooted trees and need to do
   a rooted inference procedure using software (such as PAUP*) that requires
   you to specify an outgroup when you want to score under a clock model.
"""
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--sep", dest="sep", default="\n",
        type="str",
        help="The string that separates taxa listed in the outgroup")
    (options, args) = parser.parse_args()
    if len(args) > 1:
        sys.exit("At most one argument (a newick tree string with branch lengths) can be specified. Got: %s" % ' '.join(args))
    if len(args) == 1:
        newick = args[0]
    else:
        newick = sys.stdin.readline()

    tree = dataio.trees_from_string(string=newick, schema="NEWICK")[0]
    outgroup_labels = rtree_outgroup_labels(tree)
    sep = options.sep
    sys.stdout.write("%s\n" % sep.join(outgroup_labels))


