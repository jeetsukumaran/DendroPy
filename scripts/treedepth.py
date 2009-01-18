#!/usr/bin/env python
#   Reports the root to tip distance for an ultrametric tree

import sys
from dendropy import dataio

if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-p", "--prec", dest="prec", default=0.00001, 
        type="float",
        help="The precision of the comparison that the node depth is equal regardless of tip to node path.")
    (options, args) = parser.parse_args()
    if len(args) > 1:
        sys.exit("At most one argument (a newick tree string with branch lengths) can be specified")
    if len(args) == 1:
        newick = args[0]
    else:
        newick = sys.stdin.readline()

    prec = options.prec
    tree = dataio.trees_from_string(string=newick, format="NEWICK")[0]
    node = None
    for node in tree.postorder_node_iter():
        ch = node.child_nodes()
        if len(ch) == 0:
            node.depth = 0.0
        else:
            first_child = ch[0]
            node.depth = first_child.depth + first_child.edge.length
            last_child = ch[-1]
            for nnd in ch[1:]:
                ocnd = nnd.depth + nnd.edge.length
                if abs(node.depth - ocnd) > prec:
                    sys.exit("Tree is not ultrametric (%f != %f)" % (node.depth, ocnd))
    if node is None:
        sys.exit("Empty tree encountered")
    sys.stdout.write("%f\n" % node.depth)


