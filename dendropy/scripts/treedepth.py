#!/usr/bin/env python
#   Reports the root to tip distance for an ultrametric tree

import sys
from dendropy import dataio
import dendropy


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
    dendropy.trees.add_depth_to_nodes(tree, options.prec)

    sys.stdout.write("%f\n" % tree.seed_node.depth)


