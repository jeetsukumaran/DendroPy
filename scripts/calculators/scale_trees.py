#!/usr/bin/env python
#   Reports the root to tip distance for an ultrametric tree
import sys
import StringIO
from dendropy import DataSet
from dendropy.treemanip import scale_edges

if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-m", "--multiplier", dest="multip", default=1.0,
        type="float",
        help="The multiplier used for every branch length.")
    parser.add_option("-n", "--nexus", dest="schema", action="store_const", const="NEXUS", default="NEWICK",
        help="Tree is in NEXUS schema.")
    (options, args) = parser.parse_args()


    if len(args) > 1:
        sys.exit("At most one argument (a newick tree string with branch lengths) can be specified")
    if len(args) == 1:
        s = open(args[0], 'rU')
    else:
        newick = sys.stdin.read()
        s = StringIO.StringIO(newick)

    multip = options.multip
    d = DataSet()
    
    d.read(s, schema=options.schema, rooted=True)
    if len(d.tree_lists) == 0:
        sys.exit("No trees found in file.")
    for tb in d.tree_lists:
        for tree in tb:
            scale_edges(tree, multip)
    d.write_to_stream(sys.stdout, schema=options.schema)
