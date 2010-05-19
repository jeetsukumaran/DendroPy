#!/usr/bin/env python
#   Reports the root to tip distance for an ultrametric tree

import sys
from dendropy import DataSet
import StringIO

def simple_test_tree():
    from dendropy import treegen
    t = treegen.pure_kingman(treegen.random_taxa_block(10))
    sys.stdout.write(t.compose_newick() + "\n")

if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-p", "--prec", dest="prec", default=0.00001,
        type="float",
        help="The precision of the comparison that the node depth is equal regardless of tip to node path.")
    parser.add_option("--nexus", dest="schema", action="store_const", const="NEXUS", default="NEWICK",
        help="Tree is in NEXUS schema.")
    parser.add_option("--generate-test-tree", action="store_true", dest="gen_test_tree",
        help="Generate ultrametric tree to test function.")
    (options, args) = parser.parse_args()

    if options.gen_test_tree:
        simple_test_tree()
        exit(0)

    if len(args) > 1:
        sys.exit("At most one argument (a newick tree string with branch lengths) can be specified")
    if len(args) == 1:
        s = args[0]
    else:
        newick = sys.stdin.read()
        s = StringIO.StringIO(newick)

    prec = options.prec
    d = DataSet()
    
    d.read(s, schema=options.schema, rooted=True)
    if len(d.trees_blocks) == 0:
        sys.exit("No trees found in file.")
    tree = d.trees_blocks[0][0]
    tree.calc_node_ages(attr_name='age', check_prec=options.prec)

    sys.stdout.write("%f\n" % tree.seed_node.age)


