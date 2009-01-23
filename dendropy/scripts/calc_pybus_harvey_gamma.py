#!/usr/bin/env python

if __name__ == '__main__':
    import sys
    from dendropy import dataio
    from dendropy.treecalc import pybus_harvey_gamma

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

    tree = dataio.trees_from_string(string=newick, format="NEWICK")[0]
    try:
        sys.stdout.write("%f\n" % pybus_harvey_gamma(tree))
    except ValueError as x:
        sys.exit(str(x))


