#!/usr/bin/env python

if __name__ == '__main__':
    import sys
    from dendropy import dataio
    from dendropy.treestats import pybus_harvey_gamma

    from optparse import OptionParser
    usage = """Calculates the Pybus-Harvey Gamma statistic for a tree specified as
    newick string.
"""
    (options, args) = parser.parse_args()
    if len(args) > 1:
        sys.exit("At most one argument (a newick tree string with branch lengths) can be specified. Got: %s" % ' '.join(args))
    if len(args) == 1:
        newick = args[0]
    else:
        newick = sys.stdin.readline()

    tree = dataio.trees_from_string(string=newick, schema="NEWICK")[0]
    try:
        sys.stdout.write("%f\n" % pybus_harvey_gamma(tree))
    except ValueError as x:
        sys.exit(str(x))


