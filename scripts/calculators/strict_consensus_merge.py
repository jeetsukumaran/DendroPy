#!/usr/bin/env python
import sys
import copy
import logging

from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)
from dendropy import DataSet, TaxonSet
verbose = False
from dendropy.scm import inplace_strict_consensus_merge

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-s', '--schema', dest='schema',
                        type='str', default="newick", help='The format/schema of the input data')
    parser.add_option('-g', '--gordon', dest='gordons',
                        action="store_true", default=False, help="Specify to use the Gordon's strict consensus")
    (options, args) = parser.parse_args()
    if len(args) == 0:
        sys.exit("Expecting a filename as an argument")
    schema = options.schema.upper()

    trees = []
    taxon_set = TaxonSet()
    dataset = DataSet(taxon_set=taxon_set)
    if schema == "PHYLIP":
        schema = "NEWICK"
    for f in args:
        fo = open(f, "rU")
        dataset.read(stream=fo, schema=schema)
    for tl in dataset.tree_lists:
        trees.extend(tl)

    o = inplace_strict_consensus_merge(trees, gordons_supertree=options.gordons)
    sys.stdout.write("%s;\n" % str(o))


