#!/usr/bin/env python
import sys
import copy
import logging
import itertools
from dendropy.utility.messaging import get_logger
from dendropy.treesplit import encode_splits
from dendropy.treemanip import collapse_edge
from dendropy.treecalc import symmetric_difference
from dendropy import format_split, Edge, TaxonSet, DataSet

_LOG = get_logger('scripts.long_branch_symmdiff')
verbose = False

def long_branch_symmdiff(trees_to_compare, edge_len_threshold, copy_trees=False, rooted=False):
    """Returns matrix of the symmetric_differences between trees after all
    internal edges with lengths < `edge_len_threshold` have been collapsed.

    If `copy_trees` is True then the trees will be copied first (if False, then
        the trees may will have their short edges collapsed on exit).
    """
    if copy_trees:
        tree_list = [copy.copy(i) for i in trees_to_compare]
    else:
        tree_list = list(trees_to_compare)

    n_trees = len(tree_list)
    _LOG.debug('%d Trees to compare:\n%s\n' % (n_trees, '\n'.join([str(i) for i in tree_list])))
    if n_trees < 2:
        return [0 for t in tree_list]

    f_r = []
    for tree in tree_list:
        to_collapse = []
        encode_splits(tree)
        for edge in tree.preorder_edge_iter(filter_fn=Edge.is_internal):
            elen = edge.length
            if elen is not None and elen < edge_len_threshold:
                to_collapse.append(edge)
        for edge in to_collapse:
            collapse_edge(edge)
        f_r.append(tree.is_rooted)
        tree.is_rooted = bool(rooted)
        encode_splits(tree)

    sd_row = [0]*n_trees
    sd_mat = [list(sd_row) for i in xrange(n_trees)]
    for i, tree_one in enumerate(tree_list[:-1]):
        for col_count, tree_two in enumerate(tree_list[1+i:]):
            j = i + 1 + col_count
            sd = symmetric_difference(tree_one, tree_two)
            sd_mat[i][j] = sd
            sd_mat[j][i] = sd

    if not copy_trees:
        for r, tree in itertools.izip(f_r, tree_list):
            tree.is_rooted = r
    return sd_mat


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-s', '--schema', dest='schema',
                        type='str', default="newick", help='The format/schema of the input data (default is newick)')
    parser.add_option('-c', '--cutoff', dest='cutoff',
                        type='str', default=0.0, help='The minimum edge length (any branches shorter than this will be collapsed).')
    parser.add_option('-p', '--paup-style', dest='paup',
                        action="store_true", default=False, help="Produce an output in the same format as PAUP's TreeDist command.")
    (options, args) = parser.parse_args()
    if len(args) == 0:
        sys.exit("Expecting a filename as an argument")
    schema = options.schema.upper()
    try:
        cutoff = int(options.cutoff)
    except ValueError:
        try:
            cutoff = float(options.cutoff)
        except ValueError:
            sys.exit('Expecting the cutoff to be a number found "%s"' % options.cutoff)

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

    sd_mat = long_branch_symmdiff(trees, cutoff)
    o = sys.stdout
    if options.paup:
        o.write("%s\n" % "\t".join(["tree"] + [str(1+i) for i in xrange(len(sd_mat))]))
        for n, row in enumerate(sd_mat):
            o.write("%d\t%s\n" % ((n + 1), "\t".join([str(i) for i in row[:1 + n]])))
    else:
        for row in sd_mat:
            o.write("%s\n" % "\t".join([str(i) for i in row]))


