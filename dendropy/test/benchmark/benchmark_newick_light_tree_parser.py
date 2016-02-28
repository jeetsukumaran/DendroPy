#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Benchmarking NEWICK tree parsing using a light tree data model.
"""

import sys
import os
import timeit
import argparse
from dendropy.utility.textprocessing import StringIO
from dendropy.utility import messaging
from dendropy.test.support import pathmap

from dendropy.dataio import nexusprocessing
from dendropy.dataio import newickreader

TREE_FILENAMES = [
    "APG_Angiosperms.newick",
    "GEBA.tree.newick",
    "feb032009.trees.newick",
    "Bininda-emonds_2007_mammals.newick",
    "Jetz_et_al_2012_Aves.sample.tree.newick",
    "Smith_2001_angiosperms.newick",
        ]

class Tree(object):

    def __init__(self):
        self.seed_node = Node()

    def node_factory(self):
        return Node()

class Node(object):

    def __init__(self,
            taxon=None,
            label=None,
            edge=None,
            edge_length=None):
        self.age = None
        self._edge = None
        self._child_nodes = []
        self._parent_node = None
        if edge is not None:
            self.edge = edge
        else:
            self.edge = Edge(head_node=self,
                    length=edge_length)
        self.label = label

    def add_child(self, node):
        self._child_nodes.append(node)

    ###########################################################################
    ## Hacked-in NEWICK representation.

    def as_newick_string(self, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.
        For production purposes, use the the full-fledged 'as_string()'
        method of the object.
        """
        out = io.StringIO()
        self.write_newick(out, **kwargs)
        return out.getvalue()

    def write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.  For
        production purposes, use the the full-fledged 'write_to_stream()'
        method of the object.
        """
        child_nodes = self._child_nodes
        if child_nodes:
            out.write('(')
            f_child = child_nodes[0]
            for child in child_nodes:
                if child is not f_child:
                    out.write(',')
                child.write_newick(out, **kwargs)
            out.write(')')
        if self.label is not None:
            out.write(self.label)
        if self.edge.length is not None:
            out.write(":{}".format(self.edge.length))

class Edge(object):

    def __init__(self,
            tail_node=None,
            head_node=None,
            length=None,
            rootedge=False,
            label=None):
        self.tail_node = tail_node
        self.head_node = head_node
        self.rootedge = rootedge
        self.length = length

def tree_parsing_fn_factory(src_paths, verbose=False):
    def f():
        for src_path in src_paths:
            if verbose:
                sys.stderr.write("  .. {}\n".format(src_path))
            src = open(src_path, "rU")
            nt = nexusprocessing.NexusTokenizer(src)
            np = newickreader.NewickTreeParser()
            while True:
                t = np.parse_tree_statement(nt, tree_factory=Tree)
                if t is None:
                    break
    return f

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--target-file",
            type=str,
            dest="target_files",
            default=[],
            action="append",
            help="""Path to file to be tokenized; option may be specified multiple times for multiple files. If not specified, default target set will be used.""")
    parser.add_argument("-r", "--repeat",
            type=int,
            default=10,
            help="Repeat each tokenization this number of times (default=%(default)s).")
    parser.add_argument("--delimited_output",
            action="store_true",
            default=False,
            help="Output in tab-delimited instead of aligned format")
    parser.add_argument("--delimited-output",
            action="store_true",
            default=False,
            help="Output in tab-delimited instead of aligned format")
    args = parser.parse_args()

    messenger = messaging.ConsoleMessenger(name="-benchmark")

    src_descs = []
    src_paths = []
    results = []

    if args.target_files:
        for f in args.target_files:
            ff = os.path.expanduser(os.path.expandvars(f))
            src_paths.append(ff)
            src_descs.append( ("User", f) )
    else:
        messenger.info("No sources specified: adding default benchmark target set")
        for f in TREE_FILENAMES:
            ff = pathmap.tree_source_path(f)
            src_paths.append(ff)
            src_descs.append( ("Default", f) )

    for src_path, src_desc in zip(src_paths, src_descs):
        messenger.info("Processing: '{}'".format(src_desc[1]))
        t = timeit.Timer(tree_parsing_fn_factory([src_path]))
        result = min(t.repeat(args.repeat, 1))
        messenger.info("Best time (of {} repetions): {:.10f} seconds".format(args.repeat, result))
        results.append(result)

    messenger.info("Benchmarking complete: all files processed")

    if args.delimited_output:
        result_template = "{}\t{}\t{:.10f}\n"
        header_template = "{}\t{}\t{}\n"
    else:
        max_len1 = max(len(r[0]) for r in src_descs)
        max_len2 = max(len(r[1]) for r in src_descs)
        col1 = "{{:{}}}".format(max_len1)
        col2 = "{{:{}}}".format(max_len2)
        result_template = "[" + col1 + "]  " + col2 + "  {:.10f}\n"
        header_template = col1 + "    " + col2 + "  {}\n"
    sys.stdout.write(header_template.format("Type", "File", "Seconds"))
    for result, src_desc in zip(results, src_descs):
        sys.stdout.write(result_template.format(src_desc[0], src_desc[1], result))

if __name__ == "__main__":
    main()




