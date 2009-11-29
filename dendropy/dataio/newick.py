#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Implementation of NEWICK-format data reader and writer.
"""
from cStringIO import StringIO
import re

from dendropy.utility import containers
from dendropy.utility import texttools
from dendropy.utility import iosys
from dendropy.dataio import nexustokenizer
from dendropy import dataobject

###############################################################################
## lightweight trees from NEWICK sources

def tree_source_iter(stream, **kwargs):
    """
    Iterates over a NEWICK-formatted source of trees given by file-like object
    `stream`

    All keywords defined in `dendropy.dataio.dataformat.tree_source_iter` are
    accepted.

    Note that if `encode_splits` is True, then a `taxon_set` has to be given.
    This is because adding Taxon objects to a taxon set may invalidate split
    bitmasks. Because NEWICK tree taxa are added to a TaxonSet as they are found
    on a tree, there is a strong possibility that all split bitmasks get
    invalidated in the middle of parsing a tree. To avoid this, and, more
    importantly to avoid errors downstream in client code due to this, we
    force specification of a `taxon_set` if `encode_splits` is requested.

    Other optional keyword arguments:

        - `translate_dict` should provide a dictionary mapping taxon numbers (as
           found in the source) to taxon labels (as defined in the source).
        - `rooted` specifies the default rooting interpretation of the tree (see
           `dendropy.dataio.nexustokenizer` for details).
        - `finish_node_func` is a function that will be applied to each node
           after it has been constructed.
        - `edge_len_type` specifies the type of the edge lengths (int or float)
    """
    if "taxon_set" in kwargs:
        taxon_set = kwargs["taxon_set"]
        del(kwargs["taxon_set"])
    else:
        taxon_set = None
    if "encode_splits" in kwargs and taxon_set is None:
        raise Exception('When encoding splits on trees, a pre-populated TaxonSet instance ' \
            + "must be provided using the 'taxon_set' keyword to avoid taxon/split bitmask values "\
            + "changing as new Taxon objects are added to the set.")
    newick_stream = nexustokenizer.NexusTokenizer(stream)
    i = 0
    while not newick_stream.eof:
        t = nexustokenizer.parse_tree_from_stream(newick_stream, taxon_set=taxon_set, **kwargs)
        if t is not None:
            yield t
        else:
            raise StopIteration()

###############################################################################
## split_as_newick_str

def split_as_newick_str(split, taxon_set, preserve_spaces=False):
    """
    Represents a split as a newick string.
    """
    taxlabels = [texttools.escape_nexus_token(label, preserve_spaces=preserve_spaces) for label in taxon_set.labels()]

    # do not do the root
    if split == 0 or (split == taxon_set.all_taxon_set_bitmask()):
        return "(%s)" % (",".join(taxlabels))

    idx = 0
    left = []
    right = []
    while split >= 0 and idx < len(taxlabels):
        if split & 1:
            left.append(taxlabels[idx])
        else:
            right.append(taxlabels[idx])
        idx += 1
        split = split >> 1
    assert ( len(left) + len(right) ) == len(taxlabels)
    return "((%s), (%s))" % (", ".join(left), ", ".join(right))

###############################################################################
## parse_newick_string
#
# def parse_newick_string(tree_statement, taxon_set=None, **kwargs):
#     "Processes a (SINGLE) TREE statement string."
#     stream_handle = StringIO(tree_statement)
#     stream_tokenizer = NexusTokenizer(stream_handle)
#     tree = nexustokenizer.parse_tree_from_stream(stream_tokenizer=stream_tokenizer,
#                                      taxon_set=taxon_set,
#                                      **kwargs)
#     return tree


############################################################################
## CLASS: NewickReader

class NewickReader(iosys.DataReader):
    "Implementation of DataReader for NEWICK files and strings."

    def __init__(self, **kwargs):
        """
        Recognized keywords in addition to those of `DataReader` are:

            - `default_rooting` : default root for trees read in
            - `finish_node_func` : function to be applied to each node on a
                tree as soon as it has been instantiated
            - `allow_duplicate_taxon_labels` : if True, allow duplicate labels
                on trees [False]
        """
        iosys.DataReader.__init__(self, **kwargs)
        self.stream_tokenizer = nexustokenizer.NexusTokenizer()
        self.finish_node_func = kwargs.get("finish_node_func", None)
        self.rooting_interpreter = kwargs.get("rooting_interpreter", nexustokenizer.RootingInterpreter(**kwargs))

    def read(self, stream, **kwargs):
        """
        Instantiates and returns a `DataSet` object based on the
        NEWICK-formatted contents read from the file-like object source
        `stream`.
        """
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        self.rooting_interpreter.update(**kwargs)
        if "taxon_set" in kwargs:
            if self.bound_taxon_set is not None:
                raise TypeError("Cannot specify 'taxon_set' to read() of a Reader with a bound TaxonSet")
            else:
                taxon_set = kwargs["taxon_set"]
        elif self.bound_taxon_set is not None:
            if self.bound_taxon_set not in self.dataset.taxon_sets:
                self.dataset.add_taxon_set(self.bound_taxon_set)
            taxon_set = self.bound_taxon_set
        else:
            taxon_set = self.dataset.new_taxon_set()
        tree_list = self.dataset.new_tree_list(taxon_set=taxon_set)
        kwargs["taxon_set"] = taxon_set
        kwargs["rooting_interpreter"] = self.rooting_interpreter
        for t in tree_source_iter(stream=stream, **kwargs):
            tree_list.append(t, reindex_taxa=False)
        return self.dataset

############################################################################
## CLASS: NewickWriter

class NewickWriter(iosys.DataWriter):
    "Implementation of DataWriter for NEWICK files and strings."

    def __init__(self, **kwargs):
        """
        Recognized keywords in addition to those of `DataWriter` are:

            - `edge_lengths` : if False, edges will not write edge lengths [True]
            - `internal_labels` : if False, internal labels will not be written [True]
        """
        iosys.DataWriter.__init__(self, **kwargs)
        self.edge_lengths = kwargs.get("edge_lengths", True)
        self.internal_labels = kwargs.get("internal_labels", True)
        self.preserve_spaces = kwargs.get("preserve_spaces", False)

    def write(self, stream, **kwargs):
        """
        Writes bound `DataSource` or `TaxonDomain` to a destination given
        by the file-like object `stream`.
        """
        assert self.dataset is not None, \
            "NewickWriter instance is not bound to a DataSet: no source of data"
        self.preserve_spaces = kwargs.get("preserve_spaces", self.preserve_spaces)
        for tree_list in self.dataset.tree_lists:
            if self.bound_taxon_set is None or self.bound_taxon_set is tree_list.taxon_set:
                self.write_tree_list(tree_list, stream)

    def write_tree_list(self, tree_list, stream):
        """
        Writes a `TreeList` in NEWICK format to `stream`.
        """
        if self.exclude_trees:
            return
        for tree in tree_list:
            stream.write(self.compose_node(tree.seed_node) + ';\n')

    def compose_tree(self, tree):
        "Convienience method.        "
        return self.compose_node(tree.seed_node)

    def choose_display_tag(self, node):
        """
        Based on current settings, the attributes of a node, and
        whether or not the node is a leaf, returns an appropriate tag.
        """
        if hasattr(node, 'taxon') and node.taxon:
            tag = node.taxon.label
        elif hasattr(node, 'label') and node.label:
            tag = node.label
        elif len(node.child_nodes()) == 0:
            # force label if a leaf node
            tag = node.oid
        else:
            tag = ""
        if tag:
            tag = texttools.escape_nexus_token(tag, preserve_spaces=self.preserve_spaces)
        return tag

    def compose_node(self, node):
        """
        Given a DendroPy Node, this returns the Node as a NEWICK
        statement according to the class-defined formatting rules.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            subnodes = [self.compose_node(child) for child in child_nodes]
            statement = '(' + ','.join(subnodes) + ')'
            if self.internal_labels:
                statement = statement + self.choose_display_tag(node)
            if node.edge and node.edge.length != None and self.edge_lengths:
                try:
                    statement =  "%s:%f" \
                                % (statement, float(node.edge.length))
                except ValueError:
                    statement =  "%s:%s" \
                                % (statement, node.edge.length)
            return statement
        else:
            statement = self.choose_display_tag(node)
            if node.edge and node.edge.length != None and self.edge_lengths:
                try:
                    statement =  "%s:%0.10f" \
                                % (statement, float(node.edge.length))
                except ValueError:
                    statement =  "%s:%s" \
                                % (statement, node.edge.length)
            return statement
