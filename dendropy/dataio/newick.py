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

    The following keyword arguments are recognized:

        - `taxon_set` specifies the `TaxonSet` object to be attached to the
           trees parsed and manage their taxa. If not specified, then a
           (single) new `TaxonSet` object will be created and for all the
           `Tree` objects.
        - `encode_splits` specifies whether or not split bitmasks will be
           calculated and attached to the edges.
        - `from_index` 0-based index specifying first tree to actually return

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
    if "from_index" in kwargs:
        from_index = kwargs.get("from_index")
        del(kwargs["from_index"])
    else:
        from_index = 0
    if "encode_splits" in kwargs and taxon_set is None:
        raise Exception('When encoding splits on trees, a pre-populated TaxonSet instance ' \
            + "must be provided using the 'taxon_set' keyword to avoid taxon/split bitmask values "\
            + "changing as new Taxon objects are added to the set.")
    newick_stream = nexustokenizer.NexusTokenizer(stream)
    i = 0
    while not newick_stream.eof:
        t = nexustokenizer.parse_tree_from_stream(newick_stream, taxon_set=taxon_set, **kwargs)
        if t is not None:
            if i >= from_index:
                yield t
            i += 1
        else:
            raise StopIteration()

def write_tree_list(tree_list, stream, **kwargs):
    """
    Writes out a list of trees in NEWICK-format to a destination given by
    file-like object `stream`.

    Additionally, the following keywords are recognized:

        - `edge_lengths` : if False, edges will not write edge lengths [True]
        - `internal_labels` : if False, internal labels will not be written [True]
    """
    newick_writer = NewickWriter(edge_lengths=kwargs.get("edge_lengths", True),
                                 internal_labels=kwargs.get("internal_labels", True))
    newick_writer.write_tree_list(tree_list, stream)

def read_tree_list(stream, **kwargs):
    """
    Parses a source describing a collection of trees in NEWICK format, and
    returns corresponding `TreeList` object..

    Additionally, a `TreeList` object to which to add the trees may be passed
    by `tree_list`, or, alternatively, a `TaxonSet` object with which to manage
    the taxa by `taxon_set`. Only one of `tree_list` or `taxon_set` may be
    specified. If neither is specified, then a new `TreeList` object, with its
    own associated new `TaxonSet` object will be created.
    """
    if "tree_list" in kwargs:
        assert "taxon_set" not in kwargs, \
            "Cannot specify both 'tree_list' and 'taxon_set'."
        tree_list = kwargs["tree_list"]
        del(kwargs["tree_list"])
    elif "taxon_set" in kwargs:
        tree_list = dataobject.TreeList(taxon_set=kwargs["taxon_set"])
        del(kwargs["taxon_set"])
    else:
        tree_list = dataobject.TreeList()
    for t in tree_source_iter(stream=stream, taxon_set=tree_list.taxon_set, **kwargs):
        if t is not None:
            tree_list.append(t, reindex_taxa=False)
    return tree_list

###############################################################################
## split_as_newick_str

def split_as_newick_str(split, taxon_set):
    """
    Represents a split as a newick string.
    """
    taxlabels = [texttools.escape_nexus_token(label) for label in taxon_set.labels()]

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
        self.default_rooting = kwargs.get("default_rooting", \
                nexustokenizer.RootingInterpretation.UNKNOWN_DEF_ROOTED)
        self.finish_node_func = kwargs.get("finish_node_func", None)

    def read(self, stream, **kwargs):
        """
        Instantiates and returns a `DataSet` object based on the
        NEWICK-formatted contents read from the file-like object source
        `stream`.
        """
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        if self.bound_taxon_set is not None:
            if self.bound_taxon_set not in self.dataset.taxon_sets:
                self.dataset.add_taxon_set(self.bound_taxon_set)
            taxon_set = self.bound_taxon_set
        else:
            taxon_set = self.dataset.new_taxon_set()
        tree_list = self.dataset.new_tree_list(taxon_set=taxon_set)
        if "rooted" not in kwargs:
            kwargs["rooted"] = self.default_rooting
        read_tree_list(stream=stream,
                       tree_list=tree_list,
                       **kwargs)
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

    def write(self, stream, **kwargs):
        """
        Writes bound `DataSource` or `TaxonDomain` to a destination given
        by the file-like object `stream`.
        """
        assert self.dataset is not None, \
            "NewickWriter instance is not bound to a DataSet: no source of data"
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
            return texttools.escape_nexus_token(node.taxon.label)
        elif hasattr(node, 'label') and node.label:
            return texttools.escape_nexus_token(node.label)
        elif len(node.child_nodes()) == 0:
            # force label if a leaf node
            return texttools.escape_nexus_token(node.oid)
        else:
            return ""

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
            if self.internal_labels:
                statement = self.choose_display_tag(node)
            if node.edge and node.edge.length != None and self.edge_lengths:
                try:
                    statement =  "%s:%0.10f" \
                                % (statement, float(node.edge.length))
                except ValueError:
                    statement =  "%s:%s" \
                                % (statement, node.edge.length)
            return statement
