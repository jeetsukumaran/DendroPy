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
Writing of D3-format tree to a stream.
"""

import re
import warnings
from dendropy.utility import error
from dendropy.utility import textprocessing
from dendropy.dataio import ioservice

##############################################################################
## D3Writer

class D3Writer(ioservice.DataWriter):
    """
    Formatter for D3 data.
    """

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------
        suppress_leaf_taxon_labels : boolean, default: |False|
            If |True|, then taxon labels will not be rendered for leaves.
            Default is |False|: render leaf taxon labels. See notes below for
            details.
        suppress_leaf_node_labels : boolean, default: |True|
            If |False|, then node labels (if available) will be printed for
            leaves. Defaults to |True|: do not render leaf node labels. See
            notes below for details.
        suppress_internal_taxon_labels : boolean, default: |False|
            If |True|, then taxon labels will not be printed for internal
            nodes. Default is |False|: print taxon labels for internal nodes.
            See notes below for details.
        suppress_internal_node_labels : boolean, default: |False|
            If |True|, then node labels will not be printed for internal nodes.
            Default is |False|: print node labels for internal nodes. See notes
            below for details.
        suppress_edge_lengths : boolean, default: |False|
            If |True|, will not write edge lengths. Default is |False|: edge
            lengths will be written.
        store_tree_weights : boolean, default: |False|
            If |True|, tree weights are written. Default is |False|: tree
            weights will not be written.
        suppress_annotations : boolean, default: |True|
            If |False|, metadata annotations will be written out as special
            comments. Defaults to |True|: metadata annotations will be ignored.
        suppress_item_comments : boolean, default: |True|
            If |False|, any additional comments associated with trees, nodes,
            edges, etc. will be written. Default is |True|: comments will be
            ignored.
        node_label_element_separator : string, default: ' '
            If both ``suppress_leaf_taxon_labels`` and
            ``suppress_leaf_node_labels`` are |False|, then this will be the
            string used to join them. Defaults to ' ' (space).
        node_label_compose_fn : function object or |None|, default: |None|
            If not |None|, should be a function that takes a |Node|
            object as an argument and returns the string to be used to
            represent the node in the tree statement. The return value from
            this function is used unconditionally to print a node
            representation in a tree statement, by-passing the default
            labelling function, ignoring ``suppress_leaf_taxon_labels``,
            ``suppress_leaf_node_labels=True``, ``suppress_internal_taxon_labels``,
            ``suppress_internal_node_labels``, etc. Defaults to |None|.
        edge_label_compose_fn : function object or |None|, default: |None|
            If not |None|, should be a function that takes an Edge object as
            an argument, and returns the string to be used to represent the
            edge length in the tree statement.
        real_value_format_specifier : string, default: ''
            Format specification for real/float values. Will be applied to edge
            lengths (if ``edge_label_compose_fn`` is not given) as well as
            annotations. The format specifier should be given in Python's
            string format specification mini-language. E.g. ".8f", ".4E",
            "8.4f".
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.
        """
        ioservice.DataWriter.__init__(self)
        ioservice.DataWriter.__init__(self, **kwargs)
        self.suppress_leaf_taxon_labels = kwargs.pop("suppress_leaf_taxon_labels", False)
        self.suppress_leaf_node_labels = kwargs.pop("suppress_leaf_node_labels", True)
        self.suppress_internal_taxon_labels = kwargs.pop("suppress_internal_taxon_labels", False)
        self.suppress_internal_node_labels = kwargs.pop("suppress_internal_node_labels", False)
        self.suppress_edge_lengths = kwargs.pop("suppress_edge_lengths", False)
        self.preserve_spaces = kwargs.pop("preserve_spaces", False)
        self.store_tree_weights = kwargs.pop("store_tree_weights", False)
        self.suppress_annotations = kwargs.pop("suppress_annotations", True)
        self.suppress_item_comments = kwargs.pop("suppress_item_comments", True)
        self.node_label_element_separator = kwargs.pop("node_label_element_separator", ' ')
        self.node_label_compose_fn = kwargs.pop("node_label_compose_fn", None)
        self.edge_label_compose_fn = kwargs.pop("edge_label_compose_fn", None)
        self._real_value_format_specifier = ""
        self._real_value_formatter = None
        self.real_value_format_specifier = kwargs.pop("real_value_format_specifier", self._real_value_format_specifier)
        if self.edge_label_compose_fn is None:
            self.edge_label_compose_fn = self._format_edge_length
        self.check_for_unused_keyword_arguments(kwargs)

    def _get_real_value_format_specifier(self):
        return self._real_value_format_specifier
    def _set_real_value_format_specifier(self, f):
        if f is None:
            f = ""
        self._real_value_format_specifier = f
        s = "{:" + self._real_value_format_specifier + "}"
        self._real_value_formatter = s.format
    real_value_format_specifier = property(_get_real_value_format_specifier, _set_real_value_format_specifier)

    def _format_edge_length(self, edge):
        """
        Note: instance method to allow overriding.
        """
        return self._real_value_formatter(edge.length)

    def _write(self,
            stream,
            taxon_namespaces=None,
            tree_lists=None,
            char_matrices=None,
            global_annotations_target=None):
        for tree_list in tree_lists:
            if (self.attached_taxon_namespace is not None
                    and tree_list.taxon_namespace is not self.attached_taxon_namespace):
                continue
            self._write_tree_list(stream, tree_list)

    def _write_tree_list(self, stream, tree_list):
        """
        Writes a |TreeList| in D3 schema to ``stream``.
        """
        for tree in tree_list:
            self._write_tree(stream, tree)
            stream.write("\n")

    def _write_tree(self, stream, tree):
        """
        Composes and writes ``tree`` to ``stream``.
        """
        pass

