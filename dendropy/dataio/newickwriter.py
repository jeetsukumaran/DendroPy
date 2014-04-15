#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Writing of NEWICK-format tree to a stream.
"""

import re
import warnings
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.datamodel import base
from dendropy.utility import error
from dendropy.dataio import tokenizer
from dendropy.dataio import nexusprocessing
from dendropy.dataio import ioservice

##############################################################################
## NewickReader

class NewickWriter(ioservice.DataWriter):
    """
    Formatter for NEWICK data.
    """

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------
        suppress_leaf_taxon_labels : boolean, default: `False`
            If `True`, then taxon labels will not be printed for leaves.
            Default is `False`:  print leaf taxon labels.  If both
            `suppress_leaf_taxon_labels` and `suppress_leaf_node_labels` are
            `False`, and a particular leaf node has both a label and is
            associated with taxa, then the final label for that node will be
            the concatentation of the two labels, separated by the value given
            in `nodel_label_element_separator`
        suppress_leaf_node_labels : boolean, default: `True`
            If `False`, then node labels (if available) will be printed for
            leaves. Defaults to `True`: do not print leaf node labels. Note
            that DendroPy distinguishes between *taxon* labels and *node*
            labels. In a typical NEWICK string, taxon labels are printed for
            leaf nodes, while leaf node labels are ignored (hence the default
            `True` setting to suppress leaf node labels).  If both
            `suppress_leaf_taxon_labels` and `suppress_leaf_node_labels` are
            `False`, and a particular leaf node has both a label and is
            associated with taxa, then the final label for that node will be
            the concatentation of the two labels, separated by the value given
            in `nodel_label_element_separator`
        suppress_internal_taxon_labels : boolean, default: `False`
            If `True`, then taxon labels will not be printed for internal
            nodes. Default is `False`: print taxon labels for internal nodes.
            If both `suppress_leaf_taxon_labels` and
            `suppress_leaf_node_labels` are `False`, and a particular leaf node
            has both a label and is associated with taxa, then the final label
            for that node will be the concatentation of the two labels,
            separated by the value given in `nodel_label_element_separator`
        suppress_internal_node_labels : boolean, default: `False`
            If `True`, then node labels will not be printed for internal nodes.
            Default is `False`: print node labels for internal nodes.  If both
            `suppress_leaf_taxon_labels` and `suppress_leaf_node_labels` are
            `False`, and a particular leaf node has both a label and is
            associated with taxa, then the final label for that node will be
            the concatentation of the two labels, separated by the value given
            in `nodel_label_element_separator`
        suppress_rooting : boolean, default: `False`
            If `True`, will not write rooting token ('[&R]' or '[&U]').
            Default is `False`: rooting token will be written.
        suppress_edge_lengths : boolean, default: `False`
            If `True`, will not write edge lengths. Default is `False`: edge
            lengths will be written.
        unquoted_underscores : boolean, default: `False`
            If `True`, labels with underscores will not be quoted, which will
            mean that they will be interpreted as spaces if read again ("soft"
            underscores).  If `False`, then labels with underscores
            will be quoted, resulting in "hard" underscores.  Default is
            `False`.
        preserve_spaces : boolean, default: `False`
            If `True`, spaces will not be replaced with underscores in labels
            (which means any labels containing spaces will have to be quoted).
            Default is `False`: spaces will be converted to underscores.
            False.
        store_tree_weights : boolean, default: `False`
            If `True`, tree weights are written. Default is `False`: tree
            weights will not be written.
        supppress_annotations : boolean, default: `True`
            If `False`, metadata annotations will be written out as special
            comments. Defaults to `True`: metadata annotations will be ignored.
        annotations_as_nhx : boolean, default: `False`
            If `True`, and if `suppress_annotations` is `False`, will write
            annotation as NHX statements. Default is `False`: annotations
            will not be written as NHX statements.
        suppress_item_comments : boolean, default: `True`
            If `False`, any additional comments associated with trees, nodes,
            edges, etc. will be written. Default is `True`: comments will be
            ignored.
        node_label_element_separator : string, default: ' '
            If both `suppress_leaf_taxon_labels` and
            `suppress_leaf_node_labels` are `False`, then this will be the
            string used to join them. Defaults to ' ' (space).
        node_label_compose_func : function object or `None`, default: `None`
            If not `None`, should be a function that takes a :class:`Node`
            object as an argument and returns the string to be used to
            represent the node in the tree statement. The return value from
            this function is used unconditionally to print a node
            representation in a tree statement, by-passing the default
            labelling function, ignoring `suppress_leaf_taxon_labels`,
            `suppress_leaf_node_labels=True`, `suppress_internal_taxon_labels`,
            `suppress_internal_node_labels`, etc. Defaults to `None`.
        edge_label_compose_func : function object or `None`, default: `None`
            If not `None`, should be a function that takes an Edge object as
            an argument, and returns the string to be used to represent the
            edge length in the tree statement.

        Typically, these keywords would be passed to the `write_to_path()`,
        `write_to_stream` or `as_string` arguments, when 'newick' is used as
        the schema::

            d.write_to_path('outputfile.tre',
                    schema='newick',
                    suppress_leaf_taxon_labels=False,
                    suppress_leaf_node_labels=True,
                    suppress_internal_taxon_labels=False,
                    suppress_internal_node_labels=False,
                    suppress_rooting=False,
                    suppress_edge_lengths=False,
                    unquoted_underscores=False,
                    preserve_spaces=False,
                    store_tree_weights=False,
                    suppress_annotations=True,
                    annotations_as_nhx=False,
                    suppress_item_comments=True,
                    node_label_element_separator=' ',
                    node_label_compose_func=None,
                    edge_label_compose_func=None)

        """
        legacy = {
                "internal_labels": "Use 'suppress_internal_taxon_labels' instead",
                "write_rooting": "Use 'suppress_rooting' instead",
                "quote_undescores" : "Use 'unquoted_underscores' instead",
                "annotations_as_comments": "Use 'suppress_annotations' instead"
                }
        for kw in legacy:
            if kw in kwargs:
                raise TypeError("'{}' is no longer supported: {}".format(kw, legacy[kw]))
        ioservice.DataWriter.__init__(self, **kwargs)
        self.suppress_leaf_taxon_labels = kwargs.pop("suppress_leaf_taxon_labels", False)
        self.suppress_leaf_node_labels = kwargs.pop("suppress_leaf_node_labels", True)
        self.suppress_internal_taxon_labels = kwargs.pop("suppress_internal_taxon_labels", False)
        # self.suppress_internal_taxon_labels = not kwargs.pop("internal_labels", not sef.suppress_internal_taxon_labels) # legacy
        self.suppress_internal_node_labels = kwargs.pop("suppress_internal_node_labels", False)
        self.suppress_internal_node_labels = not kwargs.pop("internal_labels", not self.suppress_internal_node_labels) # legacy
        self.suppress_rooting = kwargs.pop("suppress_rooting", False)
        # self.suppress_rooting = not kwargs.pop("write_rooting", not self.suppress_rooting) # legacy
        self.suppress_edge_lengths = kwargs.pop("suppress_edge_lengths", False)
        # self.suppress_edge_lengths = not kwargs.pop("edge_lengths", not self.suppress_edge_lengths) # legacy
        self.unquoted_underscores = kwargs.pop("unquoted_underscores", False)
        # self.unquoted_underscores = not kwargs.pop('quote_underscores', not self.unquoted_underscores) # legacy
        self.preserve_spaces = kwargs.pop("preserve_spaces", False)
        self.store_tree_weights = kwargs.pop("store_tree_weights", False)
        self.suppress_annotations = kwargs.pop("suppress_annotations", True)
        # self.suppress_annotations = not kwargs.pop("annotations_as_comments", not self.suppress_annotations) # legacy
        self.annotations_as_nhx = kwargs.pop("annotations_as_nhx", False)
        self.suppress_item_comments = kwargs.pop("suppress_item_comments", True)
        self.suppress_item_comments = not kwargs.pop("write_item_comments", not self.suppress_item_comments)
        self.node_label_element_separator = kwargs.pop("node_label_element_separator", ' ')
        self.node_label_compose_func = kwargs.pop("node_label_compose_func", None)
        self.edge_label_compose_func = kwargs.pop("edge_label_compose_func", None)
        if self.edge_label_compose_func is None:
            self.edge_label_compose_func = self._format_edge_length

    def _format_edge_length(self, edge):
        """
        Note: instance method to allow overriding.
        """
        return "{}".format(edge.length)

    def _write(self,
            stream,
            taxon_namespaces=None,
            tree_lists=None,
            char_matrices=None,
            global_annotations_target=None):
        for tree_list in tree_lists:
            self._write_tree_list(tree_list, stream)

    def _write_tree_list(self, tree_list, stream):
        """
        Writes a `TreeList` in NEWICK schema to `stream`.
        """
        for tree in tree_list:
            self._write_tree(tree, stream)
            stream.write("\n")

    def _write_tree(self, tree, stream):
        """
        Composes and writes `tree` to `stream`.
        """
        if tree.rooting_state_is_undefined or self.suppress_rooting:
            rooting = ""
        elif tree.is_rooted:
            rooting = "[&R] "
        elif not tree.is_rooted:
            rooting = "[&U] "
        else:
            rooting = ""
        if self.store_tree_weights and tree.weight is not None:
            weight = "[&W {}] ".format(tree.weight)
        else:
            weight = ""
        if not self.suppress_annotations or self.annotations_as_nhx:
            annotation_comments = nexusprocessing.format_annotation_as_comments(tree, nhx=self.annotations_as_nhx)
        else:
            annotation_comments = ""
        tree_comments = self._compose_comment_string(tree)
        stream.write("{}{}{}{}{};".format(
                rooting,
                weight,
                annotation_comments,
                tree_comments,
                self._compose_node(tree.seed_node)))

    def _compose_comment_string(self, item):
        if not self.suppress_item_comments and item.comments:
            item_comments = []
            if isinstance(item.comments, str):
                item.comments = [item.comments]
            for comment in item.comments:
                item_comments.append("[{}]".format(comment))
            item_comment_str = "".join(item_comments)
        else:
            item_comment_str = ""
        return item_comment_str

    def _choose_display_tag(self, node):
        """
        Based on current settings, the attributes of a node, and
        whether or not the node is a leaf, returns an appropriate tag.
        """
        tag = None
        if self.node_label_compose_func:
            tag = self.node_label_compose_func(node)
        else:
            tag_parts = []
            is_leaf = len(node.child_nodes()) == 0
            if is_leaf:
                if hasattr(node, 'taxon') \
                        and node.taxon \
                        and node.taxon.label is not None \
                        and not self.suppress_leaf_taxon_labels:
                    tag_parts.append(str(node.taxon.label))
                if hasattr(node, 'label') \
                        and node.label \
                        and node.label is not None \
                        and not self.suppress_leaf_node_labels:
                    tag_parts.append(str(node.label))
                if len(tag_parts) > 0:
                    tag = self.node_label_element_separator.join(tag_parts)
                else:
                    return "_" # anonymous leaf
            else:
                if hasattr(node, 'taxon') \
                        and node.taxon \
                        and node.taxon.label is not None \
                        and not self.suppress_internal_taxon_labels:
                    tag_parts.append(str(node.taxon.label))
                if hasattr(node, 'label') \
                        and node.label \
                        and node.label is not None \
                        and not self.suppress_internal_node_labels:
                    tag_parts.append(str(node.label))
                if len(tag_parts) > 0:
                    tag = self.node_label_element_separator.join(tag_parts)
                else:
                    return "" # nada
        if tag:
            tag = nexusprocessing.escape_nexus_token(tag,
                    preserve_spaces=self.preserve_spaces,
                    quote_underscores=not self.unquoted_underscores)
            return tag
        else:
            return ""

    def _compose_node(self, node):
        """
        Given a DendroPy Node, this returns the Node as a NEWICK
        statement according to the class-defined formatting rules.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            subnodes = [self._compose_node(child) for child in child_nodes]
            statement = '(' + ','.join(subnodes) + ')'
            if not (self.suppress_internal_taxon_labels and self.suppress_internal_node_labels):
                statement = statement + self._choose_display_tag(node)
            if node.edge and node.edge.length != None and not self.suppress_edge_lengths:
                statement =  "{}:{}".format(statement, self.edge_label_compose_func(node.edge))
        else:
            statement = self._choose_display_tag(node)
            if node.edge and node.edge.length != None and not self.suppress_edge_lengths:
                statement =  "{}:{}".format(statement, self.edge_label_compose_func(node.edge))
        if not self.suppress_annotations or self.annotations_as_nhx:
            node_annotation_comments = nexusprocessing.format_annotation_as_comments(node, nhx=self.annotations_as_nhx)
            edge_annotation_comments = nexusprocessing.format_annotation_as_comments(node.edge, nhx=self.annotations_as_nhx)
            statement = statement + node_annotation_comments + edge_annotation_comments
        #if self.nhx_key_to_func:
        #    nhx_to_print = []
        #    for k, v in self.nhx_key_to_func.items():
        #        r = v(node.edge)
        #        if r is not None:
        #            nhx_to_print.append("%s=%s" % (k, str(r)))
        #    if nhx_to_print:
        #        statement = statement + ('[&&NHX:%s]' % ':'.join(nhx_to_print))
        node_comment_str = self._compose_comment_string(node)
        statement += node_comment_str
        return statement