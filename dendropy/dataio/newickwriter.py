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
            If `True`, then taxon labels will not be rendered for leaves.
            Default is `False`: render leaf taxon labels. See notes below for
            details.
        suppress_leaf_node_labels : boolean, default: `True`
            If `False`, then node labels (if available) will be printed for
            leaves. Defaults to `True`: do not render leaf node labels. See
            notes below for details.
        suppress_internal_taxon_labels : boolean, default: `False`
            If `True`, then taxon labels will not be printed for internal
            nodes. Default is `False`: print taxon labels for internal nodes.
            See notes below for details.
        suppress_internal_node_labels : boolean, default: `False`
            If `True`, then node labels will not be printed for internal nodes.
            Default is `False`: print node labels for internal nodes. See notes
            below for details.
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
        Notes
        -----

        DendroPy distinguishes between *taxon* labels and *node*
        labels.

        In a NEWICK string, however, no such distinction is possible, and
        any one node can only be rendered with a single token or symbol. Thus,
        if there is more than one source of label available for a particular
        node (e.g., if both `suppress_leaf_taxon_labels` and
        `suppress_leaf_node_labels` are `False`, and a particular leaf
        node has both a taxon *and* a label associated with it), then
        the node symbol will be rendered as concatenation of the unsuppressed
        candidate labels, with each candidate label separated by the value
        given in `node_label_element_separator`. Note that this concatenated
        label requires special handling when being re-read to avoid being
        interpreted as the operational taxonomic unit concept label in
        its entirety. These defaults can all be overridden using the
        various keywords, or a custom label can be composed for the
        node by passing an appropriate function object via the
        `node_label_compose_func` argument.

        Note that, in typical NEWICK usage, labels of leaf nodes represent
        operational taxonomic unit concepts, and thus the default setting to
        render leaf taxon labels but suppress leaf node labels. Internal node
        labels, on the other hand, are typically used both to represent
        operational taxonomic unit concepts (e.g., ancestral taxa) as well as
        other concepts (e.g., support or geographic range), and thus the
        default internal node rendering is to not to suppress either the taxon
        labels or the node labels.

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
        # self.suppress_internal_node_labels = not kwargs.pop("internal_labels", not self.suppress_internal_node_labels) # legacy
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
        self.check_for_unused_keyword_arguments(kwargs)

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
            self._write_tree_list(stream, tree_list)

    def _write_tree_list(self, stream, tree_list):
        """
        Writes a `TreeList` in NEWICK schema to `stream`.
        """
        for tree in tree_list:
            self._write_tree(stream, tree)
            stream.write("\n")
        # In NEWICK format, no clear way to distinguish between
        # annotations/comments associated with tree collection and
        # annotations/comments associated with first tree. So we place them at
        # *end* of document.
        if (not self.suppress_annotations) and (hasattr(tree_list, "_annotations")):
            annotation_comments = nexusprocessing.format_item_annotations_as_comments(tree_list, nhx=self.annotations_as_nhx)
        else:
            annotation_comments = ""
        treelist_comments = self._compose_comment_string(tree_list)
        stream.write("{}{}".format(
                annotation_comments,
                treelist_comments))

    def _write_tree(self, stream, tree):
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
        if not self.suppress_annotations:
            annotation_comments = nexusprocessing.format_item_annotations_as_comments(tree, nhx=self.annotations_as_nhx)
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
                    return "" # anonymous leaf
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
                    return ""
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
        if not self.suppress_annotations:
            node_annotation_comments = nexusprocessing.format_item_annotations_as_comments(node, nhx=self.annotations_as_nhx)
            edge_annotation_comments = nexusprocessing.format_item_annotations_as_comments(node.edge, nhx=self.annotations_as_nhx)
            statement = statement + node_annotation_comments + edge_annotation_comments
        edge_comment_str = self._compose_comment_string(node.edge)
        node_comment_str = self._compose_comment_string(node)
        statement = statement + node_comment_str + edge_comment_str
        return statement
