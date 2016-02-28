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
Writing of Newick-format tree to a stream.
"""

import re
import warnings
from dendropy.utility import error
from dendropy.utility import textprocessing
from dendropy.dataio import tokenizer
from dendropy.dataio import nexusprocessing
from dendropy.dataio import ioservice

##############################################################################
## NewickReader

class NewickWriter(ioservice.DataWriter):
    """
    Formatter for Newick data.
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
        suppress_rooting : boolean, default: |False|
            If |True|, will not write rooting token ('[&R]' or '[&U]').
            Default is |False|: rooting token will be written.
        suppress_edge_lengths : boolean, default: |False|
            If |True|, will not write edge lengths. Default is |False|: edge
            lengths will be written.
        unquoted_underscores : boolean, default: |False|
            If |True|, labels with underscores will not be quoted, which will
            mean that they will be interpreted as spaces if read again ("soft"
            underscores).  If |False|, then labels with underscores
            will be quoted, resulting in "hard" underscores.  Default is
            |False|.
        preserve_spaces : boolean, default: |False|
            If |True|, spaces will not be replaced with underscores in labels
            (which means any labels containing spaces will have to be quoted).
            Default is |False|: spaces will be converted to underscores.
            False.
        store_tree_weights : boolean, default: |False|
            If |True|, tree weights are written. Default is |False|: tree
            weights will not be written.
        taxon_token_map : boolean or dict or |None|, default: |None|.
            If not |False| or |None|, a "TRANSLATE" statement will be written
            and referenced in tree statements (instead of using the taxon
            labels). If |True|, then a default translate statement will
            be used, with tokens given by the taxon indexes. If a dictionary is
            given, then the keys should be |Taxon| objects and the
            values should be the token (strings).
        suppress_annotations : boolean, default: |True|
            If |False|, metadata annotations will be written out as special
            comments. Defaults to |True|: metadata annotations will be ignored.
        annotations_as_nhx : boolean, default: |False|
            If |True|, and if ``suppress_annotations`` is |False|, will write
            annotation as NHX statements. Default is |False|: annotations
            will not be written as NHX statements.
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

        Notes
        -----

        DendroPy distinguishes between *taxon* labels and *node*
        labels.

        In a Newick string, however, no such distinction is possible, and
        any one node can only be rendered with a single token or symbol. Thus,
        if there is more than one source of label available for a particular
        node (e.g., if both ``suppress_leaf_taxon_labels`` and
        ``suppress_leaf_node_labels`` are |False|, and a particular leaf
        node has both a taxon *and* a label associated with it), then
        the node symbol will be rendered as concatenation of the unsuppressed
        candidate labels, with each candidate label separated by the value
        given in ``node_label_element_separator``. Note that this concatenated
        label requires special handling when being re-read to avoid being
        interpreted as the operational taxonomic unit concept label in
        its entirety. These defaults can all be overridden using the
        various keywords, or a custom label can be composed for the
        node by passing an appropriate function object via the
        ``node_label_compose_fn`` argument.

        Note that, in typical Newick usage, labels of leaf nodes represent
        operational taxonomic unit concepts, and thus the default setting to
        render leaf taxon labels but suppress leaf node labels. Internal node
        labels, on the other hand, are typically used both to represent
        operational taxonomic unit concepts (e.g., ancestral taxa) as well as
        other concepts (e.g., support or geographic range), and thus the
        default internal node rendering is to not to suppress either the taxon
        labels or the node labels.

        """
        ioservice.DataWriter.__init__(self)
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
        self.taxon_token_map = kwargs.pop("taxon_token_map", {})
        self.suppress_annotations = kwargs.pop("suppress_annotations", True)
        # self.suppress_annotations = not kwargs.pop("annotations_as_comments", not self.suppress_annotations) # legacy
        self.annotations_as_nhx = kwargs.pop("annotations_as_nhx", False)
        self.suppress_item_comments = kwargs.pop("suppress_item_comments", True)
        self.suppress_item_comments = not kwargs.pop("write_item_comments", not self.suppress_item_comments)
        self.node_label_element_separator = kwargs.pop("node_label_element_separator", ' ')
        self.node_label_compose_fn = kwargs.pop("node_label_compose_fn", None)
        self.edge_label_compose_fn = kwargs.pop("edge_label_compose_fn", None)
        self._real_value_format_specifier = ""
        self._real_value_formatter = None
        self.real_value_format_specifier = kwargs.pop("real_value_format_specifier", self._real_value_format_specifier)
        if self.edge_label_compose_fn is None:
            self.edge_label_compose_fn = self._format_edge_length
        self.check_for_unused_keyword_arguments(kwargs)

    def _get_taxon_tree_token(self, taxon):
        if self.taxon_token_map is None:
            self.taxon_token_map = {}
        try:
            return self.taxon_token_map[taxon]
        except KeyError:
            t = str(taxon.label)
            self.taxon_token_map[taxon] = t
            return t

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
        Writes a |TreeList| in Newick schema to ``stream``.
        """
        for tree in tree_list:
            self._write_tree(stream, tree)
            stream.write("\n")
        # In Newick format, no clear way to distinguish between
        # annotations/comments associated with tree collection and
        # annotations/comments associated with first tree. So we place them at
        # *end* of document.
        if (not self.suppress_annotations) and (hasattr(tree_list, "_annotations")):
            annotation_comments = nexusprocessing.format_item_annotations_as_comments(tree_list,
                    nhx=self.annotations_as_nhx,
                    real_value_format_specifier=self.real_value_format_specifier,
                    )
        else:
            annotation_comments = ""
        treelist_comments = self._compose_comment_string(tree_list)
        stream.write("{}{}".format(
                annotation_comments,
                treelist_comments))

    def _write_tree(self, stream, tree):
        """
        Composes and writes ``tree`` to ``stream``.
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
            annotation_comments = nexusprocessing.format_item_annotations_as_comments(tree,
                    nhx=self.annotations_as_nhx,
                    real_value_format_specifier=self.real_value_format_specifier,
                    )
        else:
            annotation_comments = ""
        tree_comments = self._compose_comment_string(tree)
        stream.write("{}{}{}{}".format(
                rooting,
                weight,
                annotation_comments,
                tree_comments,
                ))
        tree.apply(
                before_fn=lambda x: self._write_node_open(x, stream),
                after_fn=lambda x: self._write_node_close(x, stream),
                leaf_fn=lambda x: self._write_leaf(x, stream),
                )
        stream.write(";")

    def _write_node_open(self, node, out):
        if node._parent_node is None or node._parent_node._child_nodes[0] is node:
            out.write("(")
        else:
            out.write(",(")

    def _write_leaf(self, node, out):
        if not (node._parent_node is None or node._parent_node._child_nodes[0] is node):
            out.write(",")
        self._write_node_body(node, out)

    def _write_node_close(self, node, out):
        out.write(")")
        self._write_node_body(node, out)

    def _write_node_body(self, node, out):
        out.write(self._render_node_tag(node))
        if node.edge and node.edge.length != None and not self.suppress_edge_lengths:
            out.write(":{}".format(self.edge_label_compose_fn(node.edge)))
        if not self.suppress_annotations:
            node_annotation_comments = nexusprocessing.format_item_annotations_as_comments(node,
                    nhx=self.annotations_as_nhx,
                    real_value_format_specifier=self.real_value_format_specifier)
            out.write(node_annotation_comments)
            edge_annotation_comments = nexusprocessing.format_item_annotations_as_comments(node.edge,
                    nhx=self.annotations_as_nhx,
                    real_value_format_specifier=self.real_value_format_specifier)
            out.write(edge_annotation_comments)
        out.write(self._compose_comment_string(node))
        out.write(self._compose_comment_string(node.edge))

    def _compose_comment_string(self, item):
        if not self.suppress_item_comments and item.comments:
            item_comments = []
            if textprocessing.is_str_type(item.comments):
                item.comments = [item.comments]
            for comment in item.comments:
                item_comments.append("[{}]".format(comment))
            item_comment_str = "".join(item_comments)
        else:
            item_comment_str = ""
        return item_comment_str

    def _render_node_tag(self, node):
        """
        Based on current settings, the attributes of a node, and
        whether or not the node is a leaf, returns an appropriate tag.
        """
        tag = None
        if self.node_label_compose_fn:
            tag = self.node_label_compose_fn(node)
        else:
            tag_parts = []
            is_leaf = len(node.child_nodes()) == 0
            if is_leaf:
                if hasattr(node, 'taxon') \
                        and node.taxon \
                        and node.taxon.label is not None \
                        and not self.suppress_leaf_taxon_labels:
                    tag_parts.append(self._get_taxon_tree_token(node.taxon))
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
                    tag_parts.append(self._get_taxon_tree_token(node.taxon))
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

    # def _compose_node(self, node):
    #     """
    #     Given a DendroPy Node, this returns the Node as a Newick
    #     statement according to the class-defined formatting rules.
    #     """
    #     child_nodes = node.child_nodes()
    #     if child_nodes:
    #         subnodes = [self._compose_node(child) for child in child_nodes]
    #         statement = '(' + ','.join(subnodes) + ')'
    #         if not (self.suppress_internal_taxon_labels and self.suppress_internal_node_labels):
    #             statement = statement + self._render_node_tag(node)
    #         if node.edge and node.edge.length != None and not self.suppress_edge_lengths:
    #             statement =  "{}:{}".format(statement, self.edge_label_compose_fn(node.edge))
    #     else:
    #         statement = self._render_node_tag(node)
    #         if node.edge and node.edge.length != None and not self.suppress_edge_lengths:
    #             statement =  "{}:{}".format(statement, self.edge_label_compose_fn(node.edge))
    #     if not self.suppress_annotations:
    #         node_annotation_comments = nexusprocessing.format_item_annotations_as_comments(node,
    #                 nhx=self.annotations_as_nhx,
    #                 real_value_format_specifier=self.real_value_format_specifier)
    #         edge_annotation_comments = nexusprocessing.format_item_annotations_as_comments(node.edge,
    #                 nhx=self.annotations_as_nhx,
    #                 real_value_format_specifier=self.real_value_format_specifier)
    #         statement = statement + node_annotation_comments + edge_annotation_comments
    #     edge_comment_str = self._compose_comment_string(node.edge)
    #     node_comment_str = self._compose_comment_string(node)
    #     statement = statement + node_comment_str + edge_comment_str
    #     return statement
