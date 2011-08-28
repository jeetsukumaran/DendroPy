#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
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
Implementation of NEWICK-schema data reader and writer.
"""
from cStringIO import StringIO
import re

from dendropy.utility import containers
from dendropy.utility import textutils
from dendropy.utility import iosys
from dendropy.dataio import nexustokenizer
from dendropy import dataobject

###############################################################################
## lightweight trees from NEWICK sources

def tree_source_iter(stream, **kwargs):
    """
    Iterates over a NEWICK-formatted source of trees given by file-like object
    `stream`

    Note that if `encode_splits` is True, then a `taxon_set` has to be given.
    This is because adding Taxon objects to a taxon set may invalidate split
    bitmasks. Because NEWICK tree taxa are added to a TaxonSet as they are found
    on a tree, there is a strong possibility that all split bitmasks get
    invalidated in the middle of parsing a tree. To avoid this, and, more
    importantly to avoid errors downstream in client code due to this, we
    force specification of a `taxon_set` if `encode_splits` is requested.

    The following optional keyword arguments are also recognized:

        `taxon_set`
            TaxonSet object to use when reading data.

        `as_rooted=True` (or `as_unrooted=False`)
            Unconditionally interprets all trees as rooted.

        `as_unrooted=True` (or `as_rooted=False`)
            Unconditionally interprets all trees as unrooted.

        `default_as_rooted=True` (or `default_as_unrooted=False`)
            Interprets all trees as rooted if rooting not given by `[&R]`
            or `[&U]` comments.

        `default_as_unrooted=True` (or `default_as_rooted=False`)
            Interprets all trees as rooted if rooting not given by `[&R]`
            or `[&U]` comments.

        `edge_len_type`
            Specifies the type of the edge lengths (int or float).

        `extract_comment_metadata`
            If True, any 'hot comments', i.e., comments that begin with
            '&', or NHX comments associated with items will be processed
            and stored as a dictionary attribute of the object:
            "comment_metadata".

        `store_tree_weights`
            If True, process the tree weight ("[&W 1/2]") comment
            associated with each tree, if any.

        `encode_splits`
            Specifies whether or not split bitmasks will be calculated and
            attached to the edges.

        `finish_node_func`
            Is a function that will be applied to each node after it has
            been constructed.

        `case_sensitive_taxon_labels`
            If True, then taxon labels are case sensitive (different cases
            = different taxa); defaults to False.

        `allow_duplicate_taxon_labels`
            if True, allow duplicate labels on trees

        `preserve_underscores`
            If True, unquoted underscores in labels will *not* converted to
            spaces. Defaults to False: all underscores not protected by
            quotes will be converted to spaces.

        `suppress_internal_node_taxa`
            If True, internal node labels will not be treated as taxa.
            Defaults to False: internal node labels will be instantantiatd
            into Taxon objects.

        `allow_duplicate_taxon_labels`
            If True, then multiple identical taxon labels will be allowed.
            Defaults to False: treat multiple identical taxon labels as an
            error.

        `hyphens_as_tokens`
            If True, hyphens will be treated as special punctuation
            characters. Defaults to False, hyphens not treated as special
            punctuation characters.

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
    preserve_underscores = kwargs.get('preserve_underscores', False)
    hyphens_as_tokens = kwargs.get('hyphens_as_tokens', nexustokenizer.DEFAULT_HYPHENS_AS_TOKENS)
    extract_comment_metadata = kwargs.get("extract_comment_metadata", False)
    newick_stream = nexustokenizer.NexusTokenizer(stream,
                                                  preserve_underscores=preserve_underscores,
                                                  hyphens_as_tokens=hyphens_as_tokens,
                                                  extract_comment_metadata=extract_comment_metadata,
                                                  case_sensitive_taxon_labels=kwargs.get('case_sensitive_taxon_labels', False))
    while not newick_stream.eof:
        t = nexustokenizer.tree_from_token_stream(newick_stream, taxon_set=taxon_set, **kwargs)
        if t is not None:
            yield t
        else:
            raise StopIteration()

###############################################################################
## parse_newick_string
#
# def parse_newick_string(tree_statement, taxon_set=None, **kwargs):
#     "Processes a (SINGLE) TREE statement string."
#     stream_handle = StringIO(tree_statement)
#     stream_tokenizer = NexusTokenizer(stream_handle)
#     tree = nexustokenizer.tree_from_token_stream(stream_tokenizer=stream_tokenizer,
#                                      taxon_set=taxon_set,
#                                      **kwargs)
#     return tree


############################################################################
## CLASS: NewickReader

class NewickReader(iosys.DataReader):
    "Implementation of DataReader for NEWICK files and strings."

    def __init__(self, **kwargs):
        """
        __init__ recognizes the following keywords arguments:

            `taxon_set`
                TaxonSet object to use when reading data.

            `as_rooted=True` (or `as_unrooted=False`)
                Unconditionally interprets all trees as rooted.

            `as_unrooted=True` (or `as_rooted=False`)
                Unconditionally interprets all trees as unrooted.

            `default_as_rooted=True` (or `default_as_unrooted=False`)
                Interprets all trees as rooted if rooting not given by `[&R]`
                or `[&U]` comments.

            `default_as_unrooted=True` (or `default_as_rooted=False`)
                Interprets all trees as rooted if rooting not given by `[&R]`
                or `[&U]` comments.

            `edge_len_type`
                Specifies the type of the edge lengths (int or float).

            `extract_comment_metadata`
                If True, any 'hot comments', i.e., comments that begin with
                '&', or NHX comments associated with items will be processed
                and stored as a dictionary attribute of the object:
                "comment_metadata".

            `store_tree_weights`
                If True, process the tree weight ("[&W 1/2]") comment
                associated with each tree, if any.

            `encode_splits`
                Specifies whether or not split bitmasks will be calculated and
                attached to the edges.

            `finish_node_func`
                Is a function that will be applied to each node after it has
                been constructed.

            `case_sensitive_taxon_labels`
                If True, then taxon labels are case sensitive (different cases
                = different taxa); defaults to False.

            `allow_duplicate_taxon_labels`
                if True, allow duplicate labels on trees

            `preserve_underscores`
                If True, unquoted underscores in labels will *not* converted to
                spaces. Defaults to False: all underscores not protected by
                quotes will be converted to spaces.

            `suppress_internal_node_taxa`
                If True, internal node labels will not be treated as taxa.
                Defaults to False: internal node labels will be instantantiatd
                into Taxon objects.

            `allow_duplicate_taxon_labels`
                If True, then multiple identical taxon labels will be allowed.
                Defaults to False: treat multiple identical taxon labels as an
                error.

            `hyphens_as_tokens`
                If True, hyphens will be treated as special punctuation
                characters. Defaults to False, hyphens not treated as special
                punctuation characters.
        """
        iosys.DataReader.__init__(self, **kwargs)
        self.finish_node_func = kwargs.get("finish_node_func", None)
        self.rooting_interpreter = kwargs.get("rooting_interpreter", nexustokenizer.RootingInterpreter(**kwargs))
        self.hyphens_as_tokens = kwargs.get('hyphens_as_tokens', nexustokenizer.DEFAULT_HYPHENS_AS_TOKENS)
        self.encode_splits = kwargs.get('encode_splits', False)
        self.extract_comment_metadata = kwargs.get('extract_comment_metadata', False)
        self.store_tree_weights = kwargs.get('store_tree_weights', False)
        self.preserve_underscores = kwargs.get('preserve_underscores', False)
        self.suppress_internal_node_taxa = kwargs.get("suppress_internal_node_taxa", False)
        self.case_sensitive_taxon_labels = kwargs.get('case_sensitive_taxon_labels', False)
        self.edge_len_type = kwargs.get('edge_len_type', float)

    def read(self, stream):
        """
        Instantiates and returns a `DataSet` object based on the
        NEWICK-formatted contents read from the file-like object source
        `stream`.
        """
        if self.exclude_trees:
            return self.dataset
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        taxon_set = self.get_default_taxon_set()
        tree_list = self.dataset.new_tree_list(taxon_set=taxon_set)
        for t in tree_source_iter(stream=stream,
                taxon_set=taxon_set,
                rooting_interpreter=self.rooting_interpreter,
                hyphens_as_tokens=self.hyphens_as_tokens,
                extract_comment_metadata=self.extract_comment_metadata,
                store_tree_weights=self.store_tree_weights,
                encode_splits=self.encode_splits,
                preserve_underscores=self.preserve_underscores,
                suppress_internal_node_taxa=self.suppress_internal_node_taxa,
                edge_len_type=self.edge_len_type,
                case_sensitive_taxon_labels=self.case_sensitive_taxon_labels):
            tree_list.append(t, reindex_taxa=False)
        return self.dataset

############################################################################
## CLASS: NewickWriter

class NewickWriter(iosys.DataWriter):
    "Implementation of DataWriter for NEWICK files and strings."

    def __init__(self, **kwargs):
        """
        __init__ recognizes the following keywords (in addition to those of
        `DataWriter.__init__`):

            `dataset`
                Data to be written.
            `suppress_leaf_taxon_labels`
                If True, then taxon labels will not be printed for leaves.
                Default is False.
            `suppress_leaf_node_labels`
                If False, then node labels (if available) will be printed
                for leaves. Defaults to True. Note that DendroPy distinguishes
                between *taxon* labels and *node* labels. In a typical NEWICK
                string, taxon labels are printed for leaf nodes, while leaf
                node labels are ignored (hence the default 'True' setting to
                suppress leaf node labels).
            `suppress_internal_taxon_labels`
                If True, then taxon labels will not be printed for internal
                nodes.  Default is False.
                NOTE: this replaces the `internal_labels` argument which has
                been deprecated.
            `suppress_internal_node_labels`
                If True, internal node labels will not be written. Default is
                False.
                NOTE: this replaces the `internal_labels` argument which has
                been deprecated.
            `suppress_rooting`
                If True, will not write rooting statement. Default is False.
                NOTE: this replaces the `write_rooting` argument which has been
                deprecated.
            `suppress_edge_lengths`
                If True, will not write edge lengths. Default is False.
                NOTE: this replaces the `edge_lengths` argument which has been
                deprecated.
            `unquoted_underscores`
                If True, labels with underscores will not be quoted, which will
                mean that they will be interpreted as spaces if read again
                ("soft" underscores).  If False, then labels with underscores
                will be quoted, resulting in "hard" underscores.  Default is
                False.
                NOTE: this replaces the `quote_underscores` argument which has
                been deprecated.
            `preserve_spaces`
                If True, spaces not mapped to underscores in labels (which
                means any labels containing spaces will have to be
                quoted). Default is False.
                False.
            `store_tree_weights`
                If True, tree weights are written. Default is False.
            `supppress_annotations`
                If False, will write annotations as comments. Default is True.
            `annotations_as_nhx`
                If True, and if `suppress_annotations` is False, will write
                annotation as NHX statements. Default is False.
            `suppress_item_comments`
                If False, will write any additional comments. Default is True.
            `node_label_element_separator`
                If both `suppress_leaf_taxon_labels` and
                `suppress_leaf_node_labels` are False, then this will be the
                string used to join them. Defaults to ' '.
            `node_label_compose_func`
                If not None, should be a function that takes a Node object as
                an argument and returns the string to be used to represent the
                node in the tree statement. The return value
                from this function is used unconditionally to print a node
                representation in a tree statement, by-passing the default
                labelling function (and thus ignoring
                `suppress_leaf_taxon_labels`, `suppress_leaf_node_labels=True`,
                `suppress_internal_taxon_labels`, `suppress_internal_node_labels`,
                etc.). Defaults to None.
            `edge_label_compose_func`
                If not None, should be a function that takes an Edge object as
                an argument, and returns the string to be used to represent the
                edge length in the tree statement.

        Typically, these keywords would be passed to the `write_to_path()`,
        `write_to_stream` or `as_string` arguments, when 'newick' is used as
        the schema::

            d.write_to_path('data.tre', 'newick',
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
        iosys.DataWriter.__init__(self, **kwargs)

        self.suppress_leaf_taxon_labels = kwargs.get("suppress_leaf_taxon_labels", False)
        self.suppress_leaf_node_labels = kwargs.get("suppress_leaf_node_labels", True)
        self.suppress_internal_taxon_labels = kwargs.get("suppress_internal_taxon_labels", False)
        self.suppress_internal_taxon_labels = not kwargs.get("internal_labels", not self.suppress_internal_taxon_labels) # legacy
        self.suppress_internal_node_labels = kwargs.get("suppress_internal_node_labels", False)
        self.suppress_internal_node_labels = not kwargs.get("internal_labels", not self.suppress_internal_node_labels) # legacy

        self.suppress_rooting = kwargs.get("suppress_rooting", False)
        self.suppress_rooting = not kwargs.get("write_rooting", not self.suppress_rooting) # legacy

        self.suppress_edge_lengths = kwargs.get("suppress_edge_lengths", False)
        self.suppress_edge_lengths = not kwargs.get("edge_lengths", not self.suppress_edge_lengths) # legacy

        self.unquoted_underscores = kwargs.get('unquoted_underscores', False)
        self.unquoted_underscores = not kwargs.get('quote_underscores', not self.unquoted_underscores) # legacy

        self.preserve_spaces = kwargs.get("preserve_spaces", False)
        self.store_tree_weights = kwargs.get("store_tree_weights", False)

        self.suppress_annotations = kwargs.get("suppress_annotations", True)
        self.suppress_annotations = not kwargs.get("annotations_as_comments", not self.suppress_annotations) # legacy

        self.annotations_as_nhx = kwargs.get("annotations_as_nhx", False)

        self.suppress_item_comments = kwargs.get("suppress_item_comments", True)
        self.suppress_item_comments = not kwargs.get("write_item_comments", not self.suppress_item_comments)

        self.node_label_element_separator = kwargs.get("node_label_element_separator", ' ')
        self.node_label_compose_func = kwargs.get("node_label_compose_func", None)
        self.edge_label_compose_func = kwargs.get("edge_label_compose_func", None)
        if self.edge_label_compose_func is None:
            self.edge_label_compose_func = self._format_edge_length

    def write(self, stream):
        """
        Writes attached `DataSource` or `TaxonDomain` to a destination given
        by the file-like object `stream`.
        """

        if self.exclude_trees:
            return

        assert self.dataset is not None, \
            "NewickWriter instance is not attached to a DataSet: no source of data"

        for tree_list in self.dataset.tree_lists:
            if self.attached_taxon_set is None or self.attached_taxon_set is tree_list.taxon_set:
                self.write_tree_list(tree_list, stream)

    def write_tree_list(self, tree_list, stream):
        """
        Writes a `TreeList` in NEWICK schema to `stream`.
        """
        if self.exclude_trees:
            return
        for tree in tree_list:
            self.write_tree(tree, stream)

    def compose_comment_string(self, item):
        if not self.suppress_item_comments and item.comments:
            item_comments = []
            if isinstance(item.comments, str):
                item.comments = [item.comments]
            for comment in item.comments:
                item_comments.append("[%s]" % comment)
            item_comment_str = "".join(item_comments)
        else:
            item_comment_str = ""
        return item_comment_str

    def write_tree(self, tree, stream):
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
            weight = "[&W %s] " % tree.weight
        else:
            weight = ""
        if not self.suppress_annotations or self.annotations_as_nhx:
            annotation_comments = nexustokenizer.format_annotation_as_comments(tree, nhx=self.annotations_as_nhx)
        else:
            annotation_comments = ""
        tree_comments = self.compose_comment_string(tree)
        stream.write("%s%s%s%s%s;\n" % (rooting,
                weight,
                annotation_comments,
                tree_comments,
                self.compose_node(tree.seed_node)))

    def compose_tree(self, tree):
        """Convienience method."""
        stream = StringIO()
        self.write_tree(tree, stream)
        return stream.getvalue()

    def choose_display_tag(self, node):
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
            tag = textutils.escape_nexus_token(tag,
                    preserve_spaces=self.preserve_spaces,
                    quote_underscores=not self.unquoted_underscores)
            return tag
        else:
            return ""

    def _format_edge_length(self, edge):
        return "%s" % edge.length

    def compose_node(self, node):
        """
        Given a DendroPy Node, this returns the Node as a NEWICK
        statement according to the class-defined formatting rules.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            subnodes = [self.compose_node(child) for child in child_nodes]
            statement = '(' + ','.join(subnodes) + ')'
            if not (self.suppress_internal_taxon_labels and self.suppress_internal_node_labels):
                statement = statement + self.choose_display_tag(node)
            if node.edge and node.edge.length != None and not self.suppress_edge_lengths:
                statement =  "%s:%s" % (statement, self.edge_label_compose_func(node.edge))
        else:
            statement = self.choose_display_tag(node)
            if node.edge and node.edge.length != None and not self.suppress_edge_lengths:
                statement =  "%s:%s" % (statement, self.edge_label_compose_func(node.edge))
        if not self.suppress_annotations or self.annotations_as_nhx:
            node_annotation_comments = nexustokenizer.format_annotation_as_comments(node, nhx=self.annotations_as_nhx)
            edge_annotation_comments = nexustokenizer.format_annotation_as_comments(node.edge, nhx=self.annotations_as_nhx)
            statement = statement + node_annotation_comments + edge_annotation_comments
        #if self.nhx_key_to_func:
        #    nhx_to_print = []
        #    for k, v in self.nhx_key_to_func.items():
        #        r = v(node.edge)
        #        if r is not None:
        #            nhx_to_print.append("%s=%s" % (k, str(r)))
        #    if nhx_to_print:
        #        statement = statement + ('[&&NHX:%s]' % ':'.join(nhx_to_print))
        node_comment_str = self.compose_comment_string(node)
        statement += node_comment_str
        return statement
