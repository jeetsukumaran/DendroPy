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
Parsing of NEWICK-format tree from a stream.
"""

from dendropy.utility import error
from dendropy.dataio import nexusprocessing
from dendropy.dataio import ioservice

##############################################################################
## NewickTreeParser

class NewickTreeParser(object):
    """
    Builds and returns a Tree object from a stream of tokens.
    """

    class NewickTreeParserError(error.DataParseError):

        def __init__(self,
                message,
                line_num=None,
                col_num=None,
                stream=None):
            error.DataParseError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NewickTreeParserInvalidTokenError(NewickTreeParserError):

        def __init__(self,
                message,
                line_num=None,
                col_num=None,
                stream=None):
            NewickTreeParser.NewickTreeParserError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NewickTreeParserMalformedStatementError(NewickTreeParserError):

        def __init__(self,
                message,
                line_num=None,
                col_num=None,
                stream=None):
            NewickTreeParser.NewickTreeParserError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    def __init__(self):
        pass

    def parse_tree_statement(self,
            nexus_tokenizer,
            tree_factory,
            taxon_symbol_map_func):
        """
        Parses a single tree statement from a token stream and constructs a
        corresponding Tree object. Expects that the first non-comment and
        non-semi-colon token to be found, including the current token, to be
        the parenthesis that opens the tree statement. When complete, the
        current token will be the token immediately following the semi-colon,
        if any.
        """
        current_token = nexus_tokenizer.current_token
        if current_token is None:
            current_token = nexus_tokenizer.next_token()
        while current_token == ";" and not nexus_tokenizer.is_eof():
            current_token = nexus_tokenizer.require_next_token()
        if nexus_tokenizer.is_eof():
            return None
        if current_token != "(":
            raise NewickTreeParser.NewickTreeParserInvalidTokenError(
                    message="Expecting '{}' but found '{}'".format("(", current_token),
                    line_num=nexus_tokenizer.token_line_num,
                    col_num=nexus_tokenizer.token_column_num,
                    stream=nexus_tokenizer.src)
        tree = tree_factory()
        self.parse_tree_node_description(
                nexus_tokenizer=nexus_tokenizer,
                tree=tree,
                current_node=tree.seed_node,
                taxon_symbol_map_func=taxon_symbol_map_func)
        current_token = nexus_tokenizer.current_token
        while current_token == ";" and not nexus_tokenizer.is_eof():
            current_token = nexus_tokenizer.next_token()
        return tree

    def parse_tree_node_description(
            self,
            nexus_tokenizer,
            tree,
            current_node,
            taxon_symbol_map_func,
            is_internal_node=False):
        """
        Assuming that the iterator is currently sitting on a parenthesis
        that opens a node with children or the label of a leaf node, this
        will populate the node ``node`` appropriately (label, edge length,
        comments, metadata etc.) and recursively parse and add the node's
        children. When complete, the token will be the token immediately
        following the end of the node or tree statement if this is the root
        node, i.e. the token following the closing parenthesis of the node
        semi-colon terminating a tree statement.
        """
        if nexus_tokenizer.current_token == "(":
            nexus_tokenizer.require_next_token()
            node_created = False
            while True:
                if nexus_tokenizer.current_token == ",":
                    if not node_created: #184
                        # no node has been created yet: ',' designates a
                        # preceding blank node
                        current_node.add_child(tree.node_factory())
                        # do not flag node as created to allow for an extra node to be created in the event of (..,)
                    nexus_tokenizer.require_next_token()
                    while nexus_tokenizer.current_token == ",": #192
                        # another blank node
                        current_node.add_child(tree.node_factory())
                        nexus_tokenizer.require_next_token()
                        node_created = true;
                    if node_created and nexus_tokenizer.current_token == ")": #200
                        # end of node
                        current_node.add_child(new_node)
                        node_created = true;
                elif nexus_tokenizer.current_token == ")": #206
                    # end of child nodes
                    nexus_tokenizer.require_next_token()
                    break
                else: #210
                    # assume child nodes: a leaf node (if a label) or
                    # internal (if a parenthesis)
                    if nexus_tokenizer.current_token == "(":
                        is_new_internal_node = True
                    else:
                        is_new_internal_node = False
                    new_node = tree.node_factory();
                    self.parse_tree_node_description(
                            nexus_tokenizer=nexus_tokenizer,
                            tree=tree,
                            current_node=new_node,
                            taxon_symbol_map_func=taxon_symbol_map_func,
                            is_internal_node=is_new_internal_node,
                            )
                    current_node.add_child(new_node);
                    node_created = True;
        label_parsed = False
        while True:
            if nexus_tokenizer.current_token == ":": #246
                nexus_tokenizer.require_next_token()
                try:
                    edge_length = float(nexus_tokenizer.current_token)
                except ValueError:
                    ### TODO!!! handle other types of valuesNewickTreeParser.
                    raise
                current_node.edge.length = edge_length
                nexus_tokenizer.require_next_token()
            elif nexus_tokenizer.current_token == ")": #253
                # closing of parent token
                return current_node
            elif nexus_tokenizer.current_token == ";": #256
                # end of tree statement
                nexus_tokenizer.next_token()
                break
            elif nexus_tokenizer.current_token == ",": #260
                # end of this node
                return current_node
            elif nexus_tokenizer.current_token == "(": #263
                # start of another node or tree without finishing this
                # node
                raise NewickTreeParser.NewickTreeParserMalformedStatementError(
                        message="Malformed tree statement",
                        line_num=nexus_tokenizer.token_line_num,
                        col_num=nexus_tokenizer.token_column_num,
                        stream=nexus_tokenizer.src)
            else: #267
                # label
                if label_parsed: #269
                    raise NewickTreeParser.NewickTreeParserMalformedStatementError(
                            message="Expecting ':', ')', ',' or ';' after reading label but found '{}'".format(nexus_tokenizer.current_token),
                            line_num=nexus_tokenizer.token_line_num,
                            col_num=nexus_tokenizer.token_column_num,
                            stream=nexus_tokenizer.src)
                else:
                    # Label: if leaf node, then taxon; if internal node
                    # then either way
                    ### TODO!!! taxon/label handling
                    if is_internal_node:
                        current_node.label = nexus_tokenizer.current_token
                    else:
                        current_node.taxon = taxon_symbol_map_func(nexus_tokenizer.current_token)
                    label_parsed = True;
                    nexus_tokenizer.require_next_token()
        return current_node

##############################################################################
## NewickReader

class NewickReader(ioservice.DataReader):
    """
    Parser for NEWICK-formatted data.
    """

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------
        rootedness : string, {['default-unrooted'], 'default-rooted', 'force-unrooted', 'force-rooted'}
            Specifies how trees in the data source should be intepreted with
            respect to their rootedness:

                '``default-unrooted``' [default]:
                    All trees are interpreted as unrooted unless a '``[&R]``'
                    comment token explicitly specifies them as rooted.
                '``default-rooted``'
                    All trees are interpreted as rooted unless a '``[&U]``'
                    comment token explicitly specifies them as unrooted.
                '``force-unrooted``'
                    All trees are unconditionally interpreted as unrooted.
                '``force-rooted``'
                    All trees are unconditionally interpreted as rooted.

        edge_len_type : type, default: `float`
            Specifies the type of the edge lengths (`int` or `float`). Tokens
            interpreted as branch lengths will be cast to this type.
            Defaults to `float`.
        extract_comment_metadata : boolean, default: `False`
            If `True`, any comments that begin with '&' or '&&' will be parsed
            and stored as part of the annotation set of the corresponding
            object (accessible through the `annotations` attribute of the
            object). This requires that the comment contents conform to
            a particular format (NHX or BEAST: 'field = value'). If `False`,
            then the comments will not be parsed, but will be instead stored
            directly as elements of the `comments` list attribute of the
            associated object.
        store_tree_weights : boolean, default: `False`
            If `True`, process the tree weight (e.g. "``[&W 1/2]``") comment
            associated with each tree, if any. Defaults to `False`.
        encode_splits : boolean, default: `False`
            If `True`, split hash bitmasks will be calculated and attached to
            the edges.
        finish_node_func : function object, default: `None`
            If specified, this function will be applied to each node after
            it has been constructed.
        case_sensitive_taxon_labels : boolean, default: `False`
            If `True`, then taxon labels are case sensitive (e.g., "``P.regius``"
            and "``P.REGIUS``" wil be treated as different operation taxonomic
            unit concepts). Otherwise, taxon label intepretation will be made
            without regard for case.
        preserve_underscores : boolean, default: `False`
            If `True`, unquoted underscores in labels will *not* converted to
            spaces. Defaults to `False`: all underscores not protected by
            quotes will be converted to spaces.
        suppress_internal_node_taxa : boolean, default: `True`
            If `False`, internal node labels will be instantantiatd into Taxon
            objects.  Defaults to `True`: internal node labels will *not* be
            treated as taxa.
        allow_duplicate_taxon_labels : boolean, default: `False`
            If `True`, then multiple identical taxon labels will be allowed.
            Defaults to `False`: treat multiple identical taxon labels as an
            error.
        hyphens_as_tokens : boolean, default: `False`
            If `False`, hyphens are not treated as special punctuation
            characters (and thus can be used as part of labels or edge length
            values without requiring that the labels be wrapped in
            quotes). If `True`, hyphens will be treated as special
            punctuation characters.
        """
        self._parser = NewickTreeParser()

    def tree_iter(self,
            stream,
            taxon_namespace,
            tree_factory):
        """
        Iterator that yields trees in NEWICK-formatted source.

        Parameters
        ----------
        stream : file or file-like object
            A file or file-like object opened for reading.
        taxon_namespace : :class:`TaxonNamespace`
            Operational taxonomic unit namespace to use for taxon management.
        tree_factory : function object
            A function that returns a new :class:`Tree` object when called
            without arguments.

        Returns
        -------
        iter : :py:class:`collections.Iterator` [:class:`Tree`]
            An iterator yielding :class:`Tree` objects constructed based on
            data in `stream`.
        """
        nexus_tokenizer = nexusprocessing.NexusTokenizer(stream)
        taxon_symbol_mapper = nexusprocessing.NexusTaxonSymbolMapper(taxon_namespace=taxon_namespace)
        while True:
            tree = self._parser.parse_tree_statement(
                    nexus_tokenizer=nexus_tokenizer,
                    tree_factory=tree_factory,
                    taxon_symbol_map_func=taxon_symbol_mapper.lookup_taxon_symbol)
            yield tree
            if tree is None:
                raise StopIteration

    def read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            global_annotations_target=None):
        taxon_namespace = taxon_namespace_factory(label=None)
        tree_list = tree_list_factory(label=None, taxon_namespace=taxon_namespace)
        tree_factory = tree_list.new_tree
        for tree in self.tree_iter(stream=stream,
                taxon_namespace=taxon_namespace,
                tree_factory=tree_factory):
            pass
        product = self.Product(
                taxon_namespaces=None,
                tree_lists=[tree_list],
                char_matrices=None)
        return product


