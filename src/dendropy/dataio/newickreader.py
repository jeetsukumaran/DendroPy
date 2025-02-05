#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
#
##############################################################################

"""
Parsing of NEWICK-format tree from a stream.
"""

from io import StringIO
import itertools as it
from dendropy.utility import error
from dendropy.utility import deprecate
from dendropy.dataio import tokenizer
from dendropy.dataio import nexusprocessing
from dendropy.dataio import ioservice

##############################################################################
## NewickReader

class NewickReader(ioservice.DataReader):
    """
    Parser for NEWICK-formatted data.
    """

    _default_rooting_directive = None
    _default_tree_weight = 1.0

    class NewickReaderError(error.DataParseError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            error.DataParseError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NewickReaderInvalidTokenError(NewickReaderError):
        def __init__(self,
                message,
                line_num=None,
                col_num=None,
                stream=None):
            NewickReader.NewickReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NewickReaderMalformedStatementError(NewickReaderError):
        def __init__(self,
                message,
                line_num=None,
                col_num=None,
                stream=None):
            NewickReader.NewickReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NewickReaderIncompleteTreeStatementError(NewickReaderMalformedStatementError):
        def __init__(self,
                message,
                line_num=None,
                col_num=None,
                stream=None):
            NewickReader.NewickReaderMalformedStatementError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NewickReaderInvalidValueError(NewickReaderError):
        def __init__(self,
                message,
                line_num=None,
                col_num=None,
                stream=None):
            NewickReader.NewickReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NewickReaderDuplicateTaxonError(NewickReaderError):
        def __init__(self,
                message,
                line_num=None,
                col_num=None,
                stream=None):
            detailed = ("Multiple occurrences of the same taxa on trees are not"
            " supported: trees with duplicate node labels can only be"
            " processed if the labels are not parsed as operational taxonomic"
            " unit concepts but instead as simply node labels by specifying"
            " 'suppress_internal_node_taxa=True, suppress_leaf_node_taxa=True'."
            " Duplicate taxon labels: {}").format(message)
            NewickReader.NewickReaderError.__init__(self,
                    message=detailed,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    def __init__(self, **kwargs):
        """Keyword Arguments
        -----------------
        rooting : string, {['default-unrooted'], 'default-rooted', 'force-unrooted', 'force-rooted'}
            Specifies how trees in the data source should be intepreted with
            respect to their rooting:

                'default-unrooted' [default]:
                    All trees are interpreted as unrooted unless a '[&R]'
                    comment token explicitly specifies them as rooted.
                'default-rooted'
                    All trees are interpreted as rooted unless a '[&U]'
                    comment token explicitly specifies them as unrooted.
                'force-unrooted'
                    All trees are unconditionally interpreted as unrooted.
                'force-rooted'
                    All trees are unconditionally interpreted as rooted.

        edge_length_type : type, default: ``float``
            Specifies the type of the edge lengths (``int`` or ``float``). Tokens
            interpreted as branch lengths will be cast to this type.
            Defaults to ``float``.
        suppress_edge_lengths : boolean, default: |False|
            If |True|, edge length values will not be processed. If |False|,
            edge length values will be processed.
        extract_comment_metadata : boolean, default: |True|
            If |True| (default), any comments that begin with '&' or '&&' will
            be parsed and stored as part of the annotation set of the
            corresponding object (accessible through the ``annotations``
            attribute of the object). This requires that the comment
            contents conform to a particular format (NHX or BEAST: 'field =
            value'). If |False|, then the comments will not be parsed,
            but will be instead stored directly as elements of the ``comments``
            list attribute of the associated object.
        store_tree_weights : boolean, default: |False|
            If |True|, process the tree weight (e.g. "[&W 1/2]") comment
            associated with each tree, if any. Defaults to |False|.
        finish_node_fn : function object, default: |None|
            If specified, this function will be applied to each node after
            it has been constructed.
        case_sensitive_taxon_labels : boolean, default: |False|
            If |True|, then taxon labels are case sensitive (e.g., "P.regius"
            and "P.REGIUS" wil be treated as different operation taxonomic
            unit concepts). Otherwise, taxon label intepretation will be made
            without regard for case.
        preserve_underscores : boolean, default: |False|
            If |True|, unquoted underscores in labels will *not* converted to
            spaces. Defaults to |False|: all underscores not protected by
            quotes will be converted to spaces.
        suppress_internal_node_taxa : boolean, default: |True|
            If |False|, internal node labels will be instantantiated into
            |Taxon| objects. If |True|, internal node labels
            will *not* be instantantiated as strings.
        suppress_leaf_node_taxa : boolean, default: |False|
            If |False|, leaf (external) node labels will be instantantiated
            into |Taxon| objects. If |True|, leaf (external) node
            labels will *not* be instantantiated as strings.
        is_parse_jplace_tokens : boolean: |False|
            If |True|, then accept edge numbering according to the jplace
            format, as described in Matsen et. al. PLoS One, 2012
            http://dx.doi.org/10.1371/journal.pone.0031009. An instance variable
            edge_index is added to the returned tree, and an edge_number is
            added to each edge. If False [default], encountering edge labels
            raises a NewickReaderMalformedStatementError.
        is_assign_internal_labels_to_edges : boolean, default: |None|
            If |True|, internal node labels will be assigned as edge labels.
        terminating_semicolon_required : boolean, default: |True|
            If |True| [default], then a tree statement that does not end in a
            semi-colon is an error. If |False|, then no error will be raised.
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.

        """

        # base
        ioservice.DataReader.__init__(self)

        self._rooting = None
        ## (TEMPORARY and UGLY!!!!) Special handling for legacy signature
        if "as_unrooted" in kwargs or "as_rooted" in kwargs or "default_as_rooted" in kwargs or "default_as_unrooted" in kwargs:
            import collections
            legacy_kw = ("as_unrooted", "as_rooted", "default_as_rooted", "default_as_unrooted")
            legacy_kw_str = ", ".join("'{}'".format(k) for k in legacy_kw)
            if "rooting" in kwargs:
                raise ValueError("Cannot specify 'rooting' keyword argument in conjunction with any of the (legacy) keyword arguments ({}). Use 'rooting' alone.".format(legacy_kw_str))
            specs = collections.Counter(k for k in kwargs.keys() if k in legacy_kw)
            if sum(specs.values()) > 1:
                raise ValueError("Cannot specify more than one of {{ {} }} at the same time".format(legacy_kw_str))
            kw = list(specs.keys())[0]
            if kw == "as_unrooted":
                if kwargs[kw]:
                    corrected = "force-unrooted"
                else:
                    corrected = "force-rooted"
            elif kw == "as_rooted":
                if kwargs[kw]:
                    corrected = "force-rooted"
                else:
                    corrected = "force-unrooted"
            elif kw == "default_as_unrooted":
                if kwargs[kw]:
                    corrected = "default-unrooted"
                else:
                    corrected = "default-rooted"
            elif kw == "default_as_rooted":
                if kwargs[kw]:
                    corrected = "default-rooted"
                else:
                    corrected = "default-unrooted"
            deprecate.dendropy_deprecation_warning(
                    preamble="Deprecated since DendroPy 4:",
                    old_construct="{}={}".format(kw, kwargs[kw]),
                    new_construct="rooting='{}'".format(corrected),
                    stacklevel=7)
            kwargs.pop(kw)
            kwargs["rooting"] = corrected
        if "allow_duplicate_taxon_labels" in kwargs:
            raise ValueError(
                "'allow_duplicate_taxon_labels' is no longer"
                " supported: trees with duplicate node labels can only be"
                " processed if the labels are not parsed as operational taxonomic"
                " unit concepts but instead as simply node labels by specifying"
                " 'suppress_internal_node_taxa=True, suppress_leaf_node_taxa=True'."
            )
        # self.rooting = kwargs.pop("rooting", "default-unrooted")
        self.rooting = kwargs.pop("rooting", self.__class__._default_rooting_directive)
        self.edge_length_type = kwargs.pop("edge_length_type", float)
        self.suppress_edge_lengths = kwargs.pop("suppress_edge_lengths", False)
        self.extract_comment_metadata = kwargs.pop('extract_comment_metadata', True)
        self.store_tree_weights = kwargs.pop("store_tree_weights", False)
        self.default_tree_weight = kwargs.pop("default_tree_weight", self.__class__._default_tree_weight)
        self.finish_node_fn = kwargs.pop("finish_node_fn", None)
        self.case_sensitive_taxon_labels = kwargs.pop('case_sensitive_taxon_labels', False)
        self.preserve_unquoted_underscores = kwargs.pop('preserve_underscores', False)
        self.suppress_internal_node_taxa = kwargs.pop("suppress_internal_node_taxa", True)
        self.suppress_leaf_node_taxa = kwargs.pop("suppress_external_node_taxa", False) # legacy (will be deprecated)
        self.suppress_leaf_node_taxa = kwargs.pop("suppress_leaf_node_taxa", self.suppress_leaf_node_taxa)
        self.is_parse_jplace_tokens = kwargs.pop("is_parse_jplace_tokens", False)
        self.is_assign_internal_labels_to_edges = kwargs.pop("is_assign_internal_labels_to_edges", None)
        if self.is_assign_internal_labels_to_edges and not self.suppress_internal_node_taxa:
            raise ValueError("Conflicting options: cannot simultaneously assign internal labels to edges and to internal taxa")
        self.terminating_semicolon_required = kwargs.pop("terminating_semicolon_required", True)
        self.check_for_unused_keyword_arguments(kwargs)

        # per-tree book-keeping
        self._tree_statement_complete = None
        self._parenthesis_nesting_level = None
        self._seen_taxa = None

    def tree_iter(self,
            stream,
            taxon_symbol_mapper,
            tree_factory):
        """
        Iterator that yields trees in NEWICK-formatted source.

        Parameters
        ----------
        stream : file or file-like object
            A file or file-like object opened for reading.
        taxon_namespace : |TaxonNamespace|
            Operational taxonomic unit namespace to use for taxon management.
        tree_factory : function object
            A function that returns a new |Tree| object when called
            without arguments.

        Returns
        -------
        iter : :py`collections.Iterator` [|Tree|]
            An iterator yielding |Tree| objects constructed based on
            data in ``stream``.
        """
        nexus_tokenizer = nexusprocessing.NexusTokenizer(stream,
                preserve_unquoted_underscores=self.preserve_unquoted_underscores)
        while True:
            tree = self._parse_tree_statement(
                    nexus_tokenizer=nexus_tokenizer,
                    tree_factory=tree_factory,
                    taxon_symbol_map_fn=taxon_symbol_mapper.require_taxon_for_symbol)
            yield tree
            if tree is None:
                # raise StopIteration
                return

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            state_alphabet_factory=None,
            global_annotations_target=None):
        taxon_namespace = taxon_namespace_factory(label=None)
        tree_list = tree_list_factory(label=None, taxon_namespace=taxon_namespace)
        taxon_symbol_mapper = nexusprocessing.NexusTaxonSymbolMapper(
                taxon_namespace=taxon_namespace,
                enable_lookup_by_taxon_number=False,
                case_sensitive=self.case_sensitive_taxon_labels)
        tree_factory = tree_list.new_tree
        for tree in self.tree_iter(stream=stream,
                taxon_symbol_mapper=taxon_symbol_mapper,
                tree_factory=tree_factory):
            pass
        product = self.Product(
                taxon_namespaces=None,
                tree_lists=[tree_list],
                char_matrices=None)
        return product

    def _get_rooting(self):
        """
        Get rooting interpretation configuration.
        """
        return self._rooting
    def _set_rooting(self, val):
        """
        Set rooting interpretation configuration.
        """
        if val not in ("force-unrooted", "force-rooted", "default-unrooted", "default-rooted", None,):
            raise ValueError("Unrecognized rooting directive: '{}'".format(val))
        self._rooting = val
    rooting = property(_get_rooting, _set_rooting)

    def _parse_tree_statement(self,
            nexus_tokenizer,
            tree_factory,
            taxon_symbol_map_fn):
        """
        Parses a single tree statement from a token stream and constructs a
        corresponding Tree object. Expects that the first non-comment and
        non-semi-colon token to be found, including the current token, to be
        the parenthesis that opens the tree statement. When complete, the
        current token will be the token immediately following the semi-colon,
        if any.
        """
        current_token = nexus_tokenizer.current_token
        tree_comments = nexus_tokenizer.pull_captured_comments()
        while (current_token == ";" or current_token is None) and not nexus_tokenizer.is_eof():
            current_token = nexus_tokenizer.require_next_token()
            tree_comments = nexus_tokenizer.pull_captured_comments()
        if nexus_tokenizer.is_eof():
            return None
        if current_token != "(":
            # allow for possibility of single node tree, e.g.: T0:10;
            self._parenthesis_nesting_level = 0
            # raise NewickReader.NewickReaderInvalidTokenError(
            #         message="Expecting '{}' but found '{}'".format("(", current_token),
            #         line_num=nexus_tokenizer.token_line_num,
            #         col_num=nexus_tokenizer.token_column_num,
            #         stream=nexus_tokenizer.src)
        else:
            self._parenthesis_nesting_level = 1
        tree = tree_factory()
        self._process_tree_comments(tree, tree_comments, nexus_tokenizer)
        self._tree_statement_complete = False
        self._seen_taxa = set()
        self._parse_tree_node_description(
                nexus_tokenizer=nexus_tokenizer,
                tree=tree,
                current_node=tree.seed_node,
                taxon_symbol_map_fn=taxon_symbol_map_fn,
                is_internal_node=None)
        current_token = nexus_tokenizer.current_token
        if not self._tree_statement_complete:
            raise NewickReader.NewickReaderIncompleteTreeStatementError(
                    message="Incomplete or improperly-terminated tree statement (last character read was '{}' instead of a semi-colon ';')".format(nexus_tokenizer.current_token),
                    line_num=nexus_tokenizer.token_line_num,
                    col_num=nexus_tokenizer.token_column_num,
                    stream=nexus_tokenizer.src)
        self._seen_taxa = None
        self._parenthesis_nesting_level = None
        self._tree_statement_complete = None
        while current_token == ";" and not nexus_tokenizer.is_eof():
            nexus_tokenizer.clear_captured_comments()
            current_token = nexus_tokenizer.next_token()
        return tree

    def _process_tree_comments(self, tree, tree_comments, nexus_tokenizer):
        # NOTE: this also unconditionally sets the tree rootedness and
        # weighting if no comment indicating these are found; for this to work
        # in the current implementation, this method must be called once and
        # exactly once per tree.
        if not tree_comments:
            tree.is_rooted = self._parse_tree_rooting_state("")
            if self.store_tree_weights:
                tree.weight = self.default_tree_weight
            return
        rooting_token_found = False
        weighting_token_found = False
        for comment in tree_comments:
            stripped_comment = comment.strip()
            if stripped_comment in ["&u", "&U", "&r", "&R"]:
                tree.is_rooted = self._parse_tree_rooting_state(stripped_comment)
                rooting_token_found = True
            elif (self.store_tree_weights
                    and (stripped_comment.startswith("&W ") or stripped_comment.startswith("&w "))
                    ):
                try:
                    weight_expression = stripped_comment[2:]
                    if not weight_expression:
                        raise ValueError
                    we_parts = weight_expression.split("/")
                    if len(we_parts) > 2:
                        raise ValueError
                        # raise NewickReader.NewickReaderInvalidValueError(
                        #         message="Invalid tree weight expression: '{}'".format(weight_expression),
                        #         line_num=nexus_tokenizer.token_line_num,
                        #         col_num=nexus_tokenizer.token_column_num,
                        #         stream=nexus_tokenizer.src)
                    elif len(we_parts) == 2:
                        x = float(we_parts[0])
                        y = float(we_parts[1])
                        tree.weight = x/y
                    else:
                        tree.weight = float(we_parts[0])
                    weighting_token_found = True
                except ValueError:
                    exc = NewickReader.NewickReaderInvalidValueError(
                            message="Invalid tree weight expression: '{}'".format(stripped_comment),
                            line_num=nexus_tokenizer.token_line_num,
                            col_num=nexus_tokenizer.token_column_num,
                            stream=nexus_tokenizer.src)
                    exc.__context__ = None # Python 3.0, 3.1, 3.2
                    exc.__cause__ = None # Python 3.3, 3.4
                    raise exc
            elif self.extract_comment_metadata and comment.startswith("&"):
                annotations = nexusprocessing.parse_comment_metadata_to_annotations(
                    comment=comment)
                if annotations:
                    tree.annotations.update(annotations)
                else:
                    tree.comments.append(comment)
            else:
                tree.comments.append(comment)
        if not rooting_token_found:
            tree.is_rooted = self._parse_tree_rooting_state("")
        if self.store_tree_weights and not weighting_token_found:
            tree.weight = self.default_tree_weight

    def _parse_tree_rooting_state(self, rooting_comment=None):
        """
        Returns rooting state for tree with given rooting comment token, taking
        into account ``rooting`` configuration.
        """
        if self._rooting == "force-unrooted":
            return False
        elif self._rooting == "force-rooted":
            return True
        elif rooting_comment == "&R" or rooting_comment == "&r":
            return True
        elif rooting_comment == "&U" or rooting_comment == "&u":
            return False
        elif self._rooting == "default-rooted":
            return True
        elif self._rooting == "default-unrooted":
            return False
        elif self._rooting is None:
            return None
        else:
            raise TypeError("Unrecognized rooting directive: '{}'".format(self._rooting))

    def _parse_tree_node_description(
            self,
            nexus_tokenizer,
            tree,
            current_node,
            taxon_symbol_map_fn,
            is_internal_node=None):
        """
        Assuming that the iterator is currently sitting on a parenthesis that
        opens a node with children or the label of a leaf node, this will
        populate the node ``node`` appropriately (label, edge length, comments,
        metadata etc.) and recursively parse and add the node's
        children. When complete, the token will be the token immediately
        following the end of the node or tree statement if this is the root
        node, i.e. the token *following* the closing parenthesis of the node or
        the semi-colon terminating a tree statement.
        """
        current_node_comments = nexus_tokenizer.pull_captured_comments()
        if nexus_tokenizer.current_token == "(":
            # self._parenthesis_nesting_level += 1 # handled by calling code
            nexus_tokenizer.require_next_token()
            node_created = False
            for count in it.count():
                if nexus_tokenizer.current_token == ",":
                    if not node_created: #184
                        # no node has been created yet: ',' designates a
                        # preceding blank node
                        new_node = tree.node_factory()
                        nexusprocessing.process_comments_for_item(item=new_node,
                                item_comments=nexus_tokenizer.pull_captured_comments(),
                                extract_comment_metadata=self.extract_comment_metadata)
                        self._finish_node(new_node)
                        current_node.add_child(new_node)
                        ## node_created = True # do not flag node as created to allow for an extra node to be created in the event of (..,)
                    nexus_tokenizer.require_next_token()
                    while nexus_tokenizer.current_token == ",": #192
                        # another blank node
                        new_node = tree.node_factory()
                        nexusprocessing.process_comments_for_item(item=new_node,
                                item_comments=nexus_tokenizer.pull_captured_comments(),
                                extract_comment_metadata=self.extract_comment_metadata)
                        self._finish_node(new_node)
                        current_node.add_child(new_node)
                        # node_created = True; # do not flag node as created: extra node needed in the event of (..,)
                        nexus_tokenizer.require_next_token()
                    if not node_created and nexus_tokenizer.current_token == ")": #200
                        # end of node
                        new_node = tree.node_factory();
                        nexusprocessing.process_comments_for_item(item=new_node,
                                item_comments=nexus_tokenizer.pull_captured_comments(),
                                extract_comment_metadata=self.extract_comment_metadata)
                        self._finish_node(new_node)
                        current_node.add_child(new_node)
                        node_created = True;
                elif nexus_tokenizer.current_token == ")": #206
                    if count == 0:
                        # handle terminating unnamed unifurcation
                        # see https://github.com/jeetsukumaran/DendroPy/issues/76
                        new_node = tree.node_factory()
                        is_new_internal_node = False
                        self._finish_node(new_node)
                        current_node.add_child(new_node)
                    # end of child nodes
                    self._parenthesis_nesting_level -= 1
                    nexus_tokenizer.require_next_token()
                    break
                else: #210
                    # assume child nodes: a leaf node (if a label) or
                    # internal (if a parenthesis)
                    if nexus_tokenizer.current_token == "(":
                        self._parenthesis_nesting_level += 1
                        is_new_internal_node = True
                    else:
                        is_new_internal_node = False
                    new_node = tree.node_factory();
                    nexusprocessing.process_comments_for_item(item=new_node,
                            item_comments=nexus_tokenizer.pull_captured_comments(),
                            extract_comment_metadata=self.extract_comment_metadata)
                    self._parse_tree_node_description(
                            nexus_tokenizer=nexus_tokenizer,
                            tree=tree,
                            current_node=new_node,
                            taxon_symbol_map_fn=taxon_symbol_map_fn,
                            is_internal_node=is_new_internal_node,
                            )
                    current_node.add_child(new_node);
                    node_created = True;
        label_parsed = False
        self._tree_statement_complete = False
        if is_internal_node is None:
            # Initial call using ``seed_node`` does not set ``is_internal_node`` to
            # |True| or |False|, explicitly, but rather |None|. If this is the
            # case, the rest of the tree has be constructed, and we simply look
            # at whether there are children or not to determine if it is an
            # internal node. This approach allows for a single-tip tree.
            if current_node._child_nodes:
                is_internal_node = True
        if current_node_comments is None:
            current_node_comments = []
        while True:
            cc = nexus_tokenizer.pull_captured_comments()
            if cc is not None:
                current_node_comments.extend(cc)
            if nexus_tokenizer.current_token == ":": #246
                nexus_tokenizer.require_next_token()
                if not self.suppress_edge_lengths:
                    try:
                        edge_length = self.edge_length_type(nexus_tokenizer.current_token)
                    except ValueError:
                        raise NewickReader.NewickReaderMalformedStatementError(
                                message="Invalid edge length: '{}'".format(nexus_tokenizer.current_token),
                                line_num=nexus_tokenizer.token_line_num,
                                col_num=nexus_tokenizer.token_column_num,
                                stream=nexus_tokenizer.src)
                    current_node.edge.length = edge_length
                try:
                    nexus_tokenizer.require_next_token()
                except tokenizer.Tokenizer.UnexpectedEndOfStreamError as e:
                    if self.terminating_semicolon_required:
                        message = e.message + ". (Perhaps the terminating semicolon for the tree statement is missing? If so, add a semicolon to the tree statement or specify 'terminating_semicolon_required=False' to allow for missing semicolons)"
                        raise tokenizer.Tokenizer.UnexpectedEndOfStreamError(
                                message=message,
                                line_num=e.line_num,
                                col_num=e.col_num,
                                stream=e.stream)
                    else:
                        self._tree_statement_complete = True
                        break

            elif nexus_tokenizer.current_token == ")": #253
                # closing of parent token
                # self._parenthesis_nesting_level -= 1 # handled by calling code
                nexusprocessing.process_comments_for_item(item=current_node,
                        item_comments=current_node_comments,
                        extract_comment_metadata=self.extract_comment_metadata)
                self._finish_node(current_node)
                return current_node
            elif nexus_tokenizer.current_token == ";": #256
                # end of tree statement
                self._tree_statement_complete = True
                nexus_tokenizer.next_token()
                break
            elif nexus_tokenizer.current_token == ",": #260
                # end of this node
                nexusprocessing.process_comments_for_item(item=current_node,
                            item_comments=current_node_comments,
                            extract_comment_metadata=self.extract_comment_metadata)
                self._finish_node(current_node)
                return current_node
            elif nexus_tokenizer.current_token == "(": #263
                # start of another node or tree without finishing this
                # node
                self._parenthesis_nesting_level += 1
                raise NewickReader.NewickReaderMalformedStatementError(
                        message="Malformed tree statement",
                        line_num=nexus_tokenizer.token_line_num,
                        col_num=nexus_tokenizer.token_column_num,
                        stream=nexus_tokenizer.src)
            elif self.is_parse_jplace_tokens and nexus_tokenizer.current_token == '{':
                # Edge number from .jplace format
                nexus_tokenizer.require_next_token()
                edge_number = int(nexus_tokenizer.current_token)
                edge = current_node.edge
                edge.edge_number = edge_number
                try:
                    tree.edge_index.insert(edge_number, edge)
                except AttributeError:
                    tree.edge_index = []
                    tree.edge_index.insert(edge_number, edge)
                nexus_tokenizer.require_next_token() # for closing '}'
                nexus_tokenizer.require_next_token()
            else: #267
                # label
                if label_parsed: #269
                    msg = "Expecting ':'"
                    if self.is_parse_jplace_tokens:
                        msg += ", '{'"
                    msg += ", ')', ',' or ';' after reading label but found '{}'".format(nexus_tokenizer.current_token)
                    raise NewickReader.NewickReaderMalformedStatementError(
                            message=msg,
                            line_num=nexus_tokenizer.token_line_num,
                            col_num=nexus_tokenizer.token_column_num,
                            stream=nexus_tokenizer.src)
                else:
                    # Label
                    label = nexus_tokenizer.current_token
                    if ( (is_internal_node and self.suppress_internal_node_taxa)
                            or ((not is_internal_node) and self.suppress_leaf_node_taxa) ):
                        if self.is_assign_internal_labels_to_edges:
                            current_node.edge.label = label
                        else:
                            current_node.label = label
                    else:
                        node_taxon = taxon_symbol_map_fn(label)
                        if node_taxon in self._seen_taxa:
                            raise NewickReader.NewickReaderDuplicateTaxonError(
                                    message=node_taxon.label,
                                    line_num=nexus_tokenizer.token_line_num,
                                    col_num=nexus_tokenizer.token_column_num,
                                    stream=nexus_tokenizer.src)
                        self._seen_taxa.add(node_taxon)
                        current_node.taxon = node_taxon
                    label_parsed = True;
                    # nexus_tokenizer.require_next_token()
                    try:
                        nexus_tokenizer.require_next_token()
                    except tokenizer.Tokenizer.UnexpectedEndOfStreamError:
                        if self.terminating_semicolon_required:
                            raise
                        else:
                            break
        ## if we are here, we have reached the end of the tree
        if self._parenthesis_nesting_level != 0:
            raise NewickReader.NewickReaderMalformedStatementError(
                    message="Unbalanced parentheses at tree statement termination: balance index = {}".format(self._parenthesis_nesting_level),
                    line_num=nexus_tokenizer.token_line_num,
                    col_num=nexus_tokenizer.token_column_num,
                    stream=nexus_tokenizer.src)
        nexusprocessing.process_comments_for_item(item=current_node,
                item_comments=current_node_comments,
                extract_comment_metadata=self.extract_comment_metadata)
        self._finish_node(current_node)
        return current_node

    def _finish_node(self, node):
        if self.finish_node_fn is not None:
            self.finish_node_fn(node)
