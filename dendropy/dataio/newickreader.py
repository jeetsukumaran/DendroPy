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

import re
import warnings
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.datamodel import base
from dendropy.utility import error
from dendropy.dataio import nexusprocessing
from dendropy.dataio import ioservice

##############################################################################
## NewickReader

class NewickReader(ioservice.DataReader):
    """
    Parser for NEWICK-formatted data.
    """

    FIGTREE_COMMENT_FIELD_PATTERN = re.compile(r'(.+?)=({.+?,.+?}|.+?)(,|$)')
    NHX_COMMENT_FIELD_PATTERN = re.compile(r'(.+?)=({.+?,.+?}|.+?)(:|$)')

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

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------
        rooting : string, {['default-unrooted'], 'default-rooted', 'force-unrooted', 'force-rooted'}
            Specifies how trees in the data source should be intepreted with
            respect to their rooting:

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
            objects. Defaults to `True`: internal node labels will *not* be
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
            msg = StringIO()
            error.dump_stack(msg)
            warnings.warn("\n{}\nUse of keyword argument '{}={}' is deprecated; use 'rooting=\"{}\"' instead".format(msg.getvalue(), kw, kwargs[kw], corrected),
                    FutureWarning, stacklevel=4)
            kwargs.pop(kw)
            kwargs["rooting"] = corrected
        self.rooting = kwargs.pop("rooting", "default-unrooted")
        self.edge_len_type = kwargs.get("edge_len_type", float)
        self.extract_comment_metadata = kwargs.get('extract_comment_metadata', False)
        self.store_tree_weights = kwargs.get("store_tree_weights", False)
        self.encode_splits = kwargs.get("encode_splits", False)
        self.finish_node_func = kwargs.get("finish_node_func", None)
        self.case_sensitive_taxon_labels = kwargs.get('case_sensitive_taxon_labels', False)
        self.preserve_underscores = kwargs.get('preserve_underscores', False)
        self.suppress_internal_node_taxa = kwargs.get("suppress_internal_node_taxa", True)
        self.allow_duplicate_taxon_labels = kwargs.get("allow_duplicate_taxon_labels", False)
        self.hyphens_as_tokens =  kwargs.get('hyphens_as_tokens', False)

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
            tree = self._parse_tree_statement(
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

    def _get_rooting(self):
        """
        Get rooting interpretation configuration.
        """
        return self._rooting
    def _set_rooting(self, val):
        """
        Set rooting interpretation configuration.
        """
        if val not in ["force-unrooted", "force-rooted", "default-unrooted", "default-rooted"]:
            raise ValueError("Unrecognized rooting directive: '{}'".format(val))
        self._rooting = val
    rooting = property(_get_rooting, _set_rooting)

    def _parse_tree_statement(self,
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
        tree_comments = nexus_tokenizer.pull_captured_comments()
        if current_token is None:
            current_token = nexus_tokenizer.next_token()
            tree_comments = nexus_tokenizer.pull_captured_comments()
        while current_token == ";" and not nexus_tokenizer.is_eof():
            current_token = nexus_tokenizer.require_next_token()
            tree_comments = nexus_tokenizer.pull_captured_comments()
        if nexus_tokenizer.is_eof():
            return None
        if current_token != "(":
            raise NewickReader.NewickReaderInvalidTokenError(
                    message="Expecting '{}' but found '{}'".format("(", current_token),
                    line_num=nexus_tokenizer.token_line_num,
                    col_num=nexus_tokenizer.token_column_num,
                    stream=nexus_tokenizer.src)
        tree = tree_factory()
        self._process_tree_comments(tree, tree_comments)
        self._parse_tree_node_description(
                nexus_tokenizer=nexus_tokenizer,
                tree=tree,
                current_node=tree.seed_node,
                taxon_symbol_map_func=taxon_symbol_map_func)
        current_token = nexus_tokenizer.current_token
        while current_token == ";" and not nexus_tokenizer.is_eof():
            current_token = nexus_tokenizer.next_token()
        return tree

    def _process_tree_comments(self, tree, tree_comments):
        for comment in tree_comments:
            if comment in ["&u", "&U", "&r", "&R"]:
                tree.is_rooted = self._parse_tree_rooting_state(comment)
            elif comment.startswith("&W") or comment.startswith("&w"):
                if self.store_tree_weights:
                    try:
                        weight_expression = stream_tokenizer.tree_weight_comment.split(' ')[1]
                        tree.weight = eval("/".join(["float(%s)" % cv for cv in weight_expression.split('/')]))
                    except IndexError:
                        pass
                    except ValueError:
                        pass
                else:
                    # if tree weight comment is not processed,
                    # just store it
                    tree.comments.append(comment)
            elif self.extract_comment_metadata and comment.startswith("&"):
                annotations = self._parse_comment_metadata(comment)
                if annotations:
                    tree.annotations.update(annotations)
                else:
                    tree.comments.append(comment)
            else:
                tree.comments.append(comment)

    def _parse_tree_rooting_state(self, rooting_comment=None):
        """
        Returns rooting state for tree with given rooting comment token, taking
        into account `rooting` configuration.
        """
        if self._rooting == "force-unrooted":
            return False
        elif self._rooting == "force-rooted":
            return True
        elif rooting_comment == "&R" or rooting_comment == "&r":
            return True
        elif rooting_comment == "&U" or rooting_comment == "&U":
            return False
        elif self._rooting == "default-rooted":
            return True
        elif self._rooting == "default-unrooted" or self._rooting is None:
            return False
        else:
            raise TypeError("Unrecognized rooting directive: '{}'".format(self._rooting))

    def _parse_tree_node_description(
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
        current_node_comments = nexus_tokenizer.pull_captured_comments()
        if nexus_tokenizer.current_token == "(":
            nexus_tokenizer.require_next_token()
            node_created = False
            while True:
                if nexus_tokenizer.current_token == ",":
                    if not node_created: #184
                        # no node has been created yet: ',' designates a
                        # preceding blank node
                        new_node = tree.node_factory()
                        self._process_node_comments(node=new_node,
                                nexus_tokenizer=nexus_tokenizer)
                        current_node.add_child(new_node)
                        # do not flag node as created to allow for an extra node to be created in the event of (..,)
                    nexus_tokenizer.require_next_token()
                    while nexus_tokenizer.current_token == ",": #192
                        # another blank node
                        new_node = tree.node_factory()
                        self._process_node_comments(node=new_node,
                                nexus_tokenizer=nexus_tokenizer)
                        current_node.add_child(new_node)
                        nexus_tokenizer.require_next_token()
                        node_created = true;
                    if not node_created and nexus_tokenizer.current_token == ")": #200
                        # end of node
                        new_node = tree.node_factory();
                        self._process_node_comments(node=new_node,
                                nexus_tokenizer=nexus_tokenizer)
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
                    self._process_node_comments(node=new_node,
                            nexus_tokenizer=nexus_tokenizer)
                    self._parse_tree_node_description(
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
                    ### TODO!!! handle other types of valuesNewickReader.
                    raise
                current_node.edge.length = edge_length
                nexus_tokenizer.require_next_token()
            elif nexus_tokenizer.current_token == ")": #253
                # closing of parent token
                self._process_node_comments(node=current_node,
                        nexus_tokenizer=nexus_tokenizer,
                        additional_comments=current_node_comments)
                return current_node
            elif nexus_tokenizer.current_token == ";": #256
                # end of tree statement
                nexus_tokenizer.next_token()
                break
            elif nexus_tokenizer.current_token == ",": #260
                # end of this node
                self._process_node_comments(node=current_node,
                        nexus_tokenizer=nexus_tokenizer,
                        additional_comments=current_node_comments)
                return current_node
            elif nexus_tokenizer.current_token == "(": #263
                # start of another node or tree without finishing this
                # node
                raise NewickReader.NewickReaderMalformedStatementError(
                        message="Malformed tree statement",
                        line_num=nexus_tokenizer.token_line_num,
                        col_num=nexus_tokenizer.token_column_num,
                        stream=nexus_tokenizer.src)
            else: #267
                # label
                if label_parsed: #269
                    raise NewickReader.NewickReaderMalformedStatementError(
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
        self._process_node_comments(node=current_node,
                nexus_tokenizer=nexus_tokenizer,
                additional_comments=current_node_comments)
        return current_node

    def _parse_comment_metadata(comment,
            annotations=None,
            field_name_map=None,
            field_value_types=None,
            strip_leading_trailing_spaces=True):
        """
        Returns set of :class:`Annotation` objects corresponding to metadata
        given in comments.

        Parameters
        ----------
        `comment` : string
            A comment token.
        `annotations` : :class:`AnnotationSet` or `set`
            Set of :class:`Annotation` objects to which to add this annotation.
        `field_name_map` : dict
            A dictionary mapping field names (as given in the comment string)
            to strings that should be used to represent the field in the
            metadata dictionary; if not given, no mapping is done (i.e., the
            comment string field name is used directly).
        `field_value_types` : dict
            A dictionary mapping field names (as given in the comment
            string) to the value type (e.g. {"node-age" : float}.
        `strip_leading_trailing_spaces` : boolean
            Remove whitespace from comments.

        Returns
        -------
        metadata : :py:class::`set` [:class:`Annotation`]
            Set of :class:`Annotation` objects corresponding to metadata
            parsed.
        """
        if annotations is None:
            annotations = set()
        if field_name_map is None:
            field_name_map = {}
        if field_value_types is None:
            field_value_types = {}
        if comment.startswith("&&NHX:"):
            pattern = NewickReader.NHX_COMMENT_FIELD_PATTERN
            comment = comment[6:]
        elif comment.startswith("&&"):
            pattern = NewickReader.NHX_COMMENT_FIELD_PATTERN
            comment = comment[2:]
        elif comment.startswith("&"):
            pattern = NewickReader.FIGTREE_COMMENT_FIELD_PATTERN
            comment = comment[1:]
        else:
            # unrecognized metadata pattern
            return annotations
        for match_group in pattern.findall(comment):
            key, val = match_group[:2]
            if strip_leading_trailing_spaces:
                key = key.strip()
                val = val.strip()
            if key in field_value_types:
                value_type = field_value_types[key]
            else:
                value_type = None
            if val.startswith('{'):
                if value_type is not None:
                    val = [value_type(v) for v in val[1:-1].split(',')]
                else:
                    val = val[1:-1].split(',')
            else:
                if value_type is not None:
                    val = value_type(val)
            if key in field_name_map:
                key = field_name_map[key]
            annote = base.Annotation(
                    name=key,
                    value=val,
                    # datatype_hint=datatype_hint,
                    # name_prefix=name_prefix,
                    # namespace=namespace,
                    # name_is_prefixed=name_is_prefixed,
                    # is_attribute=False,
                    # annotate_as_reference=annotate_as_reference,
                    # is_hidden=is_hidden,
                    )
            annotations.add(annote)
        return annotations

    def _process_node_comments(self,
            node,
            nexus_tokenizer,
            additional_comments=None):
        node_comments = nexus_tokenizer.pull_captured_comments()
        if additional_comments is not None:
            node_comments.extend(additional_comments)
        for comment in node_comments:
            if self.extract_comment_metadata and comment.startswith("&"):
                annotations = self._parse_comment_metadata(comment)
                if annotations:
                    node.annotations.update(annotations)
                else:
                    node.comments.append(comment)
            else:
                node.comments.append(comment)
