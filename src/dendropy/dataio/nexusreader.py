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
##
##############################################################################

"""
Implementation of NEXUS-schema data reader.
"""

from dendropy.utility import error
from dendropy.dataio import ioservice
from dendropy.dataio import nexusprocessing
from dendropy.dataio import newickreader

###############################################################################
## NexusReader

class NexusReader(ioservice.DataReader):
    "Encapsulates loading and parsing of a NEXUS schema file."

    class BlockTerminatedException(Exception):
        pass

    class NexusReaderError(error.DataParseError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            error.DataParseError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NotNexusFileError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class LinkRequiredError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NoCharacterBlocksFoundError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class UndefinedBlockError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class MultipleBlockWithSameTitleError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class InvalidCharacterStateSymbolError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class InvalidContinuousCharacterValueError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class TooManyTaxaError(NexusReaderError):

        def __init__(self,
                taxon_namespace,
                max_taxa,
                label,
                line_num=None,
                col_num=None,
                stream=None):
            message = "Cannot add taxon with label '{}': Declared number of taxa ({}) already defined: {}".format(
                            label,
                            max_taxa,
                            str(["{}".format(t.label) for t in taxon_namespace]))
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class UndefinedTaxonError(NexusReaderError):

        def __init__(self,
                taxon_namespace,
                label,
                line_num=None,
                col_num=None,
                stream=None):
            message = "Taxon '{}' is not in the set of defined taxa: {}".format(
                            label,
                            str(["{}".format(t.label) for t in taxon_namespace]))
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class TooManyCharactersError(NexusReaderError):

        def __init__(self,
                max_characters,
                character,
                line_num=None,
                col_num=None,
                stream=None):
            message = "Cannot add '{}' to sequence: declared sequence length ({}) will be exceeded".format(
                    character, max_characters)
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class IncompleteBlockError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    ###########################################################################
    ## Life-cycle and Setup

    def __init__(self, **kwargs):
        """

        Keyword Arguments
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
        encode_splits : boolean, default: |False|
            If |True|, split hash bitmasks will be calculated and attached to
            the edges.
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
            into |Taxon| objects. If |True|, leaff (external) node
            labels will *not* be instantantiated as strings.
        terminating_semicolon_required : boolean, default: |True|
            If |True| [default], then a tree statement that does not end in a
            semi-colon is an error. If |False|, then no error will be raised.
        unconstrained_taxa_accumulation_mode : bool
            If |True|, then no error is raised even if the number of taxon
            names defined exceeds the number of declared taxa (as specified by
            'NTAX'). Defaults to |False|.
        automatically_substitute_missing_taxa_blocks : bool
            If |True| then, if a taxon namespace is linked to by title but is
            not given in the data file, then, if one and exactly one other
            taxon namespace has been given in the data file, this taxon
            namespace will be used; if there are multiple taxon namespaces,
            then if ``automatically_create_missing_taxa_blocks`` is |True| a
            new taxon namespace will be created, otherwise an error is raised.
            Default is |False|: if a taxon namespace is linked to by title but
            is not given in the data file, then an error is raised.
        automatically_create_missing_taxa_blocks : bool
            If |True| then taxon namespaces linked to by title but not given in
            the data file will be automatically created. If |False| taxon
            namespaces linked to by title but not given in the data file will
            result in error.
        exclude_chars : bool
            If |False|, then character data will not be read. Defaults to
            |True|: character data will be read.
        exclude_trees : bool
            If |False|, then tree data will not be read. Defaults to
            |True|: tree data will be read.
        store_ignored_blocks : bool
            If |True|, then ignored NEXUS blocks will be stored under the annotation
            (NOT attribute!) ``ignored_nexus_blocks''.
            To dereference, for e.g.: ``dataset.annotations["ignored_nexus_blocks"]``.
            Defaults to |False|: non-character and tree blocks will not be read.
        attached_taxon_namespace : |TaxonNamespace|
            Unify all operational taxonomic unit definitions in this namespace.
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.
        """

        # base
        ioservice.DataReader.__init__(self)

        # Following are NEXUS-parsing specific (i.e., not used by NEWICK
        # parsers), and need to be removed so as not to cause problems with our
        # keyword validation scheme
        self.exclude_chars = kwargs.pop("exclude_chars", False)
        self.exclude_trees = kwargs.pop("exclude_trees", False)
        self.store_ignored_blocks = kwargs.pop("store_ignored_blocks", False)
        self._data_type = kwargs.pop("data_type", "standard")
        self.attached_taxon_namespace = kwargs.pop("attached_taxon_namespace", None)

        # Following are undocumented for a GOOD reason! They are experimental and subject to change!
        self.unconstrained_taxa_accumulation_mode = kwargs.pop("unconstrained_taxa_accumulation_mode", False)
        self.automatically_create_missing_taxa_blocks = kwargs.pop("automatically_create_missing_taxa_blocks", False)
        self.automatically_substitute_missing_taxa_blocks = kwargs.pop("automatically_substitute_missing_taxa_blocks", False)

        # The following are used by NewickReader in addition to NexusReader, So
        # they are extracted/set here and then forwarded on ...
        self.preserve_underscores = kwargs.get('preserve_underscores', False)
        self.case_sensitive_taxon_labels = kwargs.get('case_sensitive_taxon_labels', False)
        self.extract_comment_metadata = kwargs.get('extract_comment_metadata', True)

        # As above, but the NEXUS format default is different from the NEWICK
        # default, so this rather convoluted approach
        # self.extract_comment_metadata = kwargs.pop('extract_comment_metadata', True)
        # kwargs["extract_comment_metadata"] = self.extract_comment_metadata

        # Create newick handler
        self.newick_reader = newickreader.NewickReader(**kwargs)

        # Set up parsing meta-variables
        self._interleave = False
        self._symbols = ""
        self._gap_char = '-'
        self._missing_char = '?'
        self._match_char = frozenset('.')
        self._file_specified_ntax = None
        self._file_specified_nchar = None
        self._nexus_tokenizer = None
        self._taxon_namespace_factory = None
        self._tree_list_factory = None
        self._char_matrix_factory = None
        self._global_annotations_target = None
        self._taxon_namespaces = []
        self._char_matrices = []
        self._tree_lists = []
        self._product = None
        self._ignored_blocks = []

    ###########################################################################
    ## Reader Implementation

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            state_alphabet_factory=None,
            global_annotations_target=None):
        """
        Instantiates and returns a DataSet object based on the
        NEXUS-formatted contents given in the file-like object ``stream``.
        """
        self._taxon_namespace_factory = taxon_namespace_factory
        self._tree_list_factory = tree_list_factory
        if self._tree_list_factory is None:
            self.exclude_trees = True
        self._char_matrix_factory = char_matrix_factory
        if self._char_matrix_factory is None:
            self.exclude_chars = True
        self._state_alphabet_factory = state_alphabet_factory
        self._global_annotations_target = global_annotations_target
        self._parse_nexus_stream(stream)
        self._product = self.Product(
                taxon_namespaces=self._taxon_namespaces,
                tree_lists=self._tree_lists,
                char_matrices=self._char_matrices)
        if self._global_annotations_target is not None and self._ignored_blocks:
            a = self._global_annotations_target.annotations.find(name="ignored_nexus_blocks")
            if a is None:
                self._global_annotations_target.annotations.add_new(
                        name="ignored_nexus_blocks",
                        value=self._ignored_blocks,
                        datatype_hint="xsd:list",
                        )
            else:
                a.extend(self._ignored_blocks)
        return self._product

    ###########################################################################
    ## Tokenizer Control

    def create_tokenizer(self, stream, **kwargs):
        self._nexus_tokenizer = nexusprocessing.NexusTokenizer(
                stream, **kwargs)
        return self._nexus_tokenizer

    def set_stream(self, stream):
        return self._nexus_tokenizer.set_stream(stream)

    ###########################################################################
    ## Book-keeping Control

    def _nexus_error(self, message, error_type=None):
        if error_type is None:
            error_type = NexusReader.NexusReaderError
        e = error_type(
                message=message,
                line_num=self._nexus_tokenizer.token_line_num,
                col_num=self._nexus_tokenizer.token_column_num,
                stream=self._nexus_tokenizer.src)
        return e

    def _too_many_taxa_error(self, taxon_namespace, label):
        e = NexusReader.TooManyTaxaError(
                taxon_namespace=taxon_namespace,
                max_taxa=self._file_specified_ntax,
                label=label,
                line_num=self._nexus_tokenizer.token_line_num,
                col_num=self._nexus_tokenizer.token_column_num,
                stream=self._nexus_tokenizer.src)
        return e

    def _undefined_taxon_error(self, taxon_namespace, label):
        e = NexusReader.UndefinedTaxonError(
                taxon_namespace=taxon_namespace,
                label=label,
                line_num=self._nexus_tokenizer.token_line_num,
                col_num=self._nexus_tokenizer.token_column_num,
                stream=self._nexus_tokenizer.src)
        return e

    def _too_many_characters_error(self, character):
        e = NexusReader.TooManyCharactersError(
                max_characters=self._file_specified_nchar,
                character=character,
                line_num=self._nexus_tokenizer.token_line_num,
                col_num=self._nexus_tokenizer.token_column_num,
                stream=self._nexus_tokenizer.src)
        return e

    def _debug_print(self, message=None, out=None):
        import sys
        if out is None:
            out = sys.stdout
        if message is None:
            message = ""
        else:
            message = " --- ({})".format(message)
        out.write("--- Current Position: Line {}, Column {}; Current token [starting at line {} and column {}]: '{}'{}\n".format(
            self._nexus_tokenizer.current_line_num,
            self._nexus_tokenizer.current_column_num,
            self._nexus_tokenizer.token_line_num,
            self._nexus_tokenizer.token_column_num,
            self._nexus_tokenizer.current_token,
            message))

    ###########################################################################
    ## Data Management

    def _new_taxon_namespace(self, title=None):
        if self.attached_taxon_namespace is not None:
            return self.attached_taxon_namespace
        taxon_namespace = self._taxon_namespace_factory(label=title)
        self._taxon_namespaces.append(taxon_namespace)
        return taxon_namespace

    def _get_taxon_namespace(self, title=None):
        if self.attached_taxon_namespace is not None:
            return self.attached_taxon_namespace
        if title is None:
            if len(self._taxon_namespaces) == 0:
                return self._new_taxon_namespace(title=title)
            elif len(self._taxon_namespaces) == 1:
                return self._taxon_namespaces[0]
            else:
                raise self._nexus_error("Multiple taxa blocks defined: require 'LINK' statement", NexusReader.LinkRequiredError)
        else:
            found = []
            for tns in self._taxon_namespaces:
                if tns.label is not None and tns.label.upper() == title.upper():
                    found.append(tns)
            if len(found) == 0:
                if self.automatically_substitute_missing_taxa_blocks:
                    if len(self._taxon_namespaces) == 1:
                        return self._taxon_namespaces[0]
                    elif not self.automatically_create_missing_taxa_blocks:
                        raise self._nexus_error("Taxa block with title '{}' not found, and multiple taxa blocks are defined for this file: unable to automatically substitute".format(title), NexusReader.UndefinedBlockError)
                if self.automatically_create_missing_taxa_blocks:
                    return self._new_taxon_namespace(title=title)
                raise self._nexus_error("Taxa block with title '{}' not found".format(title), NexusReader.UndefinedBlockError)
            elif len(found) > 1:
                raise self._nexus_error("Multiple taxa blocks with title '{}' defined".format(title), NexusReader.MultipleBlockWithSameTitleError)
            return found[0]

    def _get_taxon_symbol_mapper(self, taxon_namespace, enable_lookup_by_taxon_number=True):
        taxon_symbol_mapper = nexusprocessing.NexusTaxonSymbolMapper(
                taxon_namespace=taxon_namespace,
                enable_lookup_by_taxon_number=enable_lookup_by_taxon_number,
                case_sensitive=self.case_sensitive_taxon_labels)
        return taxon_symbol_mapper

    def _new_char_matrix(self, data_type, taxon_namespace, title=None):
        # if data_type is None:
        #     data_type = "standard"
        char_matrix = self._char_matrix_factory(
                data_type,
                taxon_namespace=taxon_namespace,
                label=title)
        self._char_matrices.append(char_matrix)
        return char_matrix

    def _new_state_alphabet(self, *args, **kwargs):
        return self._state_alphabet_factory(*args, **kwargs)

    def _get_char_matrix(self, title=None):
        if title is None:
            if len(self._char_matrices) == 1:
                return self._char_matrices[0]
            elif len(self._char_matrices) == 0:
                raise self._nexus_error("No character matrices defined", NexusReader.NoCharacterBlocksFoundError)
            else:
                raise self._nexus_error("Multiple character matrices defined: require 'LINK' statement", NexusReader.LinkRequiredError)
        else:
            found = []
            for cm in self._char_matrices:
                if cm.label.upper() == title.upper():
                    found.append(cm)
            if len(found) == 0:
                raise self._nexus_error("Character block with title '{}' not found".format(title), NexusReader.UndefinedBlockError)
            elif len(found) > 1:
                raise self._nexus_error("Multiple character blocks with title '{}' defined".format(title), NexusReader.MultipleBlockWithSameTitleError)
            return found[0]

    def _new_tree_list(self, taxon_namespace, title=None):
        tree_list = self._tree_list_factory(
                taxon_namespace=taxon_namespace,
                label=title)
        self._tree_lists.append(tree_list)
        return tree_list

    def _get_tree_list(self, title=None):
        if title is None:
            if len(self._tree_lists) == 1:
                return self._tree_lists[0]
            elif len(self._tree_lists) == 0:
                raise self._nexus_error("No tree blocks defined", NexusReader.NoCharacterBlocksFoundError)
            else:
                raise self._nexus_error("Multiple tree blocks defined: require 'LINK' statement", NexusReader.LinkRequiredError)
        else:
            found = []
            for tlst in self._tree_lists:
                if tlst.label.upper() == title.upper():
                    found.append(tlst)
            if len(found) == 0:
                raise self._nexus_error("Trees block with title '{}' not found".format(title), NexusReader.UndefinedBlockError)
            elif len(found) > 1:
                raise self._nexus_error("Multiple trees blocks with title '{}' defined".format(title), NexusReader.MultipleBlockWithSameTitleError)
            return found[0]

    ###########################################################################
    ## Main Stream Parse Driver

    def _parse_nexus_stream(self, stream):
        "Main file parsing driver."
        if self._nexus_tokenizer is None:
            self.create_tokenizer(stream,
                preserve_unquoted_underscores=self.preserve_underscores)
        else:
            self._nexus_tokenizer.set_stream(stream)
        token = self._nexus_tokenizer.next_token()
        if token.upper() != "#NEXUS":
            raise self._nexus_error("Expecting '#NEXUS', but found '{}'".format(token),
                    NexusReader.NotNexusFileError)
        while not self._nexus_tokenizer.is_eof():
            token = self._nexus_tokenizer.next_token_ucase()
            while token != None and token != 'BEGIN' and not self._nexus_tokenizer.is_eof():
                token = self._nexus_tokenizer.next_token_ucase()
            self._nexus_tokenizer.process_and_clear_comments_for_item(
                    self._global_annotations_target,
                    self.extract_comment_metadata)
            token = self._nexus_tokenizer.next_token_ucase()
            if token == 'TAXA':
                self._parse_taxa_block()
            elif token == 'CHARACTERS' or token == 'DATA':
                self._parse_characters_data_block()
            elif token == 'TREES':
                self._parse_trees_block()
            elif token in ['SETS', 'ASSUMPTIONS', 'CODONS']:
                if not self.exclude_chars:
                    self._nexus_tokenizer.skip_to_semicolon() # move past BEGIN command
                    link_title = None
                    block_title = None
                    while not (token == 'END' or token == 'ENDBLOCK') \
                            and not self._nexus_tokenizer.is_eof() \
                            and not token==None:
                        token = self._nexus_tokenizer.next_token_ucase()
                        if token == 'TITLE':
                            block_title = self._parse_title_statement()
                        elif token == "LINK":
                            link_title = self._parse_link_statement().get('characters')
                        elif token == 'CHARSET':
                            self._parse_charset_statement(block_title=block_title, link_title=link_title)
                        elif token == 'BEGIN':
                            raise self._nexus_error("'BEGIN' found without completion of previous block",
                                    NexusReader.IncompleteBlockError)
                    self._nexus_tokenizer.skip_to_semicolon() # move past END command
            elif token == 'BEGIN':
                raise self._nexus_error("'BEGIN' found without completion of previous block",
                        NexusReader.IncompleteBlockError)
            else:
                # unknown block
                if token is not None and self.store_ignored_blocks:
                    b = self._read_block_without_processing(token=token)
                    self._ignored_blocks.append(b)
                else:
                    token = self._consume_to_end_of_block(token)

    ###########################################################################
    ## TAXA BLOCK

    def _parse_taxa_block(self):
        token = ''
        self._nexus_tokenizer.allow_eof = False
        self._nexus_tokenizer.skip_to_semicolon() # move past BEGIN statement
        title = None
        taxon_namespace = None
        #while not (token == 'END' or token == 'ENDBLOCK') \
        #    and not self._nexus_tokenizer.is_eof() \
        #    and not token==None:
        while not (token == 'END' or token == 'ENDBLOCK'):
            token = self._nexus_tokenizer.next_token_ucase()
            if token == "TITLE":
                token = self._parse_title_statement()
                taxon_namespace = self._new_taxon_namespace(token)
            if token == 'DIMENSIONS':
                self._parse_dimensions_statement()
            if token == 'TAXLABELS':
                if taxon_namespace is None:
                    taxon_namespace = self._new_taxon_namespace()
                self._nexus_tokenizer.process_and_clear_comments_for_item(
                        self._global_annotations_target,
                        self.extract_comment_metadata)
                self._parse_taxlabels_statement(taxon_namespace)
        self._nexus_tokenizer.skip_to_semicolon() # move past END statement
        self._nexus_tokenizer.allow_eof = True

    def _get_taxon(self, taxon_namespace, label):
        if not self._file_specified_ntax or len(taxon_namespace) < self._file_specified_ntax:
            taxon = taxon_namespace.require_taxon(label=label,
                    is_case_sensitive=self.case_sensitive_taxon_labels)
        else:
            taxon = taxon_namespace.get_taxon(label=label,
                    is_case_sensitive=self.case_sensitive_taxon_labels)
        if taxon is None:
            raise self._too_many_taxa_error(taxon_namespace=taxon_namespace, label=label)
        return taxon

    def _parse_taxlabels_statement(self, taxon_namespace=None):
        """
        Processes a TAXLABELS command. Assumes that the file reader is
        positioned right after the "TAXLABELS" token in a TAXLABELS command.
        """
        if taxon_namespace is None:
            taxon_namespace = self._get_taxon_namespace()
        token = self._nexus_tokenizer.next_token()

        # Construct label lookup set
        # The get_taxon call is expensive for large taxon namespaces as it requires
        # a linear search. This causes significant performance penalties for loading
        # very large trees into an empty taxon namespace as each new taxon requires
        # a worst case search of the existing namespace before it can be inserted.
        # To alleviate this, we build a temporary one-time set of all the labels
        # in the taxon namespace. Now we can determine in constant-time whether
        # a label token corresponds to a new taxon that requires insertion,
        # or if an existing taxon can be fetched with get_taxon.
        label_set = set([])
        for taxon in taxon_namespace._taxa:
            if taxon_namespace.is_case_sensitive:
                label_set.add(taxon.label)
            else:
                label_set.add(taxon.lower_cased_label)

        while token != ';':
            label = token

            # Convert the token to the appropriate case to check against label set
            if taxon_namespace.is_case_sensitive:
                check_label = label
            else:
                check_label = label.lower()

            if check_label in label_set:
                taxon = taxon_namespace.get_taxon(label=label)
            else:
                if len(taxon_namespace) >= self._file_specified_ntax and not self.attached_taxon_namespace and not self.unconstrained_taxa_accumulation_mode:
                    raise self._too_many_taxa_error(taxon_namespace=taxon_namespace, label=label)
                taxon = taxon_namespace.new_taxon(label=label)

                # Add the new label to the label lookup set too
                if taxon_namespace.is_case_sensitive:
                    label_set.add(taxon.label)
                else:
                    label_set.add(taxon.lower_cased_label)

            token = self._nexus_tokenizer.next_token()
            self._nexus_tokenizer.process_and_clear_comments_for_item(taxon,
                    self.extract_comment_metadata)

    ###########################################################################
    ## LINK/TITLE PARSERS (How Mesquite handles multiple TAXA blocks)

    def _parse_title_statement(self):
        """
        Processes a MESQUITE 'TITLE' statement.
        Assumes current token is 'TITLE'
        """
        if self._nexus_tokenizer.cast_current_token_to_ucase() != "TITLE":
            raise self._nexus_error("Expecting 'TITLE' token, but instead found '{}'".format(self._nexus_tokenizer.cast_current_token_to_ucase()))
        title = self._nexus_tokenizer.require_next_token()
        sc = self._nexus_tokenizer.require_next_token()
        if sc != ";":
            raise self._nexus_error("Expecting ';' token, but instead found '{}'".format(sc))
        return title

    def _parse_link_statement(self):
        """
        Processes a MESQUITE 'LINK' statement.
        """
        # TODO: this is now pretty ugly
        # need to refactor with more abstraction
        links = {}
        token = self._nexus_tokenizer.next_token_ucase()
        while token != ';':
            if token == 'TAXA':
                token = self._nexus_tokenizer.next_token()
                if token != "=":
                    raise self._nexus_error("expecting '=' after link taxa")
                token = self._nexus_tokenizer.next_token()
                links['taxa'] = token
                token = self._nexus_tokenizer.next_token()
            if token == 'CHARACTERS':
                token = self._nexus_tokenizer.next_token()
                if token != "=":
                    raise self._nexus_error("expecting '=' after link characters")
                token = self._nexus_tokenizer.next_token()
                links['characters'] = token
                token = self._nexus_tokenizer.next_token()
        if token != ";":
            self._nexus_tokenizer.skip_to_semicolon()
        return links

    ###########################################################################
    ## CHARACTER/DATA BLOCK PARSERS AND SUPPORT

    def _parse_characters_data_block(self):
        token = self._nexus_tokenizer.cast_current_token_to_ucase()
        if token != "CHARACTERS" and token != "DATA":
            raise self._nexus_error("Expecting 'CHARACTERS' or 'DATA' token, but instead found '{}'".format(token))
        if self.exclude_chars:
            self._consume_to_end_of_block(self._nexus_tokenizer.current_token)
            return
        self._nexus_tokenizer.skip_to_semicolon() # move past BEGIN command
        block_title = None
        link_title = None
        self._data_type = "standard" # set as default
        while (token != 'END'
                and token != 'ENDBLOCK'
                and not self._nexus_tokenizer.is_eof()
                and not token==None):
            token = self._nexus_tokenizer.next_token_ucase()
            if token == 'TITLE':
                block_title = self._parse_title_statement()
            elif token == "LINK":
                link_title = self._parse_link_statement().get('taxa')
            elif token == 'DIMENSIONS':
                self._parse_dimensions_statement()
            elif token == 'FORMAT':
                self._parse_format_statement()
            elif token == 'MATRIX':
                self._parse_matrix_statement(block_title=block_title, link_title=link_title)
            elif token == 'BEGIN':
                raise self._nexus_error("'BEGIN' found without completion of previous block",
                        NexusReader.IncompleteBlockError)
            # token = self._nexus_tokenizer.cast_current_token_to_ucase()
        self._nexus_tokenizer.skip_to_semicolon() # move past END command

    def _build_state_alphabet(self, char_block, symbols):
        if self._gap_char and self._gap_char in symbols:
            symbols = [s for s in symbols if s != self._gap_char]
        sa = self._new_state_alphabet(
                fundamental_states=symbols,
                no_data_symbol=self._missing_char,
                gap_symbol=self._gap_char,
                case_sensitive=False)
        char_block.state_alphabets = [sa]
        char_block.default_state_alphabet = char_block.state_alphabets[0]

    def _parse_format_statement(self):
        """
        Processes a FORMAT command. Assumes that the file reader is
        positioned right after the "FORMAT" token in a FORMAT command.
        """
        token = self._nexus_tokenizer.require_next_token_ucase()
        while token != ';':
            if token == 'DATATYPE':
                token = self._nexus_tokenizer.require_next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.require_next_token_ucase()
                    if token == "DNA" or token == "NUCLEOTIDES":
                        self._data_type = "dna"
                    elif token == "RNA":
                        self._data_type = "rna"
                    elif token == "NUCLEOTIDE":
                        self._data_type = "nucleotide"
                    elif token == "PROTEIN":
                        self._data_type = "protein"
                    elif token == "CONTINUOUS":
                        self._data_type = "continuous"
                    else:
                        # defaults to STANDARD elif token == "STANDARD":
                        self._data_type = "standard"
                        self._symbols = "0123456789"
                else:
                    raise self._nexus_error("Expecting '=' after DATATYPE keyword")
                token = self._nexus_tokenizer.require_next_token_ucase()
            elif token == 'SYMBOLS':
                token = self._nexus_tokenizer.require_next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.require_next_token_ucase()
                    if token == '"':
                        self._symbols = ""
                        token = self._nexus_tokenizer.require_next_token_ucase()
                        while token != '"':
                            if token not in self._symbols:
                                self._symbols = self._symbols + token
                            token = self._nexus_tokenizer.require_next_token_ucase()
                    else:
                        raise self._nexus_error("Expecting '\"' before beginning SYMBOLS list")
                else:
                    raise self._nexus_error("Expecting '=' after SYMBOLS keyword")
                token = self._nexus_tokenizer.require_next_token_ucase()
            elif token == 'GAP':
                token = self._nexus_tokenizer.require_next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.require_next_token_ucase()
                    self._gap_char = token
                else:
                    raise self._nexus_error("Expecting '=' after GAP keyword")
                token = self._nexus_tokenizer.require_next_token_ucase()
            elif token == 'INTERLEAVE':
                token = self._nexus_tokenizer.require_next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.require_next_token_ucase()
                    if token.startswith("N"):
                        self._interleave = False
                    else:
                        self._interleave = True
                    token = self._nexus_tokenizer.require_next_token_ucase()
                else:
                    self._interleave = True
            elif token == 'MISSING':
                token = self._nexus_tokenizer.require_next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.require_next_token_ucase()
                    self._missing_char = token
                else:
                    raise self._nexus_error("Expecting '=' after MISSING keyword")
                token = self._nexus_tokenizer.require_next_token_ucase()
            elif token == 'MATCHCHAR':
                token = self._nexus_tokenizer.require_next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.require_next_token_ucase()
                    self._match_char = frozenset([token, token.lower()])
                else:
                    raise self._nexus_error("Expecting '=' after MISSING keyword")
                token = self._nexus_tokenizer.require_next_token_ucase()
            elif token == 'BEGIN':
                raise self._nexus_error("'BEGIN' found without completion of previous block",
                        NexusReader.IncompleteBlockError)
            else:
                token = self._nexus_tokenizer.require_next_token_ucase()

    def _parse_dimensions_statement(self):
        """
        Processes a DIMENSIONS command. Assumes that the file reader is
        positioned right after the "DIMENSIONS" token in a DIMENSIONS command.
        """
        token = self._nexus_tokenizer.require_next_token_ucase()
        while token != ';':
            if token == 'NTAX':
                token = self._nexus_tokenizer.require_next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.require_next_token_ucase()
                    if token.isdigit():
                        self._file_specified_ntax = int(token)
                    else:
                        raise self._nexus_error('Expecting numeric value for NTAX')
                else:
                    raise self._nexus_error("Expecting '=' after NTAX keyword")
            elif token == 'NCHAR':
                token = self._nexus_tokenizer.require_next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.require_next_token_ucase()
                    if token.isdigit():
                        self._file_specified_nchar = int(token)
                    else:
                        raise self._nexus_error("Expecting numeric value for NCHAR")
                else:
                    raise self._nexus_error("Expecting '=' after NCHAR keyword")
            elif token == 'BEGIN':
                raise self._nexus_error("'BEGIN' found without completion of previous block",
                        NexusReader.IncompleteBlockError)
            token = self._nexus_tokenizer.require_next_token_ucase()

    def _parse_matrix_statement(self, block_title=None, link_title=None):
        """
        Processes a MATRIX command. Assumes that the file reader
        is positioned right after the "MATRIX" token in a MATRIX command,
        and that NTAX and NCHAR have been specified accurately.
        """
        if not self._file_specified_ntax:
            raise self._nexus_error('NTAX must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        elif not self._file_specified_nchar:
            raise self._nexus_error('NCHAR must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        taxon_namespace = self._get_taxon_namespace(link_title)
        char_block = self._new_char_matrix(
                self._data_type,
                taxon_namespace=taxon_namespace,
                title=block_title)
        if self._data_type == "continuous":
            self._process_continuous_matrix_data(char_block)
        else:
            self._process_discrete_matrix_data(char_block)

    def _process_continuous_matrix_data(self, char_block):
        taxon_namespace = char_block.taxon_namespace
        token = self._nexus_tokenizer.next_token()
        first_sequence_defined = None
        if self._interleave:
            try:
                while token != ";" and not self._nexus_tokenizer.is_eof():
                    taxon = self._get_taxon(taxon_namespace=taxon_namespace, label=token)
                    self._read_continuous_character_values(char_block[taxon])
                    # if first_sequence_defined is None:
                    #     first_sequence_defined = char_block[taxon]
                    token = self._nexus_tokenizer.next_token()
            except NexusReader.BlockTerminatedException:
                token = self._nexus_tokenizer.next_token()
        else:
            while token != ';' and not self._nexus_tokenizer.is_eof():
                taxon = self._get_taxon(taxon_namespace=taxon_namespace, label=token)
                self._read_continuous_character_values(char_block[taxon])
                # if first_sequence_defined is None:
                #     first_sequence_defined = char_block[taxon]
                if len(char_block[taxon]) < self._file_specified_nchar:
                    raise self._nexus_error("Insufficient characters given for taxon '{}': expecting {} but only found {} ('{}')".format(taxon.label, self._file_specified_nchar, len(char_block[taxon]), char_block[taxon].symbols_as_string()))
                token = self._nexus_tokenizer.next_token()
        # if self._interleave:
        #     raise NotImplementedError("Continuous interleaved characters in NEXUS schema not yet supported")
        # taxon_namespace = char_block.taxon_namespace
        # token = self._nexus_tokenizer.next_token()
        # while token != ';' and not self._nexus_tokenizer.is_eof():
        #     taxon = self._get_taxon(taxon_namespace=taxon_namespace, label=token)
        #     while len(char_block[taxon]) < self._file_specified_nchar and not self._nexus_tokenizer.is_eof():
        #         # char_group = self._nexus_tokenizer.next_token(ignore_punctuation="-+")
        #         char_group = self._nexus_tokenizer.next_token()
        #         char_block[taxon].append(dataobject.CharacterDataCell(value=float(char_group)))
        #     if len(char_block[taxon]) < self._file_specified_nchar:
        #         raise self._nexus_error("Insufficient characters given for taxon '%s': expecting %d but only found %d ('%s')" \
        #             % (taxon.label, self._file_specified_nchar, len(char_block[taxon]), char_block[taxon].symbols_as_string()))
        #     token = self._nexus_tokenizer.next_token()

    def _process_discrete_matrix_data(self, char_block):
        if self._data_type == "standard":
            self._build_state_alphabet(char_block, self._symbols)
        taxon_namespace = char_block.taxon_namespace
        token = self._nexus_tokenizer.next_token()
        state_alphabet = char_block.default_state_alphabet
        first_sequence_defined = None
        if self._interleave:
            try:
                while token != ";" and not self._nexus_tokenizer.is_eof():
                    taxon = self._get_taxon(taxon_namespace=taxon_namespace, label=token)
                    self._read_character_states(char_block[taxon], state_alphabet, first_sequence_defined)
                    if first_sequence_defined is None:
                        first_sequence_defined = char_block[taxon]
                    token = self._nexus_tokenizer.next_token()
            except NexusReader.BlockTerminatedException:
                token = self._nexus_tokenizer.next_token()
        else:
            while token != ';' and not self._nexus_tokenizer.is_eof():
                taxon = self._get_taxon(taxon_namespace=taxon_namespace, label=token)
                self._read_character_states(char_block[taxon], state_alphabet, first_sequence_defined)
                if first_sequence_defined is None:
                    first_sequence_defined = char_block[taxon]
                if len(char_block[taxon]) < self._file_specified_nchar:
                    raise self._nexus_error("Insufficient characters given for taxon '%s': expecting %d but only found %d ('%s')" \
                        % (taxon.label, self._file_specified_nchar, len(char_block[taxon]), char_block[taxon].symbols_as_string()))
                token = self._nexus_tokenizer.next_token()

    def _get_state_for_multistate_tokens(self,
            state_char_seq,
            multistate_type,
            state_alphabet):
        try:
            state = state_alphabet.match_state(state_char_seq,
                    state_denomination=multistate_type)
        except KeyError:
            try:
                if multistate_type == state_alphabet.AMBIGUOUS_STATE:
                    sae = state_alphabet.new_ambiguous_state(
                            symbol=None,
                            member_state_symbols=state_char_seq)
                else:
                    sae = state_alphabet.new_polymorphic_state(
                            symbol=None,
                            member_state_symbols=state_char_seq)
            except KeyError:
                raise self._nexus_error("Unrecognized state symbols encountered in multistate sequence: '{}'".format(state_char_seq))
            else:
                return sae
        else:
            return state

    ###########################################################################
    ## TREE / TREE BLOCK PARSERS

    def _parse_tree_statement(self, tree_factory, taxon_symbol_mapper):
        """
        Processes a TREE command. Assumes that the file reader is
        positioned right after the "TREE" token in a TREE command.
        Calls on the NewickStatementParser of the trees module.
        """
        token = self._nexus_tokenizer.next_token()
        if token == '*':
            token = self._nexus_tokenizer.next_token()
        tree_name = token
        token = self._nexus_tokenizer.next_token()
        pre_tree_comments = self._nexus_tokenizer.pull_captured_comments()
        if token != '=':
            raise self._nexus_error("Expecting '=' in definition of Tree '%s' but found '%s'" % (tree_name, token))
        tree_comments = self._nexus_tokenizer.pull_captured_comments()
        # advance to '('; comments will be processed by newick reader
        self._nexus_tokenizer.next_token()
        tree = self._build_tree_from_newick_tree_string(tree_factory, taxon_symbol_mapper)
        tree.label = tree_name
        nexusprocessing.process_comments_for_item(tree, pre_tree_comments, self.extract_comment_metadata)
        nexusprocessing.process_comments_for_item(tree, tree_comments, self.extract_comment_metadata)
        # if self.extract_comment_metadata:
        #     annotations = nexustokenizer.parse_comment_metadata(tree_comments)
        #     for annote in annotations:
        #         tree.annotations.add(annote)
        #     if pre_tree_metadata_comments:
        #         pre_tree_annotations = nexustokenizer.parse_comment_metadata(pre_tree_metadata_comments)
        #         for annote in pre_annotations:
        #             tree.annotations.add(annote)
        # if tree_comments is not None and len(tree_comments) > 0:
        #     tree.comments.extend(tree_comments)
        # if self._nexus_tokenizer.current_token != ';':
        #     self._nexus_tokenizer.skip_to_semicolon()
        return tree

    def _build_tree_from_newick_tree_string(self, tree_factory, taxon_symbol_mapper):
        tree = self.newick_reader._parse_tree_statement(
                nexus_tokenizer=self._nexus_tokenizer,
                tree_factory=tree_factory,
                taxon_symbol_map_fn=taxon_symbol_mapper.require_taxon_for_symbol)
        return tree

    def _parse_translate_statement(self, taxon_namespace, taxon_symbol_mapper=None):
        """
        Processes a TRANSLATE command. Assumes that the file reader is
        positioned right after the "TRANSLATE" token in a TRANSLATE command.
        """
        token = self._nexus_tokenizer.current_token
        if taxon_symbol_mapper is None:
            taxon_symbol_mapper = self._get_taxon_symbol_mapper(taxon_namespace=taxon_namespace)
        else:
            assert taxon_symbol_mapper.taxon_namespace is taxon_namespace
        if self._file_specified_ntax is None:
            # Not yet parsed TAXA block: NEXUS file without TAXA block
            # Badly-formed NEXUS file, yet widely-found in the wild
            # Override namespace modification lock
            taxon_namespace.is_mutable = True
        while True:
            translation_token = self._nexus_tokenizer.next_token()
            if translation_token == ";" and not self._nexus_tokenizer.is_token_quoted:
                raise self._nexus_error("Expecting translation token but found ';' instead")
            translation_label = self._nexus_tokenizer.next_token()
            try:
                taxon = taxon_namespace.require_taxon(label=translation_label)
            except error.ImmutableTaxonNamespaceError:
                exc = self._undefined_taxon_error(taxon_namespace=taxon_namespace, label=translation_label)
                exc.__context__ = None # Python 3.0, 3.1, 3.2
                exc.__cause__ = None # Python 3.3, 3.4
                raise exc
            taxon_symbol_mapper.add_translate_token(translation_token, taxon)
            token = self._nexus_tokenizer.next_token() # ","
            if (not token) or (token == ';'):
                break
            if token != ',':
                raise self._nexus_error("Expecting ',' in TRANSLATE statement after definition for %s = '%s', but found '%s' instead." % (translation_token, translation_label, token))
        return taxon_symbol_mapper

    def _parse_trees_block(self):
        """
        Expectations:
            - current token: "TREES" [part of "BEGIN TREES"]
        """
        token = self._nexus_tokenizer.cast_current_token_to_ucase()
        if token != "TREES":
            raise self._nexus_error("Expecting 'TREES' token, but instead found '{}'".format(token))
        if self.exclude_trees:
            self._consume_to_end_of_block(self._nexus_tokenizer.current_token)
            return
        self._nexus_tokenizer.skip_to_semicolon() # move past "BEGIN TREES" command
        link_title = None
        taxon_namespace = None
        taxon_symbol_mapper = None
        trees_block = None
        block_title = None
        # while ((not self._nexus_tokenizer.is_eof())
        #         and self._nexus_tokenizer.current_token is not None
        #         and self._nexus_tokenixer.current_token != 'END'
        #         and self._nexus_tokenixer.current_token != 'ENDBLOCK'):
        while ((not self._nexus_tokenizer.is_eof())
                and token is not None
                and token != 'END'
                and token != 'ENDBLOCK'):
            token = self._nexus_tokenizer.next_token_ucase()
            if token == 'LINK':
                link_title = self._parse_link_statement().get("taxa")
            elif token == 'TITLE':
                block_title = self._parse_title_statement()
                token = "" # clear; repopulate at start of loop
            elif token == 'TRANSLATE':
                if taxon_namespace is None:
                    taxon_namespace = self._get_taxon_namespace(link_title)
                taxon_symbol_mapper = self._parse_translate_statement(taxon_namespace)
                token = "" # clear; repopulate at start of loop
            elif token == 'TREE':
                if taxon_namespace is None:
                    taxon_namespace = self._get_taxon_namespace(link_title)
                if taxon_symbol_mapper is None:
                    taxon_symbol_mapper = self._get_taxon_symbol_mapper(taxon_namespace=taxon_namespace)
                pre_tree_comments = self._nexus_tokenizer.pull_captured_comments()
                if trees_block is None:
                    trees_block = self._new_tree_list(taxon_namespace=taxon_namespace, title=block_title)
                # All comments leading up to the first 'TREE' statement assumed
                # to belong to the TreeList corresponding to the TREES block
                nexusprocessing.process_comments_for_item(
                        trees_block,
                        pre_tree_comments,
                        self.extract_comment_metadata)
                tree_factory = trees_block.new_tree
                while True:
                    ## After the following, the current token
                    ## will be the token immediately following
                    ## the terminating semi-colon of a tree
                    ## statement. Typically, this will be
                    ## 'TREE' if there is another tree, or
                    ## 'END'/'ENDBLOCK'.
                    tree = self._parse_tree_statement(
                            tree_factory=tree_factory,
                            taxon_symbol_mapper=taxon_symbol_mapper)
                    if self._nexus_tokenizer.is_eof() or not self._nexus_tokenizer.current_token:
                        break
                    if self._nexus_tokenizer.cast_current_token_to_ucase() != "TREE":
                        token = self._nexus_tokenizer.current_token
                        break
            elif token == 'BEGIN':
                raise self._nexus_error("'BEGIN' found without completion of previous block",
                        NexusReader.IncompleteBlockError)
        self._nexus_tokenizer.skip_to_semicolon() # move past END command

    def _parse_charset_statement(self, block_title=None, link_title=None):
        """
        Parses a character set description. Assumes token stream is positioned right after 'charset' command.
        """
        char_matrix = self._get_char_matrix(title=link_title)
        keyword = self._nexus_tokenizer.current_token
        token = self._nexus_tokenizer.next_token()
        if self._nexus_tokenizer.is_eof() or not token:
            raise self._nexus_error('Unexpected end of file or null token')
        else:
            if not token:
                raise self._nexus_error("Unexpected end of file or null token")
            else:
                charset_name = token
                token = self._nexus_tokenizer.next_token()
                if not token:
                    raise self._nexus_error("Unexpected end of file or null token")
                elif token != '=':
                    raise self._nexus_error('Expecting "=" after character set name "%s", but instead found "%s"' % (charset_name, token))
                else:
                    positions = self._parse_positions(adjust_to_zero_based=True)
                char_matrix.new_character_subset(charset_name, positions)

    def _parse_positions(self, adjust_to_zero_based=True, verify=True):
        """
        Parses a character position list. Expects next character read to be the first item in a position list.
        """
        positions = []
        # hyphens_as_tokens = self._nexus_tokenizer.hyphens_as_tokens
        # self._nexus_tokenizer.hyphens_as_tokens = True
        self._nexus_tokenizer.set_hyphens_as_captured_delimiters(True)
        token = self._nexus_tokenizer.next_token()
        max_positions = self._file_specified_nchar

        if self._nexus_tokenizer.is_eof() or not token:
            raise self._nexus_error('Unexpected end of file or null token')

        while token != ';' and token != ',' and not self._nexus_tokenizer.is_eof():
            if not token:
                break
            if token.upper() == 'ALL':
                positions = range(1, max_positions + 1)
                break
            elif token.isdigit():
                start = int(token)
                token = self._nexus_tokenizer.next_token()
                if token:
                    if token == ',' or token.isdigit() or token == ';':
                        positions.append(start)
                    elif token == '-':
                        token = self._nexus_tokenizer.next_token()
                        if token:
                            if token.isdigit() or token == '.':
                                if token == '.':
                                    end = max_positions
                                    #token = self._nexus_tokenizer.next_token()
                                else:
                                    end = int(token)
                                    #token = self._nexus_tokenizer.next_token()
                                token = self._nexus_tokenizer.next_token()
                                if token:
                                    if token == '\\' or token == '/': # (NEXUS standard only accepts '\')
                                        token = self._nexus_tokenizer.next_token()
                                        if token:
                                            if token.isdigit():
                                                step = int(token)
                                                #token = self._nexus_tokenizer.next_token()
                                            else:
                                                raise self._nexus_error('Expecting digit but found "%s".' % (token))
                                        else:
                                            raise self._nexus_error(r'Expecting other tokens after "\", but no more found.')
                                        token = self._nexus_tokenizer.next_token()
                                    else:
                                        step = 1
                                else:
                                    step = 1
                                for q in range(start, end+1, step):
                                    if q <= max_positions:
                                        positions.append(q)
                            else:
                                raise self._nexus_error('Expecting digit or ".", but found "%s".' % (token))
                        else:
                            raise self._nexus_error('Expecting other tokens after "-", but no more found.')
                    else:
                        raise self._nexus_error('Expecting digit or "all", but found "%s".' % (token))
                else:
                    positions.append(start)
        self._nexus_tokenizer.set_hyphens_as_captured_delimiters(False)
        positions = list(set(positions))
        positions.sort()
        if verify:
            for position in positions:
                if position > max_positions:
                    raise self._nexus_error("Specified position %d, but maximum position is %d" % (position, max_positions))
        if adjust_to_zero_based:
            positions = [position - 1 for position in positions]
        return positions # make unique and return

    def _consume_to_end_of_block(self, token=None):
        if token:
            token = token.upper()
        else:
            token = "DUMMY"
        while not (token == 'END' or token == 'ENDBLOCK') \
                and not self._nexus_tokenizer.is_eof() \
                and not token==None:
            self._nexus_tokenizer.skip_to_semicolon()
            token = self._nexus_tokenizer.next_token_ucase()
        return token

    def _read_block_without_processing(self, token=None):
        # used for unknown blocks we want to save
        # NOT (really) TESTED
        # Everybody else except Jeet: (REALLY) DO NOT USE!
        # Jeet: SORTA DO NOT USE WITHOUT MORE TESTING
        if token:
            token = token.upper()
        block = ["BEGIN", token]
        old_uncaptured_delimiters = self._nexus_tokenizer.uncaptured_delimiters
        old_captured_delimiters = self._nexus_tokenizer.captured_delimiters
        to_switch = "\n\r"
        for ch in to_switch:
            self._nexus_tokenizer.uncaptured_delimiters.discard(ch)
            self._nexus_tokenizer.captured_delimiters.add(ch)
        while not (token == 'END' or token == 'ENDBLOCK') \
                and not self._nexus_tokenizer.is_eof() \
                and not token==None:
            token = self._nexus_tokenizer.require_next_token()
            uctoken = token.upper()
            if uctoken == "END" or uctoken == "ENDBLOCK":
                token = uctoken
            block.append(token)
        self._nexus_tokenizer.uncaptured_delimiters = old_uncaptured_delimiters
        self._nexus_tokenizer.captured_delimiters = old_captured_delimiters
        self._nexus_tokenizer.skip_to_semicolon() # move past end
        block.append(";")
        return " ".join(block)

    def _read_character_states(self,
            character_data_vector,
            state_alphabet,
            first_sequence_defined,
            ):
        """
        Reads character sequence data substatement until the number of
        character states read is equal to ``self._file_specified_nchar`` (with
        multi-state characters, such as '(AG)' counting as a single
        state) or, if ``self._interleave`` is |True|, until an EOL is
        reached.

        Given a sequence of characters, with ambiguities denoted by
        `{<STATES>}`, this returns a list of state alphabet elements.

        For example, the following sequence:

            "ACTG(AC)GGT(CGG)(CG)GG"

        will result in a list such as:

            [<A>, <C>, <T>, <G>, <AC>, <G>, <G>, <T>, <CGG>, <CG>, <G>, <G>]

        where `<.>` is a StateIdentity object with the characters within the
        brackets as symbol(s).

        """
        if self._interleave:
            self._nexus_tokenizer.set_capture_eol(True)
        states_to_add = []
        while len(character_data_vector) + len(states_to_add) < self._file_specified_nchar:
            token = self._nexus_tokenizer.require_next_token()
            if token == "{" or token == "(":
                if token == "{":
                    # multistate_type = dataobject.StateIdentity.AMBIGUOUS_STATE
                    multistate_type = state_alphabet.AMBIGUOUS_STATE
                    closing_token = "}"
                else:
                    # multistate_type = dataobject.StateIdentity.POLYMORPHIC_STATE
                    multistate_type = state_alphabet.POLYMORPHIC_STATE
                    closing_token = ")"
                multistate_tokens = []
                while True:
                    token = self._nexus_tokenizer.require_next_token()
                    if token == closing_token:
                        break
                    multistate_tokens.append(token)
                c = "".join(multistate_tokens)
                state = self._get_state_for_multistate_tokens(c, multistate_type, state_alphabet)
                if len(character_data_vector) + len(states_to_add) == self._file_specified_nchar:
                    raise self._too_many_characters_error(c)
                states_to_add.append(state)
            elif token == "\r" or token == "\n":
                if self._interleave:
                    break
            elif token == ";":
                raise NexusReader.BlockTerminatedException
            else:
                for c in token:
                    if c in self._match_char:
                        try:
                            state = first_sequence_defined[len(character_data_vector) + len(states_to_add)]
                        except TypeError:
                            exc = self._nexus_error("Cannot dereference MATCHCHAR '{}' on first sequence".format(c), NexusReader.NexusReaderError)
                            exc.__context__ = None # Python 3.0, 3.1, 3.2
                            exc.__cause__ = None # Python 3.3, 3.4
                            raise exc
                        except IndexError:
                            exc = self._nexus_error("Cannot dereference MATCHCHAR '{}': current position ({}) exceeds length of first sequence ({})".format(c,
                                    len(character_data_vector) + len(states_to_add) + 1,
                                    len(first_sequence_defined),
                                    NexusReader.NexusReaderError))
                            exc.__context__ = None # Python 3.0, 3.1, 3.2
                            exc.__cause__ = None # Python 3.3, 3.4
                            raise exc
                    else:
                        try:
                            state = state_alphabet.full_symbol_state_map[c]
                        except KeyError:
                            exc = self._nexus_error("Unrecognized character state symbol for state alphabet '{}' ({}) : '{}'".format(
                                        state_alphabet.label,
                                        state_alphabet.__class__.__name__,
                                        c),
                                        NexusReader.InvalidCharacterStateSymbolError)
                            exc.__context__ = None # Python 3.0, 3.1, 3.2
                            exc.__cause__ = None # Python 3.3, 3.4
                            raise exc
                    if len(character_data_vector) + len(states_to_add) == self._file_specified_nchar:
                        raise self._too_many_characters_error(c)
                    states_to_add.append(state)
        if self._interleave:
            self._nexus_tokenizer.set_capture_eol(False)
        character_data_vector.extend(states_to_add)
        return character_data_vector

    def _read_continuous_character_values(self,
            character_data_vector,
            datatype=float,
            ):
        """
        Reads character sequence data substatement until the number of
        character states read is equal to ``self._file_specified_nchar`` (with
        multi-state characters, such as '(AG)' counting as a single
        state) or, if ``self._interleave`` is |True|, until an EOL is
        reached.
        """
        if self._interleave:
            self._nexus_tokenizer.set_capture_eol(True)
        while len(character_data_vector) < self._file_specified_nchar:
            token = self._nexus_tokenizer.require_next_token()
            if token == "\r" or token == "\n":
                if self._interleave:
                    break
            elif token == ";":
                raise NexusReader.BlockTerminatedException
            else:
                try:
                    state = float(token)
                except ValueError:
                    exc = self._nexus_error("Invalid value for continuous character type: '{invalid_value}'".format(datatype=datatype, invalid_value=token),
                                NexusReader.InvalidContinuousCharacterValueError)
                    exc.__context__ = None # Python 3.0, 3.1, 3.2
                    exc.__cause__ = None # Python 3.3, 3.4
                    raise exc
                    # if c in self._match_char:
                    #     try:
                    #         state = first_sequence_defined[len(character_data_vector)]
                    #     except TypeError:
                    #         exc = self._nexus_error("Cannot dereference MATCHCHAR '{}' on first sequence".format(c), NexusReader.NexusReaderError)
                    #         exc.__context__ = None # Python 3.0, 3.1, 3.2
                    #         exc.__cause__ = None # Python 3.3, 3.4
                    #         raise exc
                    #     except IndexError:
                    #         exc = self._nexus_error("Cannot dereference MATCHCHAR '{}': current position ({}) exceeds length of first sequence ({})".format(c,
                    #                 len(character_data_vector)+1,
                    #                 len(first_sequence_defined),
                    #                 NexusReader.NexusReaderError))
                    #         exc.__context__ = None # Python 3.0, 3.1, 3.2
                    #         exc.__cause__ = None # Python 3.3, 3.4
                    #         raise exc
                    # else:
                    #     try:
                    #         state = state_alphabet.full_symbol_state_map[c]
                    #     except KeyError:
                    #         exc = self._nexus_error("Unrecognized character state symbol for state alphabet '{}' ({}) : '{}'".format(
                    #                     state_alphabet.label,
                    #                     state_alphabet.__class__.__name__,
                    #                     c),
                    #                     NexusReader.InvalidCharacterStateSymbolError)
                    #         exc.__context__ = None # Python 3.0, 3.1, 3.2
                    #         exc.__cause__ = None # Python 3.3, 3.4
                    #         raise exc
                if len(character_data_vector) == self._file_specified_nchar:
                    raise self._too_many_characters_error(token)
                character_data_vector.append(state)
        if self._interleave:
            self._nexus_tokenizer.set_capture_eol(False)
        return character_data_vector
