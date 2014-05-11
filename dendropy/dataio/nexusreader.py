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
Implementation of NEXUS-schema data reader.
"""

import re
import collections
from dendropy.utility import error
from dendropy.utility import text
from dendropy.dataio import ioservice
from dendropy.dataio import nexusprocessing
from dendropy.dataio import newickreader

###############################################################################
## NexusReader

class NexusReader(ioservice.DataReader):
    "Encapsulates loading and parsing of a NEXUS schema file."

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

    class NexusReaderNotNexusFileError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NexusReaderLinkRequiredError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NexusReaderNoCharacterBlocksFoundError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NexusReaderUndefinedBlockError(NexusReaderError):
        def __init__(self, message,
                line_num=None,
                col_num=None,
                stream=None):
            NexusReader.NexusReaderError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class NexusReaderMultipleBlockWithSameTitleError(NexusReaderError):
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

        # Fallowing are NEXUS-parsing specific (i.e., not used by NEWICK parsers)
        self.exclude_chars = kwargs.pop("exclude_chars", False)
        self.exclude_trees = kwargs.pop("exclude_trees", False)
        self.data_type = kwargs.pop("data_type", "standard")

        # The following are used by NewickReader in addition to NexusReader,
        # or have different defaults. So they are extracted/set here and
        # then forwarded on ...
        self.extract_comment_metadata = kwargs.pop('extract_comment_metadata', True)
        kwargs["extract_comment_metadata"] = self.extract_comment_metadata
        self.preserve_underscores = kwargs.get('preserve_underscores', True)
        self.case_sensitive_taxon_labels = kwargs.get('case_sensitive_taxon_labels', False)

        # Create newick handler
        self.newick_reader = newickreader.NewickReader(**kwargs)

        # Set up parsing meta-variables
        self._interleave = False
        self._symbols = ""
        self._gap_char = '-'
        self._missing_char = '?'
        self._match_char = '.'
        self._file_specified_ntax = None
        self._file_specified_nchar = None
        self._nexus_tokenizer = None
        self._taxon_namespace_factory = None
        self._tree_list_factory = None
        self._char_matrix_factory = None
        self._global_annotations_target = None
        self._tree_translate_dict = {}
        self._taxon_namespaces = []
        self._char_matrices = []
        self._tree_lists = []
        self._product = None

    ###########################################################################
    ## Reader Implementation

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            global_annotations_target=None):
        """
        Instantiates and returns a DataSet object based on the
        NEXUS-formatted contents given in the file-like object `stream`.
        """
        self._nexus_tokenizer = nexusprocessing.NexusTokenizer(stream)
        self._taxon_namespace_factory = taxon_namespace_factory
        self._tree_list_factory = tree_list_factory
        self._char_matrix_factory = char_matrix_factory
        self._global_annotations_target = global_annotations_target
        self._parse_nexus_file()
        self._product = self.Product(
                taxon_namespaces=self._taxon_namespaces,
                tree_lists=self._tree_lists,
                char_matrices=self._char_matrices)
        return self._product

    ###########################################################################
    ## Book-keeping Control

    def _nexus_error(self, message, error_type=None):
        if error_type is None:
            error_type = NewickReader.NewickReaderError
        e = error_type(
                message=message,
                line_num=self._nexus_tokenizer.token_line_num,
                col_num=self._nexus_tokenizer.token_column_num,
                stream=self._nexus_tokenizer.src)
        return e

    def too_many_taxa_error(self, taxon_namespace, label):
        """
        Returns an exception object parameterized with line and
        column number values.
        """
        return self._nexus_tokenizer.too_many_taxa_error(taxon_namespace=taxon_namespace,
                max_taxa=self._file_specified_ntax,
                label=label)

    ###########################################################################
    ## Data Management

    def _new_taxon_namespace(self, title=None):
        # if self.attached_taxon_namespace is not None:
        #     return self.attached_taxon_namespace
        taxon_namespace = self._taxon_namespace_factory(label=title)
        self._taxon_namespaces.append(taxon_namespace)
        return taxon_namespace

    def _get_taxon_namespace(self, title=None):
        # if self.attached_taxon_namespace is not None:
        #     return self.attached_taxon_namespace
        if title is None:
            if len(self._taxon_namespaces) == 0:
                return self._new_taxon_namespace(title=title)
            elif len(self._taxon_namespaces) == 1:
                return self._taxon_namespaces[0]
            else:
                raise self._nexus_error("Multiple taxa blocks defined: require 'LINK' statement", NexusReader.NexusReaderLinkRequiredError)
        else:
            found = []
            for tns in self._taxon_namespaces:
                if tns.label.upper() == title.upper():
                    found.append(tns)
            if len(found) == 0:
                raise self._nexus_error("Taxa block with title '{}' not found".format(title), NexusReader.NexusReaderUndefinedBlockError)
            elif len(found) > 1:
                raise self._nexus_error("Multiple taxa blocks with title '{}' defined".format(title), NexusReader.NexusReaderMultipleBlockWithSameTitleError)
            return found[0]

    def _new_char_matrix(self, data_type, taxon_namespace, title=None):
        char_matrix = self._char_matrix_factory(
                data_type=data_type,
                taxon_namespace=taxon_namespace,
                label=title)
        self._char_matrices.append(char_matrix)
        return char_matrix

    def _get_char_matrix(self, title=None):
        if title is None:
            if len(self._char_matrices) == 1:
                return self.char_matrices[0]
            elif len(self._char_matrices) == 0:
                raise self._nexus_error("No character matrices defined", NexusReader.NexusReaderNoCharacterBlocksFoundError)
            else:
                raise self._nexus_error("Multiple character matrices defined: require 'LINK' statement", NexusReader.NexusReaderLinkRequiredError)
        else:
            found = []
            for cm in self._char_matrices:
                if cm.label.upper() == title.upper():
                    found.append(cm)
            if len(found) == 0:
                raise self._nexus_error("Character block with title '{}' not found".format(title), NexusReader.NexusReaderUndefinedBlockError)
            elif len(found) > 1:
                raise self._nexus_error("Multiple character blocks with title '{}' defined".format(title), NexusReader.NexusReaderMultipleBlockWithSameTitleError)
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
                raise self._nexus_error("No tree blocks defined", NexusReader.NexusReaderNoCharacterBlocksFoundError)
            else:
                raise self._nexus_error("Multiple tree blocks defined: require 'LINK' statement", NexusReader.NexusReaderLinkRequiredError)
        else:
            found = []
            for tlst in self._tree_lists:
                if tlst.label.upper() == title.upper():
                    found.append(tlst)
            if len(found) == 0:
                raise self._nexus_error("Character block with title '{}' not found".format(title), NexusReader.NexusReaderUndefinedBlockError)
            elif len(found) > 1:
                raise self._nexus_error("Multiple character blocks with title '{}' defined".format(title), NexusReader.NexusReaderMultipleBlockWithSameTitleError)
            return found[0]

    # def tree_source_iter(self, stream):
    #     """
    #     Iterates over a NEXUS-formatted source of trees.
    #     Only trees will be returned, and any and all character data will
    #     be skipped. The iterator will span over multiple tree blocks,
    #     but, because our NEXUS data model implementation currently does
    #     not recognize multiple taxon collection definnitions, taxa in
    #     those tree blocks will be aggregated into the same `TaxonSet` (a
    #     new one created, or the one passed to this method via the
    #     `taxon_namespace` argument). This behavior is similar to how multiple
    #     tree blocks are handled by a full NEXUS data file read.
    #     """
    #     self._reset()
    #     if self.dataset is None:
    #         self.dataset = dataobject.DataSet()
    #     self._nexus_tokenizer = nexustokenizer.NexusTokenizer(stream,
    #             preserve_underscores=self.preserve_underscores,
    #             hyphens_as_tokens=self.hyphens_as_tokens,
    #             extract_comment_metadata=self.extract_comment_metadata)
    #     token = self._nexus_tokenizer.next_token_ucase()
    #     if token.upper() != "#NEXUS":
    #         raise self._nexus_error("Expecting '#NEXUS', but found '%s'" % token)
    #     while not self._nexus_tokenizer.is_eof():
    #         token = self._nexus_tokenizer.next_token_ucase()
    #         while token != None and token != 'BEGIN' and not self._nexus_tokenizer.is_eof():
    #             token = self._nexus_tokenizer.next_token_ucase()
    #         token = self._nexus_tokenizer.next_token_ucase()
    #         if token == 'TAXA':
    #             self._parse_taxa_block()
    #         elif token == 'TREES':
    #             self._nexus_tokenizer.skip_to_semicolon() # move past BEGIN command
    #             link_title = None
    #             taxon_namespace = None
    #             self._tree_translate_dict.clear()
    #             while not (token == 'END' or token == 'ENDBLOCK') \
    #                     and not self._nexus_tokenizer.is_eof() \
    #                     and not token==None:
    #                 token = self._nexus_tokenizer.next_token_ucase()
    #                 if token == 'LINK':
    #                     link_title = self._parse_link_statement().get('taxa')
    #                 if token == 'TRANSLATE':
    #                     if not taxon_namespace:
    #                         taxon_namespace = self._get_taxon_namespace(link_title)
    #                         self._prepopulate_translate_dict(taxon_namespace)
    #                     self._parse_translate_statement(taxon_namespace)
    #                 if token == 'TREE':
    #                     if not taxon_namespace:
    #                         taxon_namespace = self._get_taxon_namespace(link_title)
    #                         self._prepopulate_translate_dict(taxon_namespace)
    #                     tree = self._parse_tree_statement(taxon_namespace)
    #                     yield tree
    #             self._nexus_tokenizer.skip_to_semicolon() # move past END command
    #         else:
    #             # unknown block
    #             while not (token == 'END' or token == 'ENDBLOCK') \
    #                 and not self._nexus_tokenizer.is_eof() \
    #                 and not token==None:
    #                 self._nexus_tokenizer.skip_to_semicolon()
    #                 token = self._nexus_tokenizer.next_token_ucase()
    #     self._reset()


    ###########################################################################
    ## Main Stream Parse Driver

    def _parse_nexus_file(self):
        "Main file parsing driver."
        token = self._nexus_tokenizer.next_token()
        if token.upper() != "#NEXUS":
            raise self._nexus_error("Expecting '#NEXUS', but found '{}'".format(token),
                    NexusReader.NexusReaderNotNexusFileError)
        else:
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
                elif token == 'CHARACTERS':
                    if not self.exclude_chars:
                        self._nexus_tokenizer.skip_to_semicolon() # move past BEGIN command
                        link_title = None
                        block_title = None
                        while not (token == 'END' or token == 'ENDBLOCK') \
                                and not self._nexus_tokenizer.is_eof() \
                                and not token==None:
                            token = self._nexus_tokenizer.next_token_ucase()
                            if token == 'TITLE':
                                token = self._nexus_tokenizer.next_token()
                                block_title = token
                            if token == "LINK":
                                link_title = self._parse_link_statement().get('taxa')
                            if token == 'DIMENSIONS':
                                self._parse_dimensions_statement()
                            if token == 'FORMAT':
                                self._parse_format_statement()
                            if token == 'MATRIX':
                                self._parse_matrix_statement(block_title=block_title, link_title=link_title)
                        self._nexus_tokenizer.skip_to_semicolon() # move past END command
                    else:
                        token = self._consume_to_end_of_block(token)
                elif token == 'DATA':
                    if not self.exclude_chars:
                        self._nexus_tokenizer.skip_to_semicolon() # move past BEGIN command
                        block_title = None
                        link_title = None
                        while not (token == 'END' or token == 'ENDBLOCK') \
                                and not self._nexus_tokenizer.is_eof() \
                                and not token==None:
                            token = self._nexus_tokenizer.next_token_ucase()
                            if token == 'TITLE':
                                token = self._nexus_tokenizer.next_token()
                                block_title = token
                            if token == "LINK":
                                link_title = self._parse_link_statement().get('taxa')
                            if token == 'DIMENSIONS':
                                self._parse_dimensions_statement()
                            if token == 'FORMAT':
                                self._parse_format_statement()
                            if token == 'MATRIX':
                                self._parse_matrix_statement(block_title=block_title, link_title=link_title)
                        self._nexus_tokenizer.skip_to_semicolon() # move past END command
                    else:
                        token = self._consume_to_end_of_block(token)
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
                                token = self._nexus_tokenizer.next_token()
                                block_title = token
                            if token == "LINK":
                                link_title = self._parse_link_statement().get('characters')
                            if token == 'CHARSET':
                                self._parse_charset_statement(block_title=block_title, link_title=link_title)
                        self._nexus_tokenizer.skip_to_semicolon() # move past END command
                else:
                    # unknown block
                    token = self._consume_to_end_of_block(token)

    ###########################################################################
    ## TAXA BLOCK / LINK PARSERS

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
                token = self._nexus_tokenizer.next_token()
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
            taxon = taxon_namespace.require_taxon(label=label, case_sensitive=self.case_sensitive_taxon_labels)
        else:
            taxon = taxon_namespace.get_taxon(label=label, case_sensitive=self.case_sensitive_taxon_labels)
        if taxon is None:
            raise self.too_many_taxa_error(taxon_namespace=taxon_namespace, label=label)
        return taxon

    def _parse_taxlabels_statement(self, taxon_namespace=None):
        """
        Processes a TAXLABELS command. Assumes that the file reader is
        positioned right after the "TAXLABELS" token in a TAXLABELS command.
        """
        if taxon_namespace is None:
            taxon_namespace = self._get_taxon_namespace()
        token = self._nexus_tokenizer.next_token()
        while token != ';':
            label = token
            # if taxon_namespace.has_taxon(label=label):
            #     pass
            # elif len(taxon_namespace) >= self._file_specified_ntax and not self.attached_taxon_namespace:
            #     raise self.too_many_taxa_error(taxon_namespace=taxon_namespace, label=label)
            # else:
            #     taxon_namespace.require_taxon(label=label)
            if len(taxon_namespace) >= self._file_specified_ntax and not self.attached_taxon_namespace:
                raise self.too_many_taxa_error(taxon_namespace=taxon_namespace, label=label)
            taxon = taxon_namespace.require_taxon(label=label)
            token = self._nexus_tokenizer.next_token()
            self._nexus_tokenizer.process_and_clear_comments_for_item(taxon,
                    self.extract_comment_metadata)

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

    def _build_state_alphabet(self, char_block, symbols):
        sa = dataobject.get_state_alphabet_from_symbols(symbols,
                gap_symbol=self._gap_char,
                missing_symbol=self._missing_char)
        char_block.state_alphabets = [sa]
        char_block.default_state_alphabet = char_block.state_alphabets[0]

    def _parse_format_statement(self):
        """
        Processes a FORMAT command. Assumes that the file reader is
        positioned right after the "FORMAT" token in a FORMAT command.
        """
        token = self._nexus_tokenizer.next_token_ucase()
        while token != ';':
            if token == 'DATATYPE':
                token = self._nexus_tokenizer.next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.next_token_ucase()
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
                        self._symbols = "01"
                else:
                    raise self._nexus_error("Expecting '=' after DATATYPE keyword")
                token = self._nexus_tokenizer.next_token_ucase()
            elif token == 'SYMBOLS':
                token = self._nexus_tokenizer.next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.next_token_ucase()
                    if token == '"':
                        self._symbols = ""
                        token = self._nexus_tokenizer.next_token_ucase()
                        while token != '"':
                            if token not in self._symbols:
                                self._symbols = self._symbols + token
                            token = self._nexus_tokenizer.next_token_ucase()
                    else:
                        raise self._nexus_error("Expecting '\"' before beginning SYMBOLS list")
                else:
                    raise self._nexus_error("Expecting '=' after SYMBOLS keyword")
                token = self._nexus_tokenizer.next_token_ucase()
            elif token == 'GAP':
                token = self._nexus_tokenizer.next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.next_token_ucase()
                    self._gap_char = token
                else:
                    raise self._nexus_error("Expecting '=' after GAP keyword")
                token = self._nexus_tokenizer.next_token_ucase()
            elif token == 'INTERLEAVE':
                token = self._nexus_tokenizer.next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.next_token_ucase()
                    if token.startswith("N"):
                        self._interleave = False
                    else:
                        self._interleave = True
                    token = self._nexus_tokenizer.next_token_ucase()
                else:
                    self._interleave = True
            elif token == 'MISSING':
                token = self._nexus_tokenizer.next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.next_token_ucase()
                    self._missing_char = token
                else:
                    raise self._nexus_error("Expecting '=' after MISSING keyword")
                token = self._nexus_tokenizer.next_token_ucase()
            elif token == 'MATCHCHAR':
                token = self._nexus_tokenizer.next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.next_token_ucase()
                    self._match_char = token
                else:
                    raise self._nexus_error("Expecting '=' after MISSING keyword")
                token = self._nexus_tokenizer.next_token_ucase()
            else:
                token = self._nexus_tokenizer.next_token_ucase()

    def _parse_dimensions_statement(self):
        """
        Processes a DIMENSIONS command. Assumes that the file reader is
        positioned right after the "DIMENSIONS" token in a DIMENSIONS command.
        """
        token = self._nexus_tokenizer.next_token_ucase()
        while token != ';':
            if token == 'NTAX':
                token = self._nexus_tokenizer.next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.next_token_ucase()
                    if token.isdigit():
                        self._file_specified_ntax = int(token)
                    else:
                        raise self._nexus_error('Expecting numeric value for NTAX')
                else:
                    raise self._nexus_error("Expecting '=' after NTAX keyword")
            elif token == 'NCHAR':
                token = self._nexus_tokenizer.next_token_ucase()
                if token == '=':
                    token = self._nexus_tokenizer.next_token_ucase()
                    if token.isdigit():
                        self._file_specified_nchar = int(token)
                    else:
                        raise self._nexus_error("Expecting numeric value for NCHAR")
                else:
                    raise self._nexus_error("Expecting '=' after NCHAR keyword")
            token = self._nexus_tokenizer.next_token_ucase()

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
                data_type=self._data_type,
                taxon_namespace=taxon_namespace,
                title=block_title)
        if self._data_type == "continuous":
            self._process_continuous_matrix_data(char_block)
        else:
            self._process_discrete_matrix_data(char_block)

    def _process_continuous_matrix_data(self, char_block):
        if self._interleave:
            raise NotImplementedError("Continuous characters in NEXUS schema not yet supported")
        taxon_namespace = char_block.taxon_namespace
        token = self._nexus_tokenizer.next_token()
        while token != ';' and not self._nexus_tokenizer.is_eof():
            taxon = self._get_taxon(taxon_namespace=taxon_namespace, label=token)
            while len(char_block[taxon]) < self._file_specified_nchar and not self._nexus_tokenizer.is_eof():
                char_group = self._nexus_tokenizer.next_token(ignore_punctuation="-+")
                char_block[taxon].append(dataobject.CharacterDataCell(value=float(char_group)))
            if len(char_block[taxon]) < self._file_specified_nchar:
                raise self._nexus_error("Insufficient characters given for taxon '%s': expecting %d but only found %d ('%s')" \
                    % (taxon.label, self._file_specified_nchar, len(char_block[taxon]), char_block[taxon].symbols_as_string()))
            token = self._nexus_tokenizer.next_token()

    def _process_discrete_matrix_data(self, char_block):
        if self._data_type == "standard":
            self._build_state_alphabet(char_block, self._symbols)
        taxon_namespace = char_block.taxon_namespace
        token = self._nexus_tokenizer.next_token()
        state_alphabet = char_block.default_state_alphabet
        if self._interleave:
            try:
                while token != ";" and not self._nexus_tokenizer.is_eof():
                    taxon = self._get_taxon(taxon_namespace=taxon_namespace, label=token)
                    self._read_character_states(char_block[taxon], state_alphabet)
                    token = self._nexus_tokenizer.next_token()
            except NexusReader.BlockTerminatedException:
                token = self._nexus_tokenizer.next_token()
        else:
            while token != ';' and not self._nexus_tokenizer.is_eof():
                taxon = self._get_taxon(taxon_namespace=taxon_namespace, label=token)
                self._read_character_states(char_block[taxon], state_alphabet)
                if len(char_block[taxon]) < self._file_specified_nchar:
                    raise self._nexus_error("Insufficient characters given for taxon '%s': expecting %d but only found %d ('%s')" \
                        % (taxon.label, self._file_specified_nchar, len(char_block[taxon]), char_block[taxon].symbols_as_string()))
                token = self._nexus_tokenizer.next_token()

    def _get_state_for_multistate_tokens(self,
            state_char_seq,
            multistate_type,
            state_alphabet):
        state = state_alphabet.match_state(state_char_seq)
        if state is not None:
            return state
        member_states = state_alphabet.get_states(symbols=state_char_seq)
        if member_states is None:
            raise self._nexus_error("Unrecognized state encountered: '{}'".format(state_char_seq))
        state = state_alphabet.match_state(symbols=[ms.symbol for ms in member_states])
        if state is not None:
            return state
        sae = state_alphabet.new_state_alphabet_element(
                symbol=None,
                multistate=multistate_type,
                member_states=member_states)
        return sae

    ###########################################################################
    ## TREE / TREE BLOCK PARSERS

    def _parse_tree_statement(self, taxon_namespace=None):
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
        if self.extract_comment_metadata:
            pre_annotations = self._nexus_tokenizer.pull_comment_metadata()
        if token != '=':
            raise self._nexus_error("Expecting '=' in definition of Tree '%s' but found '%s'" % (tree_name, token))
        tree_comments = self._nexus_tokenizer.comments
        tree = nexustokenizer.tree_from_token_stream(stream_tokenizer=self._nexus_tokenizer,
                taxon_namespace=taxon_namespace,
                translate_dict=self._tree_translate_dict,
                encode_splits=self.encode_splits,
                rooting_interpreter=self.rooting_interpreter,
                finish_node_func=self.finish_node_func,
                extract_comment_metadata=self.extract_comment_metadata,
                store_tree_weights=self.store_tree_weights,
                preserve_underscores=self.preserve_underscores,
                suppress_internal_node_taxa=self.suppress_internal_node_taxa,
                edge_len_type=self.edge_len_type,
                case_sensitive_taxon_labels=self.case_sensitive_taxon_labels)
        tree.label = tree_name
        if self.extract_comment_metadata:
            annotations = nexustokenizer.parse_comment_metadata(tree_comments)
            for annote in annotations:
                tree.annotations.add(annote)
            if pre_annotations:
                for annote in pre_annotations:
                    tree.annotations.add(annote)
        if tree_comments is not None and len(tree_comments) > 0:
            tree.comments.extend(tree_comments)
        if self._nexus_tokenizer.current_token != ';':
            self._nexus_tokenizer.skip_to_semicolon()
        return tree

    def _prepopulate_translate_dict(self, taxon_namespace):
        """
        Get default mapping of numbers to taxon labels (to be overwritten by
        a translate dictionary, if found.
        """
        for i, t in enumerate(taxon_namespace):
            self._tree_translate_dict[i+1] = t

    def _parse_translate_statement(self, taxon_namespace):
        """
        Processes a TRANSLATE command. Assumes that the file reader is
        positioned right after the "TRANSLATE" token in a TRANSLATE command.
        """
        token = self._nexus_tokenizer.current_token
        while True:
            translation_token = self._nexus_tokenizer.next_token()
            translation_label = self._nexus_tokenizer.next_token()
            self._tree_translate_dict[translation_token] = taxon_namespace.require_taxon(label=translation_label)
            token = self._nexus_tokenizer.next_token() # ","
            if (not token) or (token == ';'):
                break
            if token != ',':
                raise self._nexus_error("Expecting ',' in TRANSLATE statement after definition for %s = '%s', but found '%s' instead." % (translation_token, translation_label, token))

    def _parse_trees_block(self):
        token = 'TREES'
        if not self.exclude_trees:
            self._nexus_tokenizer.skip_to_semicolon() # move past BEGIN command
            link_title = None
            taxon_namespace = None
            trees_block = None
            block_title = None
            self._tree_translate_dict.clear()
            while not (token == 'END' or token == 'ENDBLOCK') \
                    and not self._nexus_tokenizer.is_eof() \
                    and not token==None:
                token = self._nexus_tokenizer.next_token_ucase()
                if token == 'LINK':
                    link_title = self._parse_link_statement().get("taxa")
                if token == 'TITLE':
                    token = self._nexus_tokenizer.next_token()
                    block_title = token
                if token == 'TRANSLATE':
                    if not taxon_namespace:
                        taxon_namespace = self._get_taxon_namespace(link_title)
                        self._prepopulate_translate_dict(taxon_namespace)
                    self._parse_translate_statement(taxon_namespace)
                if token == 'TREE':
                    if not taxon_namespace:
                        taxon_namespace = self._get_taxon_namespace(link_title)
                        self._prepopulate_translate_dict(taxon_namespace)
                    if not trees_block:
                        trees_block = self.dataset.new_tree_list(taxon_namespace=taxon_namespace, label=block_title)
#                    if not prepared_to_parse_trees:
#                        self._prepare_to_parse_trees(taxon_namespace)
#                        prepared_to_parse_trees = True
                    tree = self._parse_tree_statement(taxon_namespace)
                    trees_block.append(tree, reindex_taxa=False)
            self._nexus_tokenizer.skip_to_semicolon() # move past END command
        else:
            token = self._consume_to_end_of_block(token)

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
                #self.dataset.define_charset(charset_name, positions)

    def _parse_positions(self, adjust_to_zero_based=True, verify=True):
        """
        Parses a character position list. Expects next character read to be the first item in a position list.
        """
        positions = []
        hyphens_as_tokens = self._nexus_tokenizer.hyphens_as_tokens
        self._nexus_tokenizer.hyphens_as_tokens = True
        token = self._nexus_tokenizer.next_token()
        max_positions = self._file_specified_nchar

        if self._nexus_tokenizer.is_eof() or not token:
            raise self._nexus_error('Unexpected end of file or null token')

        while token != ';' and token != ',' and not self._nexus_tokenizer.is_eof():
            if token:
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
                                                raise self._nexus_error('Expecting other tokens after "\\", but no more found.')
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

        self._nexus_tokenizer.hyphens_as_tokens = hyphens_as_tokens
        positions = list(set(positions))
        positions.sort()
        if verify:
            for position in positions:
                if position > max_positions:
                    raise self._nexus_error("Specified position %d, but maximum position is %d" % (position, max_positions))
        if adjust_to_zero_based:
            positions = [position - 1 for position in positions]
        return positions # make unique and return

    def _consume_to_end_of_block(self, token):
        if token:
            token = token.upper()
        while not (token == 'END' or token == 'ENDBLOCK') \
            and not self._nexus_tokenizer.is_eof() \
            and not token==None:
            self._nexus_tokenizer.skip_to_semicolon()
            token = self._nexus_tokenizer.next_token_ucase()
        return token

    def _read_character_states(self,
            character_data_vector,
            state_alphabet,
            ):
        """
        Reads character sequence data substatement until the number of
        character states read is equal to `self._file_specified_nchar` (with
        multi-state characters, such as '(AG)' counting as a single
        state) or, if `self._interleave` is `True`, until an EOL is
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
        while len(character_data_vector) < self._file_specified_nchar:
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
                    token = self.require_next_token()
                    if token == closing_token:
                        break
                    multistate_tokens.append(token)
                c = "".join(subtokens)
                state = self._get_state_for_multistate_tokens(c, multistate_type, state_alphabet)
                character_data_vector.append(state)
                num_states_read += 1
            elif token == "\r" or token == "\n":
                if self._interleave:
                    break
            elif token == ";":
                raise NexusReader.BlockTerminatedException
            else:
                for c in token:
                    if (self._match_char is not None
                            and (c.upper() == self._match_char)):
                        # TODO! handle "."
                        raise NotImplementedError
                    else:
                        try:
                            state = state_alphabet.full_symbol_state_map[c]
                        except KeyError:
                            raise self._nexus_error("Unrecognized (single) state encountered in '{}': '{}' is not defined in {}".format("".join(char_group),
                                    char,
                                    state_alphabet.symbol_state_map.keys()))
                        character_data_vector.append(state)
        if self._interleave:
            self._nexus_tokenizer.set_capture_eol(False)
        return character_data_vector

