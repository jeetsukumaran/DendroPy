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
Pure-Python implementation of NEXUS-schema data reader.
"""

from cStringIO import StringIO
import re

from dendropy import dataobject
from dendropy.utility import textutils
from dendropy.utility import iosys
from dendropy.dataio import nexustokenizer

###############################################################################
## NexusReader

class NexusReader(iosys.DataReader):
    "Encapsulates loading and parsing of a NEXUS schema file."

    def __init__(self, **kwargs):
        """
        __init__ recognizes the following keywords (in addition to
        those of `DataReader.__init__`):

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
        self.reset()
        self.rooting_interpreter = kwargs.get("rooting_interpreter", nexustokenizer.RootingInterpreter(**kwargs))
        self.finish_node_func = kwargs.get("finish_node_func", None)
        self.allow_duplicate_taxon_labels = kwargs.get("allow_duplicate_taxon_labels", False)
        self.preserve_underscores = kwargs.get('preserve_underscores', False)
        self.suppress_internal_node_taxa = kwargs.get("suppress_internal_node_taxa", False)
        self.hyphens_as_tokens = kwargs.get('hyphens_as_tokens', nexustokenizer.DEFAULT_HYPHENS_AS_TOKENS)
        self.store_tree_weights = kwargs.get('store_tree_weights', False)
        self.extract_comment_metadata = kwargs.get('extract_comment_metadata', False)
        self.case_sensitive_taxon_labels = kwargs.get('case_sensitive_taxon_labels', False)
        self.edge_len_type = kwargs.get('edge_len_type', float)

    def read(self, stream):
        """
        Instantiates and returns a DataSet object based on the
        NEXUS-formatted contents given in the file-like object `stream`.
        """
        self.reset()
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        self._prepare_to_read_from_stream(stream)
        self._parse_nexus_file()
        return self.dataset

    def tree_source_iter(self, stream):
        """
        Iterates over a NEXUS-formatted source of trees.
        Only trees will be returned, and any and all character data will
        be skipped. The iterator will span over multiple tree blocks,
        but, because our NEXUS data model implementation currently does
        not recognize multiple taxon collection definnitions, taxa in
        those tree blocks will be aggregated into the same `TaxonSet` (a
        new one created, or the one passed to this method via the
        `taxon_set` argument). This behavior is similar to how multiple
        tree blocks are handled by a full NEXUS data file read.
        """
        self.reset()
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        self.stream_tokenizer = nexustokenizer.NexusTokenizer(stream,
                preserve_underscores=self.preserve_underscores,
                hyphens_as_tokens=self.hyphens_as_tokens,
                extract_comment_metadata=self.extract_comment_metadata)
        token = self.stream_tokenizer.read_next_token_ucase()
        if token.upper() != "#NEXUS":
            raise self.data_format_error("Expecting '#NEXUS', but found '%s'" % token)
        while not self.stream_tokenizer.eof:
            token = self.stream_tokenizer.read_next_token_ucase()
            while token != None and token != 'BEGIN' and not self.stream_tokenizer.eof:
                token = self.stream_tokenizer.read_next_token_ucase()
            token = self.stream_tokenizer.read_next_token_ucase()
            if token == 'TAXA':
                self._parse_taxa_block()
            elif token == 'TREES':
                self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
                link_title = None
                taxon_set = None
                self.tree_translate_dict.clear()
                while not (token == 'END' or token == 'ENDBLOCK') \
                        and not self.stream_tokenizer.eof \
                        and not token==None:
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token == 'LINK':
                        link_title = self._parse_link_statement().get('taxa')
                    if token == 'TRANSLATE':
                        if not taxon_set:
                            taxon_set = self._get_taxon_set(link_title)
                            self._prepopulate_translate_dict(taxon_set)
                        self._parse_translate_statement(taxon_set)
                    if token == 'TREE':
                        if not taxon_set:
                            taxon_set = self._get_taxon_set(link_title)
                            self._prepopulate_translate_dict(taxon_set)
                        tree = self._parse_tree_statement(taxon_set)
                        yield tree
                self.stream_tokenizer.skip_to_semicolon() # move past END command
            else:
                # unknown block
                while not (token == 'END' or token == 'ENDBLOCK') \
                    and not self.stream_tokenizer.eof \
                    and not token==None:
                    self.stream_tokenizer.skip_to_semicolon()
                    token = self.stream_tokenizer.read_next_token_ucase()
        self.reset()

    def reset(self):
        self.char_block_type = dataobject.StandardCharacterMatrix
        self.interleave = False
        self.symbols = ""
        self.gap_char = '-'
        self.missing_char = '?'
        self.match_char = '.'
        self.tree_translate_dict = {}
        self.taxa_blocks = {}
        self.file_specified_ntax = None
        self.file_specified_nchar = None

    def data_format_error(self, message):
        """
        Returns an exception object parameterized with line and
        column number values.
        """
        return self.stream_tokenizer.data_format_error(message)

    def too_many_taxa_error(self, taxon_set, label):
        """
        Returns an exception object parameterized with line and
        column number values.
        """
        return self.stream_tokenizer.too_many_taxa_error(taxon_set=taxon_set,
                max_taxa=self.file_specified_ntax,
                label=label)

    ###########################################################################
    ## HELPERS

    def _prepare_to_read_from_stream(self, file_obj):
        self.stream_tokenizer = nexustokenizer.NexusTokenizer(stream_handle=file_obj,
                preserve_underscores=self.preserve_underscores,
                hyphens_as_tokens=self.hyphens_as_tokens,
                extract_comment_metadata=self.extract_comment_metadata)

    def _consume_to_end_of_block(self, token):
        if token:
            token = token.upper()
        while not (token == 'END' or token == 'ENDBLOCK') \
            and not self.stream_tokenizer.eof \
            and not token==None:
            self.stream_tokenizer.skip_to_semicolon()
            token = self.stream_tokenizer.read_next_token_ucase()
        return token

    ###########################################################################
    ## DATA MANAGEMENT

    def _new_taxon_set(self, title=None):
        if self.attached_taxon_set is not None:
            return self.attached_taxon_set
        if title is None:
            title = 'DEFAULT'
        taxon_set = self.dataset.new_taxon_set(label=title)
        self.taxa_blocks[title] = taxon_set
        return taxon_set

    def _get_char_matrix(self, title=None):
        if title is None:
            if len(self.dataset.char_matrices) == 1:
                return self.dataset.char_matrices[0]
            elif len(self.dataset.char_matrices) == 0:
                raise self.data_format_error("No character matrices defined")
            else:
                raise self.data_format_error("Multiple character matrices defined: require 'LINK' statement")
        else:
            target_matrix = None
            for char_matrix in self.dataset.char_matrices:
                if char_matrix.label == title:
                    if target_matrix is not None:
                        raise self.data_format_error("Multiple character matrices with title '%s'" % title)
                    target_matrix = char_matrix
            if target_matrix is None:
                raise self.data_format_error("Character matrix with title '%s' not defined" % title)
            return target_matrix

    def _get_taxon_set(self, title=None):
        if self.attached_taxon_set is not None:
            return self.attached_taxon_set
        if title is None:
            if len(self.taxa_blocks) == 0:
                self.taxa_blocks['DEFAULT'] = self._new_taxon_set()
                return self.taxa_blocks['DEFAULT']
            elif len(self.taxa_blocks) == 1:
                return self.taxa_blocks.values()[0]
            else:
                raise self.data_format_error("Multiple taxa blocks defined: require 'LINK' statement")
        else:
            if title in self.taxa_blocks:
                return self.taxa_blocks[title]
            else:
                raise self.data_format_error("TaxaBlock with title '%s' not found" % title)

    ###########################################################################
    ## MAIN STREAM PARSE DRIVER

    def _parse_nexus_file(self):
        "Main file parsing driver."
        finish_node_func = self.finish_node_func
        token = self.stream_tokenizer.read_next_token()
        if token.upper() != "#NEXUS":
            raise self.data_format_error("Expecting '#NEXUS', but found '%s'" % token)
        else:
            while not self.stream_tokenizer.eof:
                token = self.stream_tokenizer.read_next_token_ucase()
                while token != None and token != 'BEGIN' and not self.stream_tokenizer.eof:
                    token = self.stream_tokenizer.read_next_token_ucase()
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == 'TAXA':
                    self._parse_taxa_block()
                elif token == 'CHARACTERS':
                    if not self.exclude_chars:
                        self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
                        link_title = None
                        block_title = None
                        while not (token == 'END' or token == 'ENDBLOCK') \
                                and not self.stream_tokenizer.eof \
                                and not token==None:
                            token = self.stream_tokenizer.read_next_token_ucase()
                            if token == 'TITLE':
                                token = self.stream_tokenizer.read_next_token()
                                block_title = token
                            if token == "LINK":
                                link_title = self._parse_link_statement().get('taxa')
                            if token == 'DIMENSIONS':
                                self._parse_dimensions_statement()
                            if token == 'FORMAT':
                                self._parse_format_statement()
                            if token == 'MATRIX':
                                self._parse_matrix_statement(block_title=block_title, link_title=link_title)
                        self.stream_tokenizer.skip_to_semicolon() # move past END command
                    else:
                        token = self._consume_to_end_of_block(token)
                elif token == 'DATA':
                    if not self.exclude_chars:
                        self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
                        block_title = None
                        link_title = None
                        while not (token == 'END' or token == 'ENDBLOCK') \
                                and not self.stream_tokenizer.eof \
                                and not token==None:
                            token = self.stream_tokenizer.read_next_token_ucase()
                            if token == 'TITLE':
                                token = self.stream_tokenizer.read_next_token()
                                block_title = token
                            if token == "LINK":
                                link_title = self._parse_link_statement().get('taxa')
                            if token == 'DIMENSIONS':
                                self._parse_dimensions_statement()
                            if token == 'FORMAT':
                                self._parse_format_statement()
                            if token == 'MATRIX':
                                self._parse_matrix_statement(block_title=block_title, link_title=link_title)
                        self.stream_tokenizer.skip_to_semicolon() # move past END command
                    else:
                        token = self._consume_to_end_of_block(token)
                elif token == 'TREES':
                    self._parse_trees_block()
                elif token in ['SETS', 'ASSUMPTIONS', 'CODONS']:
                    if not self.exclude_chars:
                        self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
                        link_title = None
                        block_title = None
                        while not (token == 'END' or token == 'ENDBLOCK') \
                                and not self.stream_tokenizer.eof \
                                and not token==None:
                            token = self.stream_tokenizer.read_next_token_ucase()
                            if token == 'TITLE':
                                token = self.stream_tokenizer.read_next_token()
                                block_title = token
                            if token == "LINK":
                                link_title = self._parse_link_statement().get('characters')
                            if token == 'CHARSET':
                                self._parse_charset_statement(block_title=block_title, link_title=link_title)
                        self.stream_tokenizer.skip_to_semicolon() # move past END command
                else:
                    # unknown block
                    token = self._consume_to_end_of_block(token)

        return self.dataset

    ###########################################################################
    ## TAXA BLOCK / LINK PARSERS

    def _parse_taxa_block(self):
        token = ''
        self.stream_tokenizer.allow_eof = False
        self.stream_tokenizer.skip_to_semicolon() # move past BEGIN statement
        title = None
        taxon_set = None
        #while not (token == 'END' or token == 'ENDBLOCK') \
        #    and not self.stream_tokenizer.eof \
        #    and not token==None:
        while not (token == 'END' or token == 'ENDBLOCK'):
            token = self.stream_tokenizer.read_next_token_ucase()
            if token == "TITLE":
                token = self.stream_tokenizer.read_next_token()
                taxon_set = self._new_taxon_set(token)
            if token == 'DIMENSIONS':
                self._parse_dimensions_statement()
            if token == 'TAXLABELS':
                if taxon_set is None:
                    taxon_set = self._new_taxon_set()
                self._parse_taxlabels_statement(taxon_set)
        self.stream_tokenizer.skip_to_semicolon() # move past END statement
        self.stream_tokenizer.allow_eof = True

    def _get_taxon(self, taxon_set, label):
        if not self.file_specified_ntax or len(taxon_set) < self.file_specified_ntax:
            taxon = taxon_set.require_taxon(label=label, case_insensitive=not self.case_sensitive_taxon_labels)
        else:
            taxon = taxon_set.get_taxon(label=label, case_insensitive=not self.case_sensitive_taxon_labels)
        if taxon is None:
            raise self.too_many_taxa_error(taxon_set=taxon_set, label=label)
        return taxon

    def _parse_taxlabels_statement(self, taxon_set=None):
        """
        Processes a TAXLABELS command. Assumes that the file reader is
        positioned right after the "TAXLABELS" token in a TAXLABELS command.
        """
        if taxon_set is None:
            taxon_set = self._get_taxon_set()
        token = self.stream_tokenizer.read_next_token()
        while token != ';':
            label = token
            if taxon_set.has_taxon(label=label):
                pass
            elif len(taxon_set) >= self.file_specified_ntax and not self.attached_taxon_set:
                raise self.too_many_taxa_error(taxon_set=taxon_set, label=label)
            else:
                taxon_set.require_taxon(label=label)
            token = self.stream_tokenizer.read_next_token()

    def _parse_link_statement(self):
        """
        Processes a MESQUITE 'LINK' statement.
        """
        # TODO: this is now pretty ugly
        # need to refactor with more abstraction
        links = {}
        token = self.stream_tokenizer.read_next_token_ucase()
        while token != ';':
            if token == 'TAXA':
                token = self.stream_tokenizer.read_next_token()
                if token != "=":
                    raise self.data_format_error("expecting '=' after link taxa")
                token = self.stream_tokenizer.read_next_token()
                links['taxa'] = token
                token = self.stream_tokenizer.read_next_token()
            if token == 'CHARACTERS':
                token = self.stream_tokenizer.read_next_token()
                if token != "=":
                    raise self.data_format_error("expecting '=' after link characters")
                token = self.stream_tokenizer.read_next_token()
                links['characters'] = token
                token = self.stream_tokenizer.read_next_token()
        if token != ";":
            self.stream_tokenizer.skip_to_semicolon()
        return links

    ###########################################################################
    ## CHARACTER/DATA BLOCK PARSERS AND SUPPORT

    def _build_state_alphabet(self, char_block, symbols):
        sa = dataobject.get_state_alphabet_from_symbols(symbols,
                gap_symbol=self.gap_char,
                missing_symbol=self.missing_char)
        char_block.state_alphabets = [sa]
        char_block.default_state_alphabet = char_block.state_alphabets[0]

    def _parse_format_statement(self):
        """
        Processes a FORMAT command. Assumes that the file reader is
        positioned right after the "FORMAT" token in a FORMAT command.
        """
        token = self.stream_tokenizer.read_next_token_ucase()
        while token != ';':
            if token == 'DATATYPE':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token == "DNA" or token == "NUCLEOTIDES":
                        self.char_block_type = dataobject.DnaCharacterMatrix
                    elif token == "RNA":
                        self.char_block_type = dataobject.RnaCharacterMatrix
                    elif token == "NUCLEOTIDE":
                        self.char_block_type = dataobject.NucleotideCharacterMatrix
                    elif token == "PROTEIN":
                        self.char_block_type = dataobject.ProteinCharacterMatrix
                    elif token == "CONTINUOUS":
                        self.char_block_type = dataobject.ContinuousCharacterMatrix
                    else:
                        # defaults to STANDARD elif token == "STANDARD":
                        self.char_block_type = dataobject.StandardCharacterMatrix
                        self.symbols = "01"
                else:
                    raise self.data_format_error("Expecting '=' after DATATYPE keyword")
                token = self.stream_tokenizer.read_next_token_ucase()
            elif token == 'SYMBOLS':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token == '"':
                        self.symbols = ""
                        token = self.stream_tokenizer.read_next_token_ucase()
                        while token != '"':
                            if token not in self.symbols:
                                self.symbols = self.symbols + token
                            token = self.stream_tokenizer.read_next_token_ucase()
                    else:
                        raise self.data_format_error("Expecting '\"' before beginning SYMBOLS list")
                else:
                    raise self.data_format_error("Expecting '=' after SYMBOLS keyword")
                token = self.stream_tokenizer.read_next_token_ucase()
            elif token == 'GAP':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    self.gap_char = token
                else:
                    raise self.data_format_error("Expecting '=' after GAP keyword")
                token = self.stream_tokenizer.read_next_token_ucase()
            elif token == 'INTERLEAVE':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token.startswith("N"):
                        self.interleave = False
                    else:
                        self.interleave = True
                    token = self.stream_tokenizer.read_next_token_ucase()
                else:
                    self.interleave = True
            elif token == 'MISSING':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    self.missing_char = token
                else:
                    raise self.data_format_error("Expecting '=' after MISSING keyword")
                token = self.stream_tokenizer.read_next_token_ucase()
            elif token == 'MATCHCHAR':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    self.match_char = token
                else:
                    raise self.data_format_error("Expecting '=' after MISSING keyword")
                token = self.stream_tokenizer.read_next_token_ucase()
            else:
                token = self.stream_tokenizer.read_next_token_ucase()

    def _parse_dimensions_statement(self):
        """
        Processes a DIMENSIONS command. Assumes that the file reader is
        positioned right after the "DIMENSIONS" token in a DIMENSIONS command.
        """
        token = self.stream_tokenizer.read_next_token_ucase()
        while token != ';':
            if token == 'NTAX':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token.isdigit():
                        self.file_specified_ntax = int(token)
                    else:
                        raise self.data_format_error('Expecting numeric value for NTAX')
                else:
                    raise self.data_format_error("Expecting '=' after NTAX keyword")
            elif token == 'NCHAR':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token.isdigit():
                        self.file_specified_nchar = int(token)
                    else:
                        raise self.data_format_error("Expecting numeric value for NCHAR")
                else:
                    raise self.data_format_error("Expecting '=' after NCHAR keyword")
            token = self.stream_tokenizer.read_next_token_ucase()

    def _parse_matrix_statement(self, block_title=None, link_title=None):
        """
        Processes a MATRIX command. Assumes that the file reader
        is positioned right after the "MATRIX" token in a MATRIX command,
        and that NTAX and NCHAR have been specified accurately.
        """
        if not self.file_specified_ntax:
            raise self.data_format_error('NTAX must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        elif not self.file_specified_nchar:
            raise self.data_format_error('NCHAR must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        taxon_set = self._get_taxon_set(link_title)
        char_block = self.dataset.new_char_matrix(
            char_matrix_type=self.char_block_type, \
            taxon_set=taxon_set,
            label=block_title)
        if isinstance(char_block, dataobject.ContinuousCharacterMatrix):
            self._process_continuous_matrix_data(char_block)
        else:
            self._process_discrete_matrix_data(char_block)

    def _process_continuous_matrix_data(self, char_block):
        taxon_set = char_block.taxon_set
        token = self.stream_tokenizer.read_next_token()
        while token != ';' and not self.stream_tokenizer.eof:
            taxon = self._get_taxon(taxon_set=taxon_set, label=token)
            if taxon not in char_block:
                char_block[taxon] = dataobject.CharacterDataVector(taxon=taxon)
                if self.interleave:
                    raise NotImplementedError("Continuous characters in NEXUS schema not yet supported")
                else:
                    while len(char_block[taxon]) < self.file_specified_nchar and not self.stream_tokenizer.eof:
                        char_group = self.stream_tokenizer.read_next_token(ignore_punctuation="-+")
                        char_block[taxon].append(dataobject.CharacterDataCell(value=float(char_group)))
                    if len(char_block[taxon]) < self.file_specified_nchar:
                        raise self.data_format_error("Insufficient characters given for taxon '%s': expecting %d but only found %d ('%s')" \
                            % (taxon.label, self.file_specified_nchar, len(char_block[taxon]), char_block[taxon].symbols_as_string()))
                    token = self.stream_tokenizer.read_next_token()

    def _process_discrete_matrix_data(self, char_block):
        if isinstance(char_block, dataobject.StandardCharacterMatrix):
            self._build_state_alphabet(char_block, self.symbols)
        taxon_set = char_block.taxon_set
        symbol_state_map = char_block.default_state_alphabet.symbol_state_map()
        token = self.stream_tokenizer.read_next_token()

        if self.interleave:
            try:
                while token != ";" and not self.stream_tokenizer.eof:
                    taxon = self._get_taxon(taxon_set=taxon_set, label=token)
                    if taxon not in char_block:
                        char_block[taxon] = dataobject.CharacterDataVector(taxon=taxon)
                    tokens = self.stream_tokenizer.read_statement_tokens_till_eol(ignore_punctuation="{}()")
                    if tokens is not None:
                        self._process_chars(''.join(tokens), char_block, symbol_state_map, taxon)
                    else:
                        break
                    token = self.stream_tokenizer.read_next_token()
            except nexustokenizer.NexusTokenizer.BlockTerminatedException:
                if tokens is not None:
                    self._process_chars(''.join(tokens), char_block, symbol_state_map, taxon)
                token = self.stream_tokenizer.read_next_token()
        else:
            while token != ';' and not self.stream_tokenizer.eof:
                taxon = self._get_taxon(taxon_set=taxon_set, label=token)
                if taxon not in char_block:
                    char_block[taxon] = dataobject.CharacterDataVector(taxon=taxon)
                while len(char_block[taxon]) < self.file_specified_nchar \
                        and not self.stream_tokenizer.eof:
                    char_group = self.stream_tokenizer.read_matrix_characters(self.file_specified_nchar)
                    if char_group is not None:
                        char_group = "".join(char_group)
                        self._process_chars(char_group, char_block, symbol_state_map, taxon)
                    else:
                        break
                if len(char_block[taxon]) < self.file_specified_nchar:
                    raise self.data_format_error("Insufficient characters given for taxon '%s': expecting %d but only found %d ('%s')" \
                        % (taxon.label, self.file_specified_nchar, len(char_block[taxon]), char_block[taxon].symbols_as_string()))
                token = self.stream_tokenizer.read_next_token()

    def _process_chars(self, char_group, char_block, symbol_state_map, taxon):
        if self.exclude_chars:
            return
        if not char_group:
            return
        char_group = self._parse_nexus_multistate(char_group)
        for char in char_group:
            if len(char) == 1:
                try:
                    state = symbol_state_map[char.upper()]
                except KeyError:
                    if self.match_char is not None \
                        and char.upper() == self.match_char.upper():
                            state = char_block[0][len(char_block[taxon])].value
                    else:
                        raise self.data_format_error("Unrecognized (single) state encountered in '%s': '%s' is not defined in %s" % ("".join(char_group), char, symbol_state_map.keys()))
            else:
                if hasattr(char, "open_tag"):
                    state = self._get_state_for_multistate_char(char, char_block.default_state_alphabet)
                else:
                    raise self.data_format_error("Multiple character state without multi-state mark-up: '%s'" % char)
            if state is None:
                raise self.data_format_error("Unrecognized state encountered: '%s'" % char)
            char_block[taxon].append(dataobject.CharacterDataCell(value=state))

    def _parse_nexus_multistate(self, seq):
        """
        Given a sequence of characters, with ambiguities denoted by
        `{<STATES>}`, this returns a list of characters, with unambiguous
        characters as individual elements, and the ambiguous characters in their
        own string elements. E.g.:

            "ACTG(AC)GGT(CGG)(CG)GG"

        results in:

            ['A', 'C', 'T', 'G', 'AC', 'G', 'G', 'T', 'CGG', 'CG', 'G', 'G']

        Two attributes are also added to every set of ambiguous characters,
        `open_tag` and `close_tag` with their values set to the opening and closing
        tokens.
        """
        spat = re.compile('[\(|\{].+?[\)\}]')
        mpat = re.compile('([\(|\{].+?[\)\}])')

        unambig = spat.split(seq)
        ambig = mpat.findall(seq)
        result = []
        for i in xrange(len(unambig)-1):
            a = textutils.RichString(ambig[i][1:-1])
            a.open_tag = ambig[i][0]
            a.close_tag = ambig[i][-1]
            result.extend(unambig[i])
            result.append(a)
        result.extend(unambig[-1])
        return [c for c in result if c]

    def _get_state_for_multistate_char(self, char, state_alphabet):
        state = state_alphabet.match_state(char)
        if state is not None:
            return state
        if hasattr(char, "open_tag") and char.open_tag == '{':
            multistate_type = dataobject.StateAlphabetElement.AMBIGUOUS_STATE
        elif hasattr(char, "open_tag") and char.open_tag == '(':
            multistate_type = dataobject.StateAlphabetElement.POLYMORPHIC_STATE
        else:
            return None
        member_states = state_alphabet.get_states(symbols=char)
        if member_states is None:
            return None
        state = state_alphabet.match_state(symbols=[ms.symbol for ms in member_states])
        if state is not None:
            return state
        sae = dataobject.StateAlphabetElement(symbol=None,
            multistate=multistate_type,
            member_states=member_states)
        state_alphabet.append(sae)
        return sae

    ###########################################################################
    ## TREE / TREE BLOCK PARSERS

#    def _prepare_to_parse_trees(self, taxon_set):
#        self.tree_translate_dict = {}
#        self.tax_label_lookup = {}
#        for n, t in enumerate(taxon_set):
#            self.tree_translate_dict[str(n + 1)] = t
#        # add labels second so that labels have priority over number
#        for n, t in enumerate(taxon_set):
#            l = t.label
#            self.tree_translate_dict[l] = t
#            self.tax_label_lookup[l] = t
#            if self.encode_splits:
#                ti = taxon_set.index(t)
#                t.split_bitmask = (1 << ti)

    def _parse_tree_statement(self, taxon_set=None):
        """
        Processes a TREE command. Assumes that the file reader is
        positioned right after the "TREE" token in a TREE command.
        Calls on the NewickStatementParser of the trees module.
        """
        token = self.stream_tokenizer.read_next_token()
        if token == '*':
            token = self.stream_tokenizer.read_next_token()
        tree_name = token
        token = self.stream_tokenizer.read_next_token()
        if token != '=':
            raise self.data_format_error("Expecting '=' in definition of Tree '%s' but found '%s'" % (tree_name, token))
        tree_comments = self.stream_tokenizer.comments
        tree = nexustokenizer.tree_from_token_stream(stream_tokenizer=self.stream_tokenizer,
                taxon_set=taxon_set,
                translate_dict=self.tree_translate_dict,
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
        if tree_comments is not None and len(tree_comments) > 0:
            tree.comments.extend(tree_comments)
        if self.stream_tokenizer.current_token != ';':
            self.stream_tokenizer.skip_to_semicolon()
        return tree

    def _prepopulate_translate_dict(self, taxon_set):
        """
        Get default mapping of numbers to taxon labels (to be overwritten by
        a translate dictionary, if found.
        """
        for i, t in enumerate(taxon_set):
            self.tree_translate_dict[i+1] = t

    def _parse_translate_statement(self, taxon_set):
        """
        Processes a TRANSLATE command. Assumes that the file reader is
        positioned right after the "TRANSLATE" token in a TRANSLATE command.
        """
        token = self.stream_tokenizer.current_token
        while True:
            translation_token = self.stream_tokenizer.read_next_token()
            translation_label = self.stream_tokenizer.read_next_token()
            self.tree_translate_dict[translation_token] = taxon_set.require_taxon(label=translation_label)
            token = self.stream_tokenizer.read_next_token() # ","
            if (not token) or (token == ';'):
                break
            if token != ',':
                raise self.data_format_error("Expecting ',' in TRANSLATE statement after definition for %s = '%s', but found '%s' instead." % (translation_token, translation_label, token))

    def _parse_trees_block(self):
        token = 'TREES'
        if not self.exclude_trees:
            self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
            link_title = None
            taxon_set = None
            trees_block = None
            block_title = None
            self.tree_translate_dict.clear()
            while not (token == 'END' or token == 'ENDBLOCK') \
                    and not self.stream_tokenizer.eof \
                    and not token==None:
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == 'LINK':
                    link_title = self._parse_link_statement().get("taxa")
                if token == 'TITLE':
                    token = self.stream_tokenizer.read_next_token()
                    block_title = token
                if token == 'TRANSLATE':
                    if not taxon_set:
                        taxon_set = self._get_taxon_set(link_title)
                        self._prepopulate_translate_dict(taxon_set)
                    self._parse_translate_statement(taxon_set)
                if token == 'TREE':
                    if not taxon_set:
                        taxon_set = self._get_taxon_set(link_title)
                        self._prepopulate_translate_dict(taxon_set)
                    if not trees_block:
                        trees_block = self.dataset.new_tree_list(taxon_set=taxon_set, label=block_title)
#                    if not prepared_to_parse_trees:
#                        self._prepare_to_parse_trees(taxon_set)
#                        prepared_to_parse_trees = True
                    tree = self._parse_tree_statement(taxon_set)
                    trees_block.append(tree, reindex_taxa=False)
            self.stream_tokenizer.skip_to_semicolon() # move past END command
        else:
            token = self._consume_to_end_of_block(token)

    def _parse_charset_statement(self, block_title=None, link_title=None):
        """
        Parses a character set description. Assumes token stream is positioned right after 'charset' command.
        """
        char_matrix = self._get_char_matrix(title=link_title)
        keyword = self.stream_tokenizer.current_token
        token = self.stream_tokenizer.read_next_token()
        if self.stream_tokenizer.eof or not token:
            raise self.data_format_error('Unexpected end of file or null token')
        else:
            if not token:
                raise self.data_format_error("Unexpected end of file or null token")
            else:
                charset_name = token
                token = self.stream_tokenizer.read_next_token()
                if not token:
                    raise self.data_format_error("Unexpected end of file or null token")
                elif token != '=':
                    raise self.data_format_error('Expecting "=" after character set name "%s", but instead found "%s"' % (charset_name, token))
                else:
                    positions = self._parse_positions(adjust_to_zero_based=True)
                char_matrix.new_character_subset(charset_name, positions)
                #self.dataset.define_charset(charset_name, positions)

    def _parse_positions(self, adjust_to_zero_based=True, verify=True):
        """
        Parses a character position list. Expects next character read to be the first item in a position list.
        """
        positions = []
        hyphens_as_tokens = self.stream_tokenizer.hyphens_as_tokens
        self.stream_tokenizer.hyphens_as_tokens = True
        token = self.stream_tokenizer.read_next_token()
        max_positions = self.file_specified_nchar

        if self.stream_tokenizer.eof or not token:
            raise self.data_format_error('Unexpected end of file or null token')

        while token != ';' and token != ',' and not self.stream_tokenizer.eof:
            if token:
                if token.upper() == 'ALL':
                    positions = range(1, max_positions + 1)
                    break
                elif token.isdigit():
                    start = int(token)
                    token = self.stream_tokenizer.read_next_token()
                    if token:
                        if token == ',' or token.isdigit() or token == ';':
                            positions.append(start)
                        elif token == '-':
                            token = self.stream_tokenizer.read_next_token()
                            if token:
                                if token.isdigit() or token == '.':
                                    if token == '.':
                                        end = max_positions
                                        #token = self.stream_tokenizer.read_next_token()
                                    else:
                                        end = int(token)
                                        #token = self.stream_tokenizer.read_next_token()
                                    token = self.stream_tokenizer.read_next_token()
                                    if token:
                                        if token == '\\' or token == '/': # (NEXUS standard only accepts '\')
                                            token = self.stream_tokenizer.read_next_token()
                                            if token:
                                                if token.isdigit():
                                                    step = int(token)
                                                    #token = self.stream_tokenizer.read_next_token()
                                                else:
                                                    raise self.data_format_error('Expecting digit but found "%s".' % (token))
                                            else:
                                                raise self.data_format_error('Expecting other tokens after "\\", but no more found.')
                                            token = self.stream_tokenizer.read_next_token()
                                        else:
                                            step = 1
                                    else:
                                        step = 1
                                    for q in range(start, end+1, step):
                                        if q <= max_positions:
                                            positions.append(q)
                                else:
                                    raise self.data_format_error('Expecting digit or ".", but found "%s".' % (token))
                            else:
                                raise self.data_format_error('Expecting other tokens after "-", but no more found.')
                        else:
                            raise self.data_format_error('Expecting digit or "all", but found "%s".' % (token))
                    else:
                        positions.append(start)

        self.stream_tokenizer.hyphens_as_tokens = hyphens_as_tokens
        positions = list(set(positions))
        positions.sort()
        if verify:
            for position in positions:
                if position > max_positions:
                    raise self.data_format_error("Specified position %d, but maximum position is %d" % (position, max_positions))
        if adjust_to_zero_based:
            positions = [position - 1 for position in positions]
        return positions # make unique and return


