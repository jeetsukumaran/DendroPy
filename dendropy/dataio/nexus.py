#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Implementation of NEXUS-format data reader and writer.
"""

from cStringIO import StringIO
import re

from dendropy import dataobject
from dendropy.utility import texttools
from dendropy.utility import iosys
from dendropy.dataio import nexustokenizer
from dendropy.dataio import newick

###############################################################################
## tree_source_iter

def tree_source_iter(stream, **kwargs):
    """
    Iterates over a NEXUS-formatted source of trees given by file-like object
    `stream`

    The following optional keyword arguments are recognized:

        - `taxon_set` specifies the `TaxonSet` object to be attached to the
           trees parsed and manage their taxa. If not specified, then a
           (single) new `TaxonSet` object will be created and for all the
           `Tree` objects.
        - `encode_splits` specifies whether or not split bitmasks will be
           calculated and attached to the edges.
        - `translate_dict` should provide a dictionary mapping taxon numbers (as
           found in the source) to taxon labels (as defined in the source).
        - `rooted` specifies the default rooting interpretation of the tree (see
           `dendropy.dataio.nexustokenizer` for details).
        - `finish_node_func` is a function that will be applied to each node
           after it has been constructed.
        - `edge_len_type` specifies the type of the edge lengths (int or float)
        - `from_index` 0-based index specifying first tree to actually return

    Only trees will be returned, and any and all character data will
    be skipped. The iterator will span over multiple tree blocks,
    but, because our NEXUS data model implementation currently does
    not recognize multiple taxon collection definnitions, taxa in
    those tree blocks will be aggregated into the same `TaxonSet` (a
    new one created, or the one passed to this method via the
    `taxon_set` argument). This behavior is similar to how multiple
    tree blocks are handled by a full NEXUS data file read.
    """
    reader = NexusReader(**kwargs)
    for i, tree in enumerate(reader.tree_source_iter(stream, **kwargs)):
        yield tree

def generalized_tree_source_iter(stream, **kwargs):
    """
    Diagnoses and handles both NEXUS and NEWICK files.
    """
    stream_tokenizer = nexustokenizer.NexusTokenizer(stream)
    token = stream_tokenizer.read_next_token_ucase()
    format = None
    if token == "#NEXUS":
        format = "nexus"
    else:
        if token == "(":
            format = "newick"
    try:
        stream_tokenizer.stream_handle.seek(0)
    except IOError:
        raise TypeError("File format of non-random access source (such as stdin) must be specified in advance.")
    if format == "nexus":
        return tree_source_iter(stream, **kwargs)
    elif format == "newick":
        return newick.tree_source_iter(stream, **kwargs)
    else:
        raise TypeError("Cannot diagnose file format based on first token found: '%s' (looking for '#NEXUS' or '(')")

###############################################################################
## NexusReader

class NexusReader(iosys.DataReader):
    "Encapsulates loading and parsing of a NEXUS format file."

    def __init__(self, **kwargs):
        """
        Recognized keywords in addition to those of `DataReader` are:

            - `default_rooting` : default root for trees read in
            - `finish_node_func` : function to be applied to each node on a
               tree as soon as it has been instantiated
            - `allow_duplicate_taxon_labels` : if True, allow duplicate labels
              on trees [False]
        """
        iosys.DataReader.__init__(self, **kwargs)
        self.stream_tokenizer = nexustokenizer.NexusTokenizer()
        self.default_rooting = kwargs.get("default_rooting", nexustokenizer.RootingInterpretation.UNKNOWN_DEF_ROOTED)
        self.finish_node_func = kwargs.get("finish_node_func", None)
        self.allow_duplicate_taxon_labels = kwargs.get("allow_duplicate_taxon_labels", False)
        self.reset()

    def read(self, stream, **kwargs):
        """
        Instantiates and returns a DataSet object based on the
        NEXUS-formatted contents given in the file-like object `stream`.

        """
        self.reset()
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        self._prepare_to_read_from_stream(stream)
        self._parse_nexus_file()
        self.reset()
        return self.dataset

    def tree_source_iter(self, stream, **kwargs):
        """
        Iterates over a NEXUS-formatted source of trees.

        The following optional keyword arguments are recognized:

            - `encode_splits` specifies whether or not split bitmasks will be
               calculated and attached to the edges.
            - `translate_dict` should provide a dictionary mapping taxon numbers (as
               found in the source) to taxon labels (as defined in the source).
            - `rooted` specifies the default rooting interpretation of the tree (see
               `dendropy.dataio.nexustokenizer` for details).
            - `finish_node_func` is a function that will be applied to each node
               after it has been constructed.
            - `edge_len_type` specifies the type of the edge lengths (int or float)

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
        if "taxon_set" in kwargs:
            self._current_taxon_set = kwargs["taxon_set"]
        self.stream_tokenizer = nexustokenizer.NexusTokenizer(stream)
        token = self.stream_tokenizer.read_next_token_ucase()
        if token != "#NEXUS":
            raise self.data_format_error("Expecting '#NEXUS', but found '%s'" % token)
        while not self.stream_tokenizer.eof:
            token = self.stream_tokenizer.read_next_token_ucase()
            while token != None and token != 'BEGIN' and not self.stream_tokenizer.eof:
                token = self.stream_tokenizer.read_next_token_ucase()
            token = self.stream_tokenizer.read_next_token_ucase()
            if token == 'TAXA':
                self._parse_taxa_block()
            elif token == 'TREES':
                self._prepare_to_parse_trees()
                self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
                while not (token == 'END' or token == 'ENDBLOCK') \
                    and not self.stream_tokenizer.eof \
                    and not token==None:
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token == 'TRANSLATE':
                        self._parse_translate_statement()
                    if token == 'TREE':
                        tree = self._parse_tree_statement()
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
        self.char_block_type = dataobject.StandardCharacterArray
        self.interleave = False
        self.symbols = "012"
        self.gap_char = '-'
        self.missing_char = '?'
        self.match_char = '.'
        self.tree_translate_dict = {}
        self.tax_label_lookup = {}
        self._current_taxon_set = None

    def data_format_error(self, message):
        """
        Returns an exception object parameterized with line and
        column number values.
        """
        return self.stream_tokenizer.data_format_error(message)

    ###########################################################################
    ## HELPERS

    def _prepare_to_read_from_stream(self, file_obj):
        self.stream_tokenizer = nexustokenizer.NexusTokenizer()
        self.stream_tokenizer.stream_handle = file_obj

    def _consume_to_end_of_block(self, token):
        while not (token == 'END' or token == 'ENDBLOCK') \
            and not self.stream_tokenizer.eof \
            and not token==None:
            self.stream_tokenizer.skip_to_semicolon()
            token = self.stream_tokenizer.read_next_token_ucase()
        return token

    ###########################################################################
    ## DATA MANAGEMENT

    def _get_current_taxon_set(self):
        """
        In the future, we may allow for multiple taxon sets, in which
        case this will return the most-recently read taxon set or
        retrieve a particular one based on a LINK or some such statement.
        """
        if self._current_taxon_set is None:
            if self.bound_taxon_set is None:
                self._current_taxon_set = self.dataset.new_taxon_set()
            else:
                self._current_taxon_set = self.bound_taxon_set
                if self._current_taxon_set not in self.dataset.taxon_sets:
                    self.dataset.taxon_sets.add(self._current_taxon_set)
        return self._current_taxon_set

    def _set_current_taxon_set(self, taxon_set):
        self._current_taxon_set = taxon_set
        if self._current_taxon_set not in self.dataset.taxon_sets:
            self.dataset.taxon_sets.add(self._current_taxon_set)

    current_taxon_set = property(_get_current_taxon_set, _set_current_taxon_set)

    def _get_new_taxon_set(self):
        """For now, we restrict NEXUS data to single taxon sets."""
        return self.current_taxon_set

    ###########################################################################
    ## MAIN STREAM PARSE DRIVER

    def _parse_nexus_file(self):
        "Main file parsing driver."
        finish_node_func = self.finish_node_func
        self.reset()
        token = self.stream_tokenizer.read_next_token()
        if token != "#NEXUS":
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
                        while not (token == 'END' or token == 'ENDBLOCK') \
                                and not self.stream_tokenizer.eof \
                                and not token==None:
                            token = self.stream_tokenizer.read_next_token_ucase()
                            if token == 'DIMENSIONS':
                                self._parse_dimensions_statement()
                            if token == 'FORMAT':
                                self._parse_format_statement()
                            if token == 'MATRIX':
                                self._parse_matrix_statement()
                        self.stream_tokenizer.skip_to_semicolon() # move past END command
                    else:
                        token = self._consume_to_end_of_block(token)
                elif token == 'DATA':
                    self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
                    while not (token == 'END' or token == 'ENDBLOCK') \
                            and not self.stream_tokenizer.eof \
                            and not token==None:
                        token = self.stream_tokenizer.read_next_token_ucase()
                        if token == 'DIMENSIONS':
                            self._parse_dimensions_statement()
                        if token == 'FORMAT':
                            self._parse_format_statement()
                        if token == 'MATRIX':
                            self._parse_matrix_statement()
                    self.stream_tokenizer.skip_to_semicolon() # move past END command
                elif token == 'TREES':
                    self._parse_trees_block()
                else:
                    # unknown block
                    token = self._consume_to_end_of_block(token)

        return self.dataset

    ###########################################################################
    ## TAXA BLOCK PARSERS

    def _parse_taxa_block(self):
        token = ''
        self.stream_tokenizer.skip_to_semicolon() # move past BEGIN statement
        while not (token == 'END' or token == 'ENDBLOCK') \
            and not self.stream_tokenizer.eof \
            and not token==None:
            token = self.stream_tokenizer.read_next_token_ucase()
            if token == 'DIMENSIONS':
                self._parse_dimensions_statement()
            if token == 'TAXLABELS':
                self._parse_taxlabels_statement()
        self.stream_tokenizer.skip_to_semicolon() # move past END statement

    def _parse_taxlabels_statement(self):
        """
        Processes a TAXLABELS command. Assumes that the file reader is
        positioned right after the "TAXLABELS" token in a TAXLABELS command.
        """
        token = self.stream_tokenizer.read_next_token()
        while token != ';':
            label = token
            if self.current_taxon_set.has_taxon(label=label):
                pass
            elif len(self.current_taxon_set) >= self.file_specified_ntax:
                raise self.data_format_error("Cannot add '%s':" % label \
                                      + " Declared number of taxa (%d) already defined: %s" % (self.file_specified_ntax,
                                          str([("%s" % t.label) for t in self.current_taxon_set])))
                self.current_taxon_set.new_taxon(label=label)
            token = self.stream_tokenizer.read_next_token()

    ###########################################################################
    ## CHARACTER/DATA BLOCK PARSERS AND SUPPORT

    def _build_state_alphabet(self, char_block, symbols):
        sa = dataobject.StateAlphabet()
        for symbol in symbols:
            sa.append(dataobject.StateAlphabetElement(symbol=symbol))
        if self.missing_char:
            sa.append(dataobject.StateAlphabetElement(symbol=self.missing_char,
                                           multistate=dataobject.StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=sa.get_states(symbols=symbols)))
        if self.gap_char:
            sa.append(dataobject.StateAlphabetElement(symbol=self.gap_char,
                                           multistate=dataobject.StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=sa.get_states(symbols=symbols)))
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
                        self.char_block_type = dataobject.DnaCharacterArray
                    elif token == "RNA":
                        self.char_block_type = dataobject.RnaCharacterArray
                    elif token == "PROTEIN":
                        self.char_block_type = dataobject.ProteinCharacterArray
                    else:
                        # defaults to STANDARD elif token == "STANDARD":
                        self.char_block_type = dataobject.StandardCharacterArray
                        self.symbols = "12"
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

    def _parse_matrix_statement(self):
        """
        Processes a MATRIX command. Assumes that the file reader
        is positioned right after the "MATRIX" token in a MATRIX command,
        and that NTAX and NCHAR have been specified accurately.
        """
        if not self.file_specified_ntax:
            raise self.data_format_error('NTAX must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        elif not self.file_specified_nchar:
            raise self.data_format_error('NCHAR must be defined by DIMENSIONS command to non-zero value before MATRIX command')

        char_block = self.dataset.new_char_array(char_array_type=self.char_block_type, \
            taxon_set=self.current_taxon_set)

        if isinstance(char_block, dataobject.StandardCharacterArray):
            self._build_state_alphabet(char_block, self.symbols)

        symbol_state_map = char_block.default_state_alphabet.symbol_state_map()

        token = self.stream_tokenizer.read_next_token()
        while token != ';' and not self.stream_tokenizer.eof:
            taxon = self.current_taxon_set.require_taxon(label=token)
            if taxon not in char_block:
                if not self.exclude_chars:
                    char_block[taxon] = dataobject.CharacterDataVector(taxon=taxon)
            if self.interleave:
                char_group = StringIO()
                while self.stream_tokenizer.current_file_char != '\n' \
                        and self.stream_tokenizer.current_file_char != '\r' \
                        and not self.stream_tokenizer.eof \
                        and token != ";":
                    token = self.stream_tokenizer.read_next_token(ignore_punctuation="{}()")
                    char_group.write(token)
                char_group = char_group.getvalue()
                self._process_chars(char_group, char_block, symbol_state_map, taxon)
                token = self.stream_tokenizer.read_next_token(ignore_punctuation="{}()")
            else:
                while len(char_block[taxon]) < self.file_specified_nchar and not self.stream_tokenizer.eof:
                    char_group = self.stream_tokenizer.read_next_token(ignore_punctuation="{}()")
                    self._process_chars(char_group, char_block, symbol_state_map, taxon)
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
                    state = symbol_state_map[char]
                except KeyError:
                    raise self.data_format_error("Unrecognized (single) state encountered:'%s' is not defined in %s" % (char, symbol_state_map.keys()))
            else:
                if hasattr(char, "open_tag"):
                    state = self._get_state_for_multistate_char(char, char_block.default_state_alphabet)
            if state is None:
                raise self.data_format_error("Unrecognized state encountered:'%s'" % char)
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
            a = texttools.RichString(ambig[i][1:-1])
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

    def _prepare_to_parse_trees(self):
            self.tree_translate_dict = {}
            self.tax_label_lookup = {}
            for n, t in enumerate(self.current_taxon_set):
                self.tree_translate_dict[str(n + 1)] = t
            # add labels second so that numbers have priority over number
            for n, t in enumerate(self.current_taxon_set):
                l = t.label
                self.tree_translate_dict[l] = t
                self.tax_label_lookup[l] = t
                if self.encode_splits:
                    ti = self.current_taxon_set.index(t)
                    t.clade_mask = (1 << ti)

    def _parse_tree_statement(self):
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

        rooted = self.default_rooting
        if rooted == nexustokenizer.RootingInterpretation.UNKNOWN_DEF_ROOTED \
                or rooted == nexustokenizer.RootingInterpretation.UNKNOWN_DEF_UNROOTED:
            for c in self.stream_tokenizer.comments:
                if c == '&U' or c == '&u':
                    rooted = nexustokenizer.RootingInterpretation.UNROOTED
                    break
                elif c == '&R' or c == '&r':
                    rooted = nexustokenizer.RootingInterpretation.ROOTED
                    break
        tree = nexustokenizer.parse_tree_from_stream(stream_tokenizer=self.stream_tokenizer,
                taxon_set=self.current_taxon_set,
                translate_dict=self.tree_translate_dict,
                encode_splits=self.encode_splits,
                rooted=rooted,
                finish_node_func=self.finish_node_func)
        tree.label = tree_name

        if self.stream_tokenizer.current_token != ';':
            self.stream_tokenizer.skip_to_semicolon()
        return tree

    def _parse_translate_statement(self):
        """
        Processes a TRANSLATE command. Assumes that the file reader is
        positioned right after the "TRANSLATE" token in a TRANSLATE command.
        """
        token = self.stream_tokenizer.current_token
        while True:
            translation_token = self.stream_tokenizer.read_next_token()
            translation_label = self.stream_tokenizer.read_next_token()
            t = self.tax_label_lookup.get(translation_label)
            if t is None:
                t = self.current_taxon_set.require_taxon(label=translation_label)
            self.tree_translate_dict[translation_token] = t

            token = self.stream_tokenizer.read_next_token() # ","
            if (not token) or (token == ';'):
                break
            if token != ',':
                raise self.data_format_error("Expecting ',' in TRANSLATE statement after definition for %s = '%s', but found '%s' instead." % (translation_token, translation_label, token))

    def _parse_trees_block(self):
        token = 'TREES'
        if not self.exclude_trees:
            trees_block = self.dataset.new_tree_list(taxon_set=self.current_taxon_set)
            self._prepare_to_parse_trees()
            self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
            while not (token == 'END' or token == 'ENDBLOCK') \
                and not self.stream_tokenizer.eof \
                and not token==None:
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == 'TRANSLATE':
                    self._parse_translate_statement()
                if token == 'TREE':
                    tree = self._parse_tree_statement()
                    trees_block.append(tree, reindex_taxa=False)
            self.stream_tokenizer.skip_to_semicolon() # move past END command
        else:
            token = self.consume_to_end_of_block(token)

###############################################################################
## NexusWriter

class NexusWriter(iosys.DataWriter):
    "Implements the DataWriter interface for handling NEXML files."

    def __init__(self, **kwargs):
        """
        Recognized keywords in addition to those of `DataReader` are:

            - `simple` : if True, write in simple NEXUS format, i.e. in a
              single "DATA" block, instead of separate "TAXA" and "CHARACTER"
              blocks. [False]
            - `taxa_block` : if False, do not write a "TAXA" block [True]
            - `write_rooting` : if False, do not write a comment before each
              tree indicating its rooting state [True]
            - `edge_lengths` : if False, edges will not write edge lengths [True]
            - `internal_labels` : if False, internal labels will not be written [True]
            - `comment` : list of lines of text to be added as comments to the
              file
        """
        iosys.DataWriter.__init__(self, **kwargs)
        self.simple = kwargs.get("simple", False)
        self.exclude_taxa = kwargs.get("exclude_taxa", True)
        self.is_write_rooting = kwargs.get("write_rooting", True)
        self.is_write_edge_lengths = kwargs.get("edge_lengths", True)
        self.is_write_internal_labels = kwargs.get("internal_labels", True)
        self.spaces_to_underscore = kwargs.get("spaces_to_underscore", False)
        self.comment = kwargs.get("comment", [])

    def write(self, stream, **kwargs):
        """
        Writes bound `DataSource` or `TaxonDomain` to the file-like object
        `stream`.
        """
        assert self.dataset is not None, \
            "NexusWriter instance is not bound to a DataSet: no source of data"
        stream.write('#NEXUS\n\n')
        if self.comment is not None:
            if isinstance(self.comment, list):
                for line in self.comment:
                    if line.strip().replace("\n", "").replace("\r", ""):
                        stream.write("[ %s ]\n" % line)
                    else:
                        stream.write("\n")
                stream.write("\n")
            else:
                stream.write("[ %s ]\n\n" % self.comment)
        if (( (not self.exclude_chars) and self.dataset.char_arrays) \
                or ( (not self.exclude_trees) and self.dataset.tree_lists)) \
                and (not self.simple) \
                and (not self.exclude_taxa):
            for taxon_set in self.dataset.taxon_sets:
                if self.bound_taxon_set is None or taxon_set is self.bound_taxon_set:
                    self.write_taxa_block(taxon_set, stream=stream)
        if not self.exclude_chars:
            for char_array in self.dataset.char_arrays:
                if self.bound_taxon_set is None or char_array.taxon_set is self.bound_taxon_set:
                    self.write_char_block(char_array=char_array, stream=stream)
        if not self.exclude_trees:
            for tree_list in self.dataset.tree_lists:
                if self.bound_taxon_set is None or tree_list.taxon_set is self.bound_taxon_set:
                    self.write_trees_block(tree_list=tree_list, stream=stream)

    def write_taxa_block(self, taxon_set, stream):
        block = []
        block.append('begin taxa;')
        block.append('    dimensions ntax=%d;' % len(taxon_set))
        block.append('    taxlabels')
        for taxon in taxon_set:
            block.append('        %s' % texttools.escape_nexus_token(taxon.label, spaces_to_underscore=self.spaces_to_underscore))
        block.append('  ;')
        block.append('end;\n\n')
        stream.write('\n'.join(block))

    def write_trees_block(self, tree_list, stream):
        block = []
        newick_writer = newick.NewickWriter(edge_lengths=self.is_write_edge_lengths,
            internal_labels=self.is_write_internal_labels)
        block.append('begin trees;')
        for treeidx, tree in enumerate(tree_list):
            if tree.label:
                tree_name = tree.label
            else:
                tree_name = str(treeidx)
            newick_str = newick_writer.compose_node(tree.seed_node)
            if tree.is_rooted and self.is_write_rooting:
                rooting = "[&R] "
            elif not tree.is_rooted and self.is_write_rooting:
                rooting = "[&U] "
            else:
                rooting = ""
            block.append('    tree %s = %s%s;' % (texttools.escape_nexus_token(tree_name, spaces_to_underscore=self.spaces_to_underscore),
                rooting,
                newick_str))
        block.append('end;\n\n')
        stream.write('\n'.join(block))

    def write_char_block(self, char_array, stream, simple_nexus=False):
        nexus = []
        taxlabels = [texttools.escape_nexus_token(taxon.label, spaces_to_underscore=self.spaces_to_underscore) \
                for taxon in char_array.taxon_set]
        max_label_len = max([len(label) for label in taxlabels])
        nchar = max([len(seq) for seq in char_array.values()])
        if simple_nexus:
            nexus.append('begin data;')
            ntaxstr = "ntax=%d" % len(taxlabels)
        else:
            nexus.append('begin characters;')
            ntaxstr = ""
        nexus.append('    dimensions %s nchar=%d;' % (ntaxstr, nchar))
        nexus.append('    format %s;' % self.compose_format_terms(char_array))
        nexus.append('    matrix')
        state_string_map = {}
        for taxon in char_array.taxon_set:
            seq_vec = char_array[taxon]
            seq = StringIO()
            for cell in seq_vec:
                state = cell.value
                assert state is not None, "Undefined state encountered in character sequence."
                try:
                    seq.write(state_string_map[state])
                except:
                    if state.symbol is not None:
                        state_string_map[state] = state.symbol
                    elif state.multistate == dataobject.StateAlphabetElement.AMBIGUOUS_STATE:
                        state_string_map[state] = "{%s}" % ("".join(state.fundamental_symbols))
                    elif state.multistate == dataobject.StateAlphabetElement.POLYMORPHIC_STATE:
                        state_string_map[state] = "(%s)" % ("".join(state.fundamental_symbols))
                    else:
                        raise Exception("Could not match character state to symbol: '%s'." % state)
                    seq.write(state_string_map[state])
            nexus.append('%s    %s' % (texttools.escape_nexus_token(taxon.label, spaces_to_underscore=self.spaces_to_underscore).ljust(max_label_len), seq.getvalue()))
        nexus.append('    ;')
        nexus.append('end;\n\n')
        stream.write('\n'.join(nexus))

    def compose_format_terms(self, char_array):
        format = []
        if isinstance(char_array, dataobject.DnaCharacterArray):
            format.append("datatype=dna")
            format.append("gap=- missing=? matchchar=.")
        elif isinstance(char_array, dataobject.RnaCharacterArray):
            format.append("datatype=rna")
            format.append("gap=- missing=? matchchar=.")
        elif isinstance(char_array, dataobject.ProteinCharacterArray):
            format.append("datatype=protein")
            format.append("gap=- missing=? matchchar=.")
        else:
            format.append("datatype=standard")

            fundamental_symbols = set()
            for state_alphabet in char_array.state_alphabets:
                for s in state_alphabet.fundamental_states():
                    if s.symbol is not None:
                        fundamental_symbols.add(s.symbol)
                    else:
                        raise Exception("Could not match character state to symbol: '%s'." % s)
            format.append('symbols="%s"' % "".join(fundamental_symbols))

            equates = set()
            for state_alphabet in char_array.state_alphabets:
                for a in state_alphabet.ambiguous_states():
                    if a.symbol == "?":
                        format.append("missing=?")
                    elif a.symbol == "-":
                        format.append("gap=-")
                    else:
                        if a.symbol is not None:
                            equates.append("%s={%s}" % (a.symbol, "".join(a.fundamental_symbols())))

            for state_alphabet in char_array.state_alphabets:
                for p in state_alphabet.polymorphic_states():
                    if p.symbol is not None:
                        equates.append("%s=(%s)" % (p.symbol, "".join(p.fundamental_symbols())))

            if equates:
                format.append('equate="%s"' % equates)

        return ' '.join(format)

