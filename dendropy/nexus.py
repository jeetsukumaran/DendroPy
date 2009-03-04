#! /usr/bin/env python

############################################################################
##  nexus.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2007 Jeet Sukumaran and Mark T. Holder.
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
############################################################################

"""
This module wraps routines needed for reading and writing trees (and data) in
NEXUS format. Note that the current parser supports a SINGLE taxa block / 
character block / data block / trees block model. Behavior with multiple taxa 
blocks, data blocks, trees blocks etc. is undefined.
"""

import re
import os
import sys
from cStringIO import StringIO

from dendropy import datasets
from dendropy import taxa
from dendropy import trees
from dendropy import characters   
from dendropy import get_logger
_LOG = get_logger("dendropy.nexus")


class RootingInterpretation:
    UNKNOWN_DEF_ROOTED = 0
    UNKNOWN_DEF_UNROOTED = 1
    ROOTED = 2
    UNROOTED = 3


def build_state_alphabet(char_block, symbols, missing):
    sa = characters.StateAlphabet()
    for symbol in symbols:
        sa.append(characters.StateAlphabetElement(symbol=symbol))
    if self.missing_char:
        sa.append(characters.StateAlphabetElement(symbol=missing,
                                       multistate=characters.StateAlphabetElement.AMBIGUOUS_STATE,
                                       member_states=sa.get_states(symbols=symbols)))
    char_block.state_alphabets = [sa]
    char_block.default_state_alphabet = char_block.state_alphabets[0]
    return char_block

#######################################################################
## NESTED CLASSES

class SyntaxException(Exception):

    def __init__(self, row=None, column=None, message=None):
        self.row = row
        self.column = column
        self.message = message

    def __str__(self):
        if self.row is None:
            t = ""
        else:
            t =  " IN LINE %d" % self.row
        return 'ERROR PARSING FILE%s: %s' % (t, self.message)

############################################################################
## Standard Tree Iterator
      
def iterate_over_trees(file_obj=None, taxa_block=None, dataset=None, **kwargs):
    """
    Generator to iterate over trees in data file.
    Primary goal is to be memory efficient, storing no more than one tree
    at a time. Speed might have to be sacrificed for this!
    """
    if dataset is None:
        dataset = datasets.Dataset()
    stream_tokenizer = PurePythonNexusStreamTokenizer(file_obj)
    token = stream_tokenizer.read_next_token_ucase()
    if token == "#NEXUS":
        file_format = "NEXUS"
    else:
        stream_tokenizer.stream_handle.seek(0)
        file_format = "NEWICK"
    for tree in dataset.iterate_over_trees(file_obj, taxa_block=taxa_block, format=file_format, **kwargs):
        yield tree

############################################################################
## Universal Nex-ish Readers

def read_dataset(file_obj, dataset=None):
    """
    Note: due to usage of seek(), does not work on stream sources 
    (such as stdin). In this case, client must either buffer stream to string
    (e.g. file_obj=StringIO(sys.stdin.read())) or use direct methods.
    """
    stream_tokenizer = PurePythonNexusStreamTokenizer(file_obj)
    token = stream_tokenizer.read_next_token_ucase()
    stream_tokenizer.stream_handle.seek(0)
    if token == "#NEXUS":
        reader = NexusReader()
    else: # assume NEWICK
        reader = NewickReader()
    return reader.read_dataset(file_obj=file_obj, dataset=dataset)

def read_trees(file_obj, dataset=None):
    dataset = read_dataset(file_obj=file_obj, dataset=dataset)
    return dataset.trees_blocks

############################################################################
## support functions

def split_to_newick(split, taxa_block):
    """
    Represents a split as a newick string.
    """
    taxlabels = taxa_block.labels()
    
    # do not do the root
    if split == 0 or (split == taxa_block.all_taxa_bitmask()):
        return "(%s)" % (",".join(taxlabels))
    
    idx = 0
    left = []
    right = []
    while split >= 0 and idx < len(taxlabels):
        if split & 1:
            left.append(taxlabels[idx])
        else:
            right.append(taxlabels[idx])
        idx += 1
        split = split >> 1
    assert ( len(left) + len(right) ) == len(taxlabels)
    return "((%s), (%s))" % (", ".join(left), ", ".join(right))

def map_to_iupac_ambiguity_code(states):
    """
    Given a sequence of characters, maps ambiguities given the form of
    `AG` to IUPAC codes (e.g., `AC -> R').
    """
    if len(states) == 1:
        return states[0]
    states = [state.upper() for state in states]
    if states.count('A') and states.count('C') and states.count('G') and (states.count('T') or states.count('U')):
        return 'N'
    if states.count('A') and states.count('C') and states.count('G'):
        return 'V'
    if states.count('A') and states.count('C') and (states.count('T') or states.count('U')):
        return 'H'
    if states.count('A') and states.count('G') and (states.count('T') or states.count('U')):
        return 'D'
    if states.count('C') and states.count('G') and (states.count('T') or states.count('U')):
        return 'B'
    if states.count('A') and states.count('C'):
        return 'M'
    if states.count('A') and states.count('G'):
        return 'R'
    if states.count('C') and (states.count('T') or states.count('U')):
        return 'W'
    if states.count('C') and states.count('G'):
        return 'S'
    if states.count('C') and (states.count('T') or states.count('U')):
        return 'Y'
    if states.count('G') and (states.count('T') or states.count('U')):
        return 'K'
    raise Exception('Unrecognized characters in "%s"' % states)

def parse_sequence_iupac_ambiguities(seq):
    """
    Given a sequence of characters, with ambiguities denoted by
    `{<STATES>}`, this returns a sequence of characters with the
    ambiguities mapped to the IUPAC codes.
    """
    if isinstance(seq, list):
      as_list = True
      seq = ''.join(seq)
    else:
      as_list = False
    result = ""
    pattern = re.compile('{(.*?)}')
    pos = 0
    match = pattern.search(seq, pos)
    if match:
        while match:
          if pos > 0:
            result = result + seq[pos-1:match.start()]
          else:
            result = result + seq[pos:match.start()]
          result = result + map_to_iupac_ambiguity_code(match.group(1))
          pos = match.end() + 1
          if pos < len(seq) - 1:
            match = pattern.search(seq, pos)
          else:
            match = None
        result = result + seq[pos-1:]
    else:
      return seq
    if as_list:
      return [char for char in result]
    else:
      return result
      
############################################################################
##  Fundamental NEWICK Tree parsers

def parse_newick_string(tree_statement, 
                        taxa_block=None, 
                        translate_dict=None, 
                        encode_splits=False, 
                        rooted=RootingInterpretation.UNKNOWN_DEF_UNROOTED,
                        finish_node_func=None):
    "Processes a (SINGLE) TREE statement string."
    stream_handle = StringIO(tree_statement)
    stream_tokenizer = PurePythonNexusStreamTokenizer(stream_handle)
    tree = parse_newick_tree_stream(stream_tokenizer=stream_tokenizer, 
                                     taxa_block=taxa_block,
                                     translate_dict=translate_dict, 
                                     encode_splits=encode_splits,
                                     rooted=rooted,
                                     finish_node_func=finish_node_func)
    return tree    

class StrToTaxon(object):
    def __init__(self, taxa, translate_dict=None):
        self.taxa = taxa
        self.translate = translate_dict or {}
    def get_taxon(self, label, taxon_required=True):
        v = self.translate.get(label)
        if v is not None:
            return v
        #_LOG.debug("translate_dict miss: %s " % label)
        t = self.taxa.get_taxon(label=label, taxon_required=taxon_required)
        if t is not None:
            self.translate[label] = t #@this could lead to problems when we support multiple taxon blocks, but now it'll speed thing up
        return t
    def index(self, t):
        return self.taxa.index(t)

def parse_newick_tree_stream(stream_tokenizer, 
                             taxa_block=None,
                             translate_dict=None,
                             encode_splits=False,
                             rooted=RootingInterpretation.UNKNOWN_DEF_UNROOTED,
                             finish_node_func=None):
    """
    Processes a (SINGLE) TREE statement. Assumes that the input stream is
    located at the beginning of the statement (i.e., the first
    parenthesis that defines the tree).
    """
    token = stream_tokenizer.read_next_token()
    if not token:
        return None    
    tree = trees.Tree()
    if taxa_block is None:
        taxa_block = taxa.TaxaBlock()    
    tree.taxa_block = taxa_block

    if encode_splits:
        if rooted == RootingInterpretation.UNKNOWN_DEF_ROOTED or rooted == RootingInterpretation.UNKNOWN_DEF_UNROOTED:
            r = stream_tokenizer.tree_rooted_comment()
            if r is not None:
                if r:
                    rooted = RootingInterpretation.ROOTED
                else:
                    rooted = RootingInterpretation.UNROOTED


        if rooted == RootingInterpretation.UNKNOWN_DEF_ROOTED or RootingInterpretation.ROOTED:
            tree.split_edges = {}
        else:
            d = utils.NormalizedBitmaskDict(mask=taxa_block.all_taxa_bitmask())
            tree.split_edges = d
        split_map = tree.split_edges

    stt = StrToTaxon(taxa_block, translate_dict)

    tree.seed_node = trees.Node()
    curr_node = tree.seed_node
    if encode_splits:
        curr_node.edge.clade_mask = 0L

    while True:
        if not token or token == ';':
            if curr_node is not tree.seed_node:
                raise stream_tokenizer.syntax_exception('Unbalanced parentheses -- not enough ")" characters found in tree description')
            if encode_splits:
                split_map[curr_node.edge.clade_mask] = curr_node.edge
            break
        if token == '(':
            tmp_node = trees.Node()
            if encode_splits:
                tmp_node.edge.clade_mask = 0L
            curr_node.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
        elif token == ',':
            tmp_node = trees.Node()
            if curr_node.is_leaf() and not curr_node.taxon:
                raise stream_tokenizer.syntax_exception('Missing taxon specifier in a tree -- found either a "(," or ",," construct.')
            p = curr_node.parent_node
            if not p:
                raise stream_tokenizer.syntax_exception('Comma found one the "outside" of a newick tree description')
                if encode_splits:
                    tmp_node.edge.clade_mask = 0L
                    p.edge.clade_mask |= curr_node.edge.clade_mask
            if finish_node_func is not None:
                finish_node_func(curr_node, tree)
            p.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
        else:
            if token == ')':
                if curr_node.is_leaf() and not curr_node.taxon:
                    raise stream_tokenizer.syntax_exception('Missing taxon specifier in a tree -- found either a "(," or ",," construct.')
                p = curr_node.parent_node
                if not p:
                    raise stream_tokenizer.syntax_exception('Unbalanced parentheses -- too many ")" characters found in tree description')
                if encode_splits:
                    cm |= curr_node.edge.clade_mask
                    p.edge.clade_mask = cm
                    split_map[cm] = curr_node.edge
                if finish_node_func is not None:
                    finish_node_func(curr_node, tree)
                curr_node = p
            else:
                is_leaf = curr_node.is_leaf()
                t = stt.get_taxon(label=token, taxon_required=is_leaf)
                if t is None:
                    curr_node.label = token
                else:
                    curr_node.taxon = t
                    if encode_splits:
                        try:
                            cm = t.clade_mask
                        except:
                            cm = 1 << (stt.index(t))
                        e = curr_node.edge
                        e.clade_mask = cm
                        split_map[cm] = e

            token = stream_tokenizer.read_next_token()    
            if token == ':':
                edge_length_str = stream_tokenizer.read_next_token(ignore_punctuation='-')
                try:
                    curr_node.edge.length = float(edge_length_str)
                except ValueError:
                    curr_node.edge.length = edge_length_str
                token = stream_tokenizer.read_next_token()
    return tree
     
############################################################################
## CLASS: PurePythonNexusStreamTokenizer

class PurePythonNexusStreamTokenizer(object):
    "Encapsulates reading NEXUS/NEWICK tokens from file."

    #######################################################################
    ## STATIC METHODS

    punctuation = '\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>'
    whitespace = ' \0\t\n\r'
    multi_char_needs_quoting = re.compile(' \0\t\n\r\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>_')
    single_char_needs_quoting = re.compile(' \0\t\n\r\[\'_')
    def is_whitespace(char):
        return char in PurePythonNexusStreamTokenizer.whitespace

    is_whitespace = staticmethod(is_whitespace)

    def has_whitespace(s):
        return re.search('['+PurePythonNexusStreamTokenizer.whitespace+']', s) != None

    has_whitespace = staticmethod(has_whitespace)

    def is_punctuation(char):
        return char in PurePythonNexusStreamTokenizer.punctuation

    is_punctuation = staticmethod(is_punctuation)

    def has_punctuation(s):
        return re.search('['+PurePythonNexusStreamTokenizer.punctuation+']', s) != None

    has_punctuation = staticmethod(has_punctuation)

    def is_whitespace_or_punctuation(char):
        return PurePythonNexusStreamTokenizer.is_whitespace(char) or PurePythonNexusStreamTokenizer.is_punctuation(char)

    is_whitespace_or_punctuation = staticmethod(is_whitespace_or_punctuation)

    def has_whitespace_or_punctuation(s):
        return PurePythonNexusStreamTokenizer.has_whitespace(s) or PurePythonNexusStreamTokenizer.has_punctuation(s)

    has_whitespace_or_punctuation = staticmethod(has_whitespace_or_punctuation)

    def validate_identifier(label):
        if PurePythonNexusStreamTokenizer.has_whitespace_or_punctuation(label):
            if (label[0] == "'" and label[1] == "'") or label[0] != "'":
                label = "'" + label
            if (label[-1] == "'" and label[-2] == "'") or label[-1] != "'":
                label = label + "'"
        return label

    validate_identifier = staticmethod(validate_identifier)


    #######################################################################
    ## INSTANCE METHODS

    def __init__(self, stream_handle=None):
        self._reset()
        if stream_handle:
            self.stream_handle = stream_handle
            
    def _reset(self):
        self.stream_handle = None
        self.current_file_char = None
        self.current_token = None
        self.eof = False
        self.current_line_number = 1
        self.current_col_number = 1
        self.previous_file_char = None

    def _get_current_file_char(self):
        "Returns the current character from the file stream."
        if self.__current_file_char == None:
            self.__current_file_char = self.read_next_char()
        return self.__current_file_char

    def _set_current_file_char(self, new_char):
        self.__current_file_char = new_char

    current_file_char = property(_get_current_file_char, _set_current_file_char)

    def tree_rooted_comment(self):
        """This is as hack it is called at the beginning of reading a tree. 
        Returns True if &R in comments, False in &U and None if no command comment
        is found """
        for c in self.comments:
            if c == '&U' or c == '&u':
                return False
            if c == '&R' or c == '&r':
                return True
        return False

    def read_next_char(self):
        """
        Advances the file stream cursor to the next character and returns
        it.
        """
        if self.stream_handle:
            read_char = self.stream_handle.read(1) # returns empty string if EOF
            if read_char == '':
                self.eof = True
            else:
                if self.previous_file_char == '\n':
                    self.current_line_number = self.current_line_number + 1
                    self.current_col_number = 0
                self.previous_file_char = self.__current_file_char
                self.current_col_number += 1
            self.current_file_char = read_char
            return self.__current_file_char
        return None

    def _raw_read_next_char(self):
            read_char = self.stream_handle.read(1) # returns empty string if EOF
            if read_char == '':
                raise StopIteration()
            if self.previous_file_char == '\n':
                self.current_line_number = self.current_line_number + 1
                self.current_col_number = 0
            self.previous_file_char = self.__current_file_char
            self.current_col_number += 1
            self.current_file_char = read_char
            return self.__current_file_char
        
    def skip_comment(self):
        """
        Reads characters from the file until current comment block (and
        any nested comment block terminates. Assumes current cursor
        position is on first character inside comment block.
        """
        cmt_body = StringIO()
        c = self.current_file_char
#         assert c == '['
        c = self.read_next_char()
        nesting = 1
        while not self.eof:
            if c == ']':
                if nesting == 1:
                    break
                nesting -= 1
            if c == '[':
                nesting += 1
            cmt_body.write(c)
            c = self.read_next_char()
        self.comments.append(cmt_body.getvalue())
        self.read_next_char()

    def read_noncomment_character(self):
        """
        Gets the first character outside a comment block from the
        current file stream position, inclusive.
        """
        if self.current_file_char == '[':
            self.skip_comment()
        return self.current_file_char
    noncomment_file_char = property(read_noncomment_character)

    def skip_to_significant_character(self):
        "Advances to the first non-whitespace character outside a comment block."
        while (PurePythonNexusStreamTokenizer.is_whitespace(self.current_file_char) or self.current_file_char=='[') and not self.eof:
            if self.current_file_char=='[':
                self.read_noncomment_character()
            else:
                self.read_next_char()
        return self.current_file_char


    def read_next_token(self, ignore_punctuation=None):
        """
        Reads the next token in the file stream. A token in this context is any word or punctuation character
        outside of a comment block.
        """
        self.comments = []
        if ignore_punctuation == None:
            ignore_punctuation = []
        self.current_token = None
        if self.eof:
            return None
        c = self.skip_to_significant_character()
        if self.eof:
            return None
        if c == "'":
            token = StringIO()
            try:
                fastfunc = self._raw_read_next_char
                c = fastfunc()
                while True:
                    if c == "'":
                        c = self.read_next_char()
                        if c == "'":
                            token.write("'")
                            c = fastfunc()
                        else:
                            break
                    else:
                        token.write(c)
                        c = fastfunc()
            except StopIteration:
                self.eof = True
                raise self.syntax_exception("Unexpected end of file inside quoted token")
            tokenstr = token.getvalue()
        else:
            quick_check = PurePythonNexusStreamTokenizer.is_punctuation
            if quick_check(c) and c not in ignore_punctuation:
                if c == '_':
                    c = ' '
                self.current_token = c
                self.read_next_char()
                return c
            token = StringIO()
            fget = self.stream_handle.read
            quick_check = PurePythonNexusStreamTokenizer.is_whitespace_or_punctuation

            if ignore_punctuation:
                while True:
                    if quick_check(c):
                        if c == '[':
                            self.skip_comment()
                            c = self.current_file_char
                            continue
                        if not c in ignore_punctuation:
                            break
                    if c == '_':
                        c = ' '
                    token.write(c)
                    prev = c
                    c = fget(1)
                    if not c:
                        break
            else:
                while True:
                    if quick_check(c):
                        if c == '[':
                            self.skip_comment()
                            c = self.current_file_char
                            continue
                        break
                    if c == '_':
                        c = ' '
                    token.write(c)
                    prev = c
                    c = fget(1)
                    if not c:
                        break
            tokenstr = token.getvalue()
            self.current_col_number += len(tokenstr)      
            self.current_file_char = c
            self.previous_file_char = prev
        
        self.current_token = tokenstr
        return tokenstr

    def read_next_token_ucase(self, ignore_punctuation=()):
        """
        Reads the next token in the file stream, upper-casing it
        before returning it.
        """
        t = self.read_next_token(ignore_punctuation=ignore_punctuation)
        if t != None:
            return t.upper()
        else:
            return None

    def skip_to_semicolon(self):
        "Advances the file stream cursor to the next semi-colon."
        token = self.read_next_token()
        while token != ';' and not self.eof and token != None:
            token = self.read_next_token()
            pass
    def syntax_exception(self, message):
            """
            Returns an exception object parameterized with line and
            column number values.
            """
            return SyntaxException(row=self.current_line_number, 
                                   column=self.current_col_number, 
                                   message=message)
############################################################################
## CLASS: NexusReader

class PurePythonNexusReader(datasets.Reader):
    "Encapsulates loading and parsing of a NEXUS format file."

    def __init__(self):
        datasets.Reader.__init__(self)
        self.stream_tokenizer = PurePythonNexusStreamTokenizer()
        self.encode_splits = False
        self.default_rooting = RootingInterpretation.UNKNOWN_DEF_ROOTED
        self.finish_node_func = None
        self._reset()

    def _reset(self):
        self.dataset = datasets.Dataset()
        self.char_block_type = characters.StandardCharactersBlock
        self.interleave = False
        self.symbols = "012"
        self.gap_char = '-'
        self.missing_char = '?'
        self.match_char = '.'
        self.tree_translate_dict = {}
        self.tax_label_lookup = {}

    def _prepare_to_read_file(self, file_obj):
        self.stream_tokenizer = PurePythonNexusStreamTokenizer()
        self.stream_tokenizer.stream_handle = file_obj

    def read_dataset(self, file_obj, dataset=None):
        """
        Instantiates and returns a DataSet object based on the
        NEXUS-formatted contents read from the file descriptor object
        `file_obj`.
        """
        #self._reset()
        self._prepare_to_read_file(file_obj)
        return self._parse_nexus_file(dataset)

    ## Class-specific ##

    def _parse_nexus_file(self, dataset=None):
        "Main file parsing driver."
        finish_node_func = self.finish_node_func
        self._reset()
        if dataset is not None:
            self.dataset = dataset
        token = self.stream_tokenizer.read_next_token_ucase()
        if token != "#NEXUS":
            raise self.syntax_exception('Expecting "#NEXUS", but found "%s"' % token)
        else:
            while not self.stream_tokenizer.eof:
                token = self.stream_tokenizer.read_next_token_ucase()
                while token != None and token != 'BEGIN' and not self.stream_tokenizer.eof:
                    token = self.stream_tokenizer.read_next_token_ucase()
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == 'TAXA':
                    self._parse_taxa_block()
                elif token == 'CHARACTERS':
                    if self.include_characters:
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

    def _consume_to_end_of_block(self, token):
        while not (token == 'END' or token == 'ENDBLOCK') \
            and not self.stream_tokenizer.eof \
            and not token==None:
            self.stream_tokenizer.skip_to_semicolon()
            token = self.stream_tokenizer.read_next_token_ucase()
        return token            
        
    def syntax_exception(self, message):
        """
        Returns an exception object parameterized with line and
        column number values.
        """
        return self.stream_tokenizer.syntax_exception(message)

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
                        self.char_block_type = characters.DnaCharactersBlock
                    elif token == "RNA":
                        self.char_block_type = characters.RnaCharactersBlock
                    elif token == "PROTEIN":
                        self.char_block_type = characters.ProteinCharactersBlock
                    else:                        
                        # defaults to STANDARD elif token == "STANDARD":
                        self.char_block_type = characters.StandardCharactersBlock
                        self.symbols = "12"
                else:
                    raise self.syntax_exception('Expecting "=" after DATATYPE keyword')
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
                        raise self.syntax_exception('Expecting \'"\' before beginning SYMBOLS list')
                else:
                    raise self.syntax_exception('Expecting "=" after SYMBOLS keyword')
            elif token == 'GAP':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    self.gap_char = token
                else:
                    raise self.syntax_exception('Expecting "=" after GAP keyword')
            elif token == 'INTERLEAVE':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token.upper().startswith("N"):
                        self.interleave = False
                    else:
                        self.interleave = True
                else:
                    self.interleave = True                    
            elif token == 'MISSING':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    self.missing_char = token
                else:
                    raise self.syntax_exception('Expecting "=" after MISSING keyword')
            elif token == 'MATCHCHAR':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    self.match_char = token
                else:
                    raise self.syntax_exception('Expecting "=" after MISSING keyword')
            token = self.stream_tokenizer.read_next_token_ucase()
         
    def _parse_tree_statement(self, taxa_block):
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
            raise self.syntax_exception('Expecting "=" in definition of Tree "%s" but found "%s"' % (tree_name, token))

        rooted = self.default_rooting
        if rooted == RootingInterpretation.UNKNOWN_DEF_ROOTED or rooted == RootingInterpretation.UNKNOWN_DEF_UNROOTED:
            for c in self.stream_tokenizer.comments:
                if c == '&U' or c == '&u':
                    rooted = RootingInterpretation.UNROOTED
                    break
                elif c == '&R' or c == '&r':
                    rooted = RootingInterpretation.ROOTED
                    break
        tree = parse_newick_tree_stream(stream_tokenizer=self.stream_tokenizer,
                                        taxa_block=taxa_block,
                                        translate_dict=self.tree_translate_dict,
                                        encode_splits=self.encode_splits,
                                        rooted=rooted,
                                        finish_node_func=self.finish_node_func)
        tree.label = tree_name

        if self.stream_tokenizer.current_token != ';':
            self.stream_tokenizer.skip_to_semicolon()
        return tree

    def _parse_translate_statement(self, taxa):
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
                t = taxa.get_taxon(label=translation_label)
            self.tree_translate_dict[translation_token] = t

            token = self.stream_tokenizer.read_next_token() # ","
            if (not token) or (token == ';'):
                break
            if token != ',':
                raise self.syntax_exception('Expecting "," in TRANSLATE statement after definition for %s = "%s", but found "%s" instead' % (translation_token, translation_label, token))

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
                        raise self.syntax_exception('Expecting numeric value for NTAX')
                else:
                    raise self.syntax_exception('Expecting "=" after NTAX keyword')
            elif token == 'NCHAR':
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == '=':
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token.isdigit():
                        self.file_specified_nchar = int(token)
                    else:
                        raise self.syntax_exception('Expecting numeric value for NCHAR')
                else:
                    raise self.syntax_exception('Expecting "=" after NCHAR keyword')
            token = self.stream_tokenizer.read_next_token_ucase()

    def _get_default_taxa_block(self, taxa_block=None):
        if taxa_block is None:
            if len(self.dataset.taxa_blocks) == 0:
                taxa_block = taxa.TaxaBlock()
                self.dataset.add_taxa_block(taxa_block=taxa_block)
            else:
                taxa_block = self.dataset.taxa_blocks[0]
        else:
            if taxa_block not in self.dataset.taxa_blocks:
                self.dataset.taxa_blocks.append(taxa_block)
        return taxa_block                          
          
    def _parse_taxlabels_statement(self, taxa_block=None):
        """
        Processes a TAXLABELS command. Assumes that the file reader is
        positioned right after the "TAXLABELS" token in a TAXLABELS command.
        """
        taxa_block = self._get_default_taxa_block(taxa_block)                
        token = self.stream_tokenizer.read_next_token()
        while token != ';':
            label = token
            taxa_block.add_taxon(label=label)        
            token = self.stream_tokenizer.read_next_token()

    def _build_state_alphabet(self, char_block, symbols):
        return build_state_alphabet(char_block, symbols, missing=self.missing_char)

    def _parse_matrix_statement(self, taxa_block=None):
        """
        Processes a MATRIX command. Assumes that the file reader
        is positioned right after the "MATRIX" token in a MATRIX command,
        and that NTAX and NCHAR have been specified accurately.
        """
        taxa_block = self._get_default_taxa_block(taxa_block)
        if not self.file_specified_ntax:
            raise self.syntax_exception('NTAX must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        elif not self.file_specified_nchar:
            raise self.syntax_exception('NCHAR must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        else:
            char_block = self.char_block_type()
            char_block.taxa_block = taxa_block
            if isinstance(char_block, characters.StandardCharactersBlock):
                char_block = self._build_state_alphabet(char_block, self.symbols)
            self.dataset.add_char_block(char_block=char_block)
            symbol_state_map = char_block.default_state_alphabet.symbol_state_map()
            if True: # future: trap and handle no labels, transpose etc.
                token = self.stream_tokenizer.read_next_token()
                while token != ';' and not self.stream_tokenizer.eof:
                    taxon = taxa_block.get_taxon(label=token)
                    if taxon not in char_block:
                        if self.include_characters:
                            char_block[taxon] = characters.CharacterDataVector(taxon=taxon)
                    if self.interleave:
                        while self.stream_tokenizer.current_file_char != '\n' and self.stream_tokenizer.current_file_char != '\r':
                            if self.stream_tokenizer.current_file_char not in [' ', '\t'] and self.include_characters:
                                state = symbol_state_map[self.stream_tokenizer.current_file_char]
                                char_block[taxon].append(characters.CharacterDataCell(value=state))
                            self.stream_tokenizer.read_next_char()
                    else:
                        while len(char_block[taxon]) < self.file_specified_nchar and not self.stream_tokenizer.eof:
                            char_group = self.stream_tokenizer.read_next_token()
                            if self.include_characters:
                                char_group = parse_sequence_iupac_ambiguities(char_group)
                                for char in char_group:
                                    state = symbol_state_map[char]
                                    char_block[taxon].append(characters.CharacterDataCell(value=state))
                    token = self.stream_tokenizer.read_next_token()
            else:
                ## TODO: NO LABELS/TRANSPOSED ##
                pass
                
    def iterate_over_trees(self, file_obj=None, taxa_block=None, dataset=None):
        """
        Generator to iterate over trees in data file.
        Primary goal is to be memory efficient, storing no more than one tree
        at a time. Speed might have to be sacrificed for this!
        """
        if dataset is None:
            dataset = datasets.Dataset()
        if taxa_block is None:
            taxa_block = taxa.TaxaBlock()
        if not (taxa_block in dataset.taxa_blocks):
            dataset.taxa_blocks.append(taxa_block)
        stream_tokenizer = PurePythonNexusStreamTokenizer(file_obj)
        self.stream_tokenizer = stream_tokenizer
        while not self.stream_tokenizer.eof:
            token = self.stream_tokenizer.read_next_token_ucase()
            while token != None and token != 'BEGIN' and not self.stream_tokenizer.eof:
                token = self.stream_tokenizer.read_next_token_ucase()
            token = self.stream_tokenizer.read_next_token_ucase()
            if token == 'TAXA':
                self._parse_taxa_block(taxa_block)
            elif token == 'TREES':
                self._prepare_to_parse_trees(taxa_block)
                self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
                while not (token == 'END' or token == 'ENDBLOCK') \
                    and not self.stream_tokenizer.eof \
                    and not token==None:
                    token = self.stream_tokenizer.read_next_token_ucase()
                    if token == 'TRANSLATE':
                        self._parse_translate_statement(taxa_block)
                    if token == 'TREE':
                        tree = self._parse_tree_statement(taxa_block)
                        yield tree
                self.stream_tokenizer.skip_to_semicolon() # move past END command
            else:
                # unknown block
                while not (token == 'END' or token == 'ENDBLOCK') \
                    and not self.stream_tokenizer.eof \
                    and not token==None:
                    self.stream_tokenizer.skip_to_semicolon()
                    token = self.stream_tokenizer.read_next_token_ucase()

    def _parse_taxa_block(self, taxa_block=None):
        token = ''
        self.stream_tokenizer.skip_to_semicolon() # move past BEGIN statement
        while not (token == 'END' or token == 'ENDBLOCK') \
            and not self.stream_tokenizer.eof \
            and not token==None:
            token = self.stream_tokenizer.read_next_token_ucase()
            if token == 'DIMENSIONS':
                self._parse_dimensions_statement()
            if token == 'TAXLABELS':
                self._parse_taxlabels_statement(taxa_block)
        self.stream_tokenizer.skip_to_semicolon() # move past END statement

    def _parse_trees_block(self):
        token = 'TREES'
        if self.include_trees:
            trees_block = trees.TreesBlock()
            trees_block.taxa_block = self._get_default_taxa_block()
            self.dataset.add_trees_block(trees_block=trees_block)
            self._prepare_to_parse_trees(trees_block.taxa_block)
            self.stream_tokenizer.skip_to_semicolon() # move past BEGIN command
            while not (token == 'END' or token == 'ENDBLOCK') \
                and not self.stream_tokenizer.eof \
                and not token==None:
                token = self.stream_tokenizer.read_next_token_ucase()
                if token == 'TRANSLATE':
                    self._parse_translate_statement(trees_block.taxa_block)
                if token == 'TREE':
                    tree = self._parse_tree_statement(trees_block.taxa_block)
                    trees_block.append(tree)
            self.stream_tokenizer.skip_to_semicolon() # move past END command
        else:
            token = self.consume_to_end_of_block(token)

    def _prepare_to_parse_trees(self, taxa_block):

            self.tree_translate_dict = {}
            self.tax_label_lookup = {}
            for n, t in enumerate(taxa_block):
                self.tree_translate_dict[str(n + 1)] = t
            # add labels second so that numbers have priority over number
            for n, t in enumerate(taxa_block):
                l = t.label
                self.tree_translate_dict[l] = t
                self.tax_label_lookup[l] = t
                if self.encode_splits:
                    ti = taxa_block.index(t)
                    t.clade_mask = (1 << ti)


############################################################################
## CLASS: NexusWriter



class NexusWriter(datasets.Writer):
    "Implements the DataWriter interface for handling NEXML files."

    def __init__(self, simple=False):
        "Calls the base class constructor."
        datasets.Writer.__init__(self)
        self.simple = simple
        self.write_rooting = True
        self.comment = []

    def escape_token(label):
        if PurePythonNexusStreamTokenizer.has_whitespace(label) or PurePythonNexusStreamTokenizer.has_punctuation(label):
            return "'" + label + "'"
        else:
            return label
    escape_token = staticmethod(escape_token)

    def write_dataset(self, dataset, dest):
        """
        Writes dataset to a full NEXUS
        document.
        """
        dest.write('#NEXUS\n\n')

        if self.comment is not None:
            if isinstance(self.comment, list):
                for line in self.comment:
                    if line.strip().replace("\n", "").replace("\r", ""):
                        dest.write("[ %s ]\n" % line)
                    else:
                        dest.write("\n")
                dest.write("\n")
            else:
                dest.write("[ %s ]\n\n" % self.comment)
        if (dataset.char_blocks or dataset.trees_blocks) and not self.simple:
            self.write_taxa_block(taxa_block=dataset.taxa_blocks[0], dest=dest)
        if dataset.char_blocks:
            self.write_char_block(char_block=dataset.char_blocks[0], dest=dest, simple_nexus=self.simple)
        if dataset.trees_blocks:
            self.write_trees_block(trees_block=dataset.trees_blocks[0], dest=dest)

    def write_trees_block(self, trees_block, dest):
        block = []
        newick_writer = NewickWriter()
        block.append('BEGIN TREES;')
        for treeidx, tree in enumerate(trees_block):
            if tree.label:
                tree_name = tree.label
            else:
                tree_name = str(treeidx)
            newick_str = newick_writer.compose_node(tree.seed_node)
            if tree.is_rooted and self.write_rooting:
                rooting = "[&R] "
            elif not tree.is_rooted and self.write_rooting:
                rooting = "[&U] "
            else:
                rooting = ""
            block.append('    tree %s = %s%s;' % (self.escape_token(tree_name), rooting, newick_str))
        block.append('END;\n\n')
        dest.write('\n'.join(block))

    def write_taxa_block(self, taxa_block, dest):
        block = []
        block.append('BEGIN TAXA;')
        block.append('    DIMENSIONS NTAX=%d;' % len(taxa_block))
        block.append('    TAXLABELS')
        for taxon in taxa_block:
            block.append('        %s' % self.escape_token(taxon.label))
        block.append('  ;')
        block.append('END;\n\n')
        dest.write('\n'.join(block))

    def compose_format_terms(self, char_block):
        format = []
        if isinstance(char_block, characters.DnaCharactersBlock):
            format.append("DATATYPE=DNA")
        elif isinstance(char_block, characters.RnaCharactersBlock):
            format.append("DATATYPE=RNA")
        elif isinstance(char_block, characters.ProteinCharactersBlock):
            format.append("DATATYPE=PROTEIN")
        else:
            format.append("DATATYPE=STANDARD")
            format.append('SYMBOLS="%s"' % ''.join([str(a) for a in char_block.default_state_alphabet]))
        format.append("GAP=- MISSING=? MATCHCHAR=.")
        return ' '.join(format)

    def write_char_block(self, char_block, dest, simple_nexus=False):
        nexus = []
        taxlabels = [self.escape_token(taxon.label) for taxon in char_block]
        max_label_len = max([len(label) for label in taxlabels])
        nchar = max([len(seq) for seq in char_block.values()])
        if simple_nexus:
            nexus.append('BEGIN DATA;')
            ntaxstr = "NTAX=%d" % len(taxlabels)
        else:
            nexus.append('BEGIN CHARACTERS;')
            ntaxstr = ""
        nexus.append('    DIMENSIONS %s NCHAR=%d;' % (ntaxstr, nchar))
        nexus.append('    FORMAT %s;' % self.compose_format_terms(char_block))
        nexus.append('    MATRIX')
        taxa = char_block.keys()
        taxa.sort()
        for taxon in taxa:
            seq_vec = char_block[taxon]
            seq = ''.join([str(seq_sym) for seq_sym in seq_vec])
#             seq.replace('~','-')
            seq = seq.ljust(nchar, '-')
            nexus.append('%s    %s' % (self.escape_token(taxon.label).ljust(max_label_len), seq))
        nexus.append('    ;')
        nexus.append('END;\n\n')
        dest.write('\n'.join(nexus))
                
############################################################################
## CLASS: NewickReader        
        
class NewickReader(datasets.Reader):
    "Implementation of TreeReader for NEWICK files and strings."
    
    def __init__(self):
        """
        This has changed so much so many times that any documentation
        I put in here will probably be obselete in a matter of hours
        so this comment is all you are getting.
        """
        datasets.Reader.__init__(self)
        self.encode_splits = False
        self.default_rooting = RootingInterpretation.UNKNOWN_DEF_UNROOTED
        self.finish_node_func = None

        
    def read_dataset(self, file_obj, dataset=None):
        """
        Instantiates and returns a DataSet object based on the
        NEWICK-formatted contents read from the file descriptor object
        `file_obj`.
        """
        taxa_block = None
        need_to_add_taxa = True
        if dataset is None:
            dataset = datasets.Dataset()
        elif len(dataset.taxa_blocks) == 1:
            taxa_block = dataset.taxa_blocks[0]
            need_to_add_taxa = False
        if taxa_block is None:
            taxa_block = taxa.TaxaBlock()
        trees_block = self.read_trees(file_obj, taxa_block=taxa_block)
        if need_to_add_taxa:
            dataset.add_taxa_block(taxa_block=taxa_block)
        dataset.add_trees_block(trees_block=trees_block, taxa_block=taxa_block)
        return dataset

    def read_trees(self, file_obj=None, taxa_block=None):
        """
        Instantiates and returns a TreesBlock object based
        on the Newick-formatted contents read from the file
        descriptor object `file_obj`.
        """
        trees_block = trees.TreesBlock()                    
        if taxa_block is not None:
            trees_block.taxa_block = taxa_block
        stream_tokenizer = PurePythonNexusStreamTokenizer(file_obj)
        tree = parse_newick_tree_stream(stream_tokenizer, taxa_block=taxa_block, encode_splits=self.encode_splits, rooted=self.default_rooting, finish_node_func=self.finish_node_func)
        if taxa_block is None:
            taxa_block = tree.taxa_block
            trees_block.taxa_block = taxa_block
        while tree is not None:
            trees_block.append(tree)
            tree = parse_newick_tree_stream(stream_tokenizer, taxa_block=taxa_block, encode_splits=self.encode_splits, rooted=self.default_rooting, finish_node_func=self.finish_node_func)
        return trees_block

    def iterate_over_trees(self, file_obj=None, taxa_block=None, dataset=None):
        """
        Generator to iterate over trees in data file.
        Primary goal is to be memory efficient, storing no more than one tree
        at a time. Speed might have to be sacrificed for this!
        """
        if dataset is None:
            dataset = datasets.Dataset() or dataset
        if taxa_block is None:
            taxa_block = taxa.TaxaBlock()
        if not (taxa_block in dataset.taxa_blocks):
            dataset.taxa_blocks.append(taxa_block)
        stream_tokenizer = PurePythonNexusStreamTokenizer(file_obj)
        token = stream_tokenizer.read_next_token_ucase()
        stream_tokenizer.stream_handle.seek(0)
        while not stream_tokenizer.eof:
            tree = parse_newick_tree_stream(stream_tokenizer=stream_tokenizer, 
                                                 taxa_block=taxa_block, 
                                                 translate_dict=None,
                                                 encode_splits=self.encode_splits,
                                                 rooted=self.default_rooting,
                                                 finish_node_func=self.finish_node_func) 
            if tree:
                yield tree

############################################################################
## CLASS: NewickWriter

class NewickWriter(datasets.Writer):
    """
    Handles representation and serialization of a DendroPy Tree object
    in NEWICK format.
    """

    def __init__(self, **kwargs):
        """
        Instantiates the object, setting default for various
        formatting/representation options.
        """
        self.edge_lengths = True
        self.internal_labels = True
        self.support_as_edge_lengths = False
        self.support_as_labels = False
        self.support_as_percentages = False
        self.support_decimals = None

    def write_dataset(self, dataset, dest):
        """
        Writes a DataSet object to a full document-level
        representation of the format being implemented by the
        deriving class. 
        """
        for trees_block in dataset.trees_blocks:
            for tree in trees_block:
                dest.write(self.compose_node(tree.seed_node) + ';\n')                                
                
    def compose_taxlabel(self, label):
        if re.search('[' + PurePythonNexusStreamTokenizer.whitespace + PurePythonNexusStreamTokenizer.punctuation + ']', label) != None:
            return "'" + label + "'"
        else:
            return label    
    
    def compose_tree(self, tree):
        "Convienience method.        "
        return self.compose_node(tree.seed_node)

    def choose_display_tag(self, node):
        """
        Based on current settings, the attributes of a node, and
        whether or not the node is a leaf, returns an appropriate tag.
        """
        if hasattr(node, 'taxon') and node.taxon:
            return self.compose_taxlabel(node.taxon.label)
        elif hasattr(node, 'label') and node.label:
            return self.compose_taxlabel(node.label)
        elif len(node.child_nodes()) == 0:
            # force label if a leaf node
            return self.compose_taxlabel(node.oid)
        else:
            return ""
        
    def compose_node(self, node):
        """
        Given a DendroPy Node, this returns the Node as a NEWICK
        statement according to the class-defined formatting rules.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            subnodes = [self.compose_node(child) for child in child_nodes]
            statement = '(' + ','.join(subnodes) + ')'
            if self.internal_labels:
                statement = statement + self.choose_display_tag(node)
            if node.edge and node.edge.length != None and self.edge_lengths:
                try:
                    statement =  "%s:%f" \
                                % (statement, float(node.edge.length))
                except ValueError:
                    statement =  "%s:%s" \
                                % (statement, node.edge.length)
            return statement
        else:
            if self.internal_labels:
                statement = self.choose_display_tag(node)
            if node.edge and node.edge.length != None and self.edge_lengths:
                try:
                    statement =  "%s:%0.10f" \
                                % (statement, float(node.edge.length))
                except ValueError:
                    statement =  "%s:%s" \
                                % (statement, node.edge.length)
            return statement
            
try:
    import nclwrapper
except:
    NexusStreamTokenizer = PurePythonNexusStreamTokenizer
    NexusReader = PurePythonNexusReader
else:
    from threading import Thread, Event

    _LOG.debug("Using NCLNexusReader")


    class NCLTreeStream(nclwrapper.NxsTreeStream):
        """Simple thread-safe class that waits for `need_tree_event', and signals the
        presence of a new tree by `ready_event`"""
        def __init__(self, need_tree_event, ready_event):
            self.need_tree_event = need_tree_event 
            self.ready_event = ready_event
            self.tree_tokens = None
            self.taxa_block =  None
            nclwrapper.NxsTreeStream.__init__(self)

        def handleTree(self, ftd, tb):
            t = ftd.GetTreeTokens()
            rooted_flag = ftd.IsRooted()
            self.need_tree_event.wait()
            self.need_tree_event.clear()
            try:
                self.ncl_taxa_block =  tb.GetTaxaBlockPtr()
                tb_iid = self.ncl_taxa_block.GetInstanceIdentifierString()
                self.tree_tokens =  t
                self.rooted_flag = rooted_flag
            except Exception, v:
                _LOG.warn("Exception: %s" % str(v))
            self.ready_event.set()
            return False

    class NCLTreeStreamThread(Thread):
        def __init__(self, file_path, need_tree_event, ready_event, format="NEXUS", **kwargs):
            """Subclass of thread ,that uses a NCLTreeStream to get trees
            from NCL one at a time"""
            self.nts = NCLTreeStream(need_tree_event, ready_event)
            self.file_path = file_path
            self.format = format
            self.exception = None
            self.reader = nclwrapper.MultiFormatReader()
            self.reader.cullIdenticalTaxaBlocks(True)
            Thread.__init__(self, 
                            group=None, 
                            target=None, 
                            name=None, 
                            args=tuple(), 
                            kwargs=dict(**kwargs))
        def run(self):
            self.done = False
            try:
                self.nts.ReadFilepath(self.file_path, self.format, self.reader)
            except ValueError, v:
                self.exception = v
            self.nts.need_tree_event.wait()
            self.nts.need_tree_event.clear()
            self.done = True
            self.nts.tree_tokens = None
            self.nts.taxa_block = None
            self.nts.ready_event.set()

    def _ncl_datatype_enum_to_dendropy(d):
        e = nclwrapper.NxsCharactersBlock
        if d == e.dna or d == e.nucleotide:
            return characters.DnaCharactersBlock
        if d == e.rna:
            return characters.RnaCharactersBlock
        if d == e.protein:
            return characters.ProteinCharactersBlock
        if (d == e.continuous) or (d == e.mixed) or (d == e.codon):
            s = d == e.continuous and "continuous" or (d == e.mixed and "mixed" or "codon")
            raise NotImplementedError("%s datatype not supported" % s)

    class ListOfTokenIterator(object):
        def __init__(self, tokens):
            self.tokens_iter = iter(tokens)
            self.eof = False
            self.queued = None
            self.tree_rooted = None

        def tree_rooted_comment(self):
            "This is a hack, and only works if you have just one tree in the token stream"
            return self.tree_rooted

        def __iter__(self):
            return self
        def read_next_token(self, ignore_punctuation=None):
            if not self.eof:
                try:
                    return self.tokens_iter.next()
                except StopIteration:
                    self.eof = True
            return ""

        def read_next_token_ucase(self):
            t = self.read_next_token()
            if t:
                return t.upper()
        def syntax_exception(self, msg):
            return SyntaxException(message=msg)
        
    class NCLBasedReader(datasets.Reader):
        "Encapsulates loading and parsing of a NEXUS format file."
    
        def __init__(self, format="NEXUS"):
            datasets.Reader.__init__(self)
            self.purePythonReader = PurePythonNexusReader()
            self.encode_splits = False
            self.default_rooting = RootingInterpretation.UNKNOWN_DEF_ROOTED
            self.finish_node_func = None
            self.format = format
            self._prev_taxa_block = None
            self.ncl_taxa_to_native = {}
            self._taxa_to_fill = None
        def _get_fp(self, file_obj):
            "Returns filepath and True if the file that `file_obj` refers to exists on the filesystem"
            try:
                n = file_obj.name
                use_ncl = os.path.exists(n)
                return n, use_ncl
            except AttributeError:
                return "", False
            
        def read_dataset(self, file_obj, dataset=None):
            """
            Instantiates and returns a DataSet object based on the
            NEXUS-formatted contents read from the file descriptor object
            `file_obj`.
            """
            n, use_ncl = self._get_fp(file_obj)
            if not use_ncl:
                self.purePythonReader.encode_splits = self.encode_splits
                self.purePythonReader.default_rooting = self.default_rooting
                self.purePythonReader.finish_node_func = self.finish_node_func
                return self.purePythonReader.read_dataset(file_obj, dataset=dataset)
            return self.read_filepath_into_dataset(n, dataset=dataset)

        def _ncl_characters_block_to_native(self, taxa_block, ncl_cb):
            """
            Processes a FORMAT command. Assumes that the file reader is
            positioned right after the "FORMAT" token in a FORMAT command.
            """
            raw_matrix = ncl_cb.GetRawDiscreteMatrixRef()
            if ncl_cb.IsMixedType():
                _LOG.warn("Mixed datatype character blocks are not supported in Dendropy.  Skipping...")
                return None
            char_block_type = _ncl_datatype_enum_to_dendropy(ncl_cb.GetDataType())
            mapper = ncl_cb.GetDatatypeMapperForCharRef(0)
            symbols = mapper.GetSymbols()
            state_codes_mapping = mapper.GetPythonicStateVectors()
                
            char_block = char_block_type()
            char_block.taxa_block = taxa_block
            if isinstance(char_block, characters.StandardCharactersBlock):
                char_block = self.build_state_alphabet(char_block, symbols, '?')
            symbol_state_map = char_block.default_state_alphabet.symbol_state_map()
            
            ncl_numeric_code_to_state = []
            for s in symbols:
                ncl_numeric_code_to_state.append(symbol_state_map[s])
            for sc in state_codes_mapping[len(symbols):-2]:
                search = set()
                for fundamental_state in sc:
                    search.add(ncl_numeric_code_to_state[fundamental_state])
                found = False
                for sym, state in symbol_state_map.iteritems():
                    ms = state.member_states
                    if ms:
                        possible = set(ms)
                        if possible == search:
                            found = True
                            ncl_numeric_code_to_state.append(state)
                            break
                if not found:
                    raise ValueError("NCL datatype cannot be coerced into datatype because ambiguity code for %s is missing " % str(search))
            ncl_numeric_code_to_state[-2] = symbol_state_map['-']
            ncl_numeric_code_to_state[-1] = symbol_state_map['?']
                
                
            
            assert (len(raw_matrix) == len(taxa_block))
            for row_ind, taxon in enumerate(taxa_block):
                v = characters.CharacterDataVector(taxon=taxon)
                raw_row = raw_matrix[row_ind]
                char_block[taxon] = v
                if self.include_characters:
                    for c in raw_row:
                        state = ncl_numeric_code_to_state[c]
                        v.append(characters.CharacterDataCell(value=state))
    
            #dataset.characters_blocks.append(char_block)
            supporting_exsets = False
            supporting_charset_exsets = False
    
            if supporting_exsets:
                s = ncl_cb.GetExcludedIndexSet()
                print "Excluded chars =", str(nclwrapper.NxsSetReader.GetSetAsVector(s))
            if supporting_charset_exsets:
                nab = m.GetNumAssumptionsBlocks(ncl_cb)
                for k in xrange(nab):
                    a = m.GetAssumptionsBlock(ncl_cb, k)
                    cs = a.GetCharSetNames()
                    print "CharSets have the names " , str(cs)
            return char_block

        def read_filepath_into_dataset(self, file_path, dataset=None):
            if dataset is None:
                dataset = datasets.Dataset()
            self._taxa_to_fill = None
            m = nclwrapper.MultiFormatReader()
            m.cullIdenticalTaxaBlocks(True)

            self._register_taxa_context(m, dataset.taxa_blocks)

            m.ReadFilepath(file_path, self.format)
            
            num_taxa_blocks = m.GetNumTaxaBlocks()
            for i in xrange(num_taxa_blocks):
                ncl_tb = m.GetTaxaBlock(i)
                taxa_block = self._ncl_taxa_block_to_native(ncl_tb)
                dataset.add_taxa_block(taxa_block=taxa_block)
                #nab = m.GetNumAssumptionsBlocks(ncl_tb)
                #for k in xrange(nab):
                #    a = m.GetAssumptionsBlock(ncl_tb, k)
                #    cs = a.GetTaxSetNames()
                #    print "TaxSets have the names " , str(cs)
                num_char_blocks = m.GetNumCharactersBlocks(ncl_tb)
                for j in xrange(num_char_blocks):
                    ncl_cb = m.GetCharactersBlock(ncl_tb, j)
                    char_block = self._ncl_characters_block_to_native(taxa_block, ncl_cb)
                    if char_block:
                        dataset.add_char_block(char_block=char_block)
                ntrb = m.GetNumTreesBlocks(ncl_tb)
                for j in xrange(ntrb):
                    trees_block = trees.TreesBlock()
                    trees_block.taxa_block = taxa_block
                    ncl_trb = m.GetTreesBlock(ncl_tb, j)
                    for k in xrange(ncl_trb.GetNumTrees()):
                        ftd = ncl_trb.GetFullTreeDescription(k)
                        tokens = ftd.GetTreeTokens()
                        rooted_flag = ftd.IsRooted()
                        t = self._ncl_tree_tokens_to_native_tree(ncl_tb, taxa_block, tokens, rooted_flag=rooted_flag)
                        if t:
                            trees_block.append(t)
                    dataset.add_trees_block(trees_block=trees_block)
            return dataset


        def iterate_over_trees(self, file_obj=None, taxa_block=None, dataset=None):
            """
            Generator to iterate over trees in data file.
            Primary goal is to be memory efficient, storing no more than one tree
            at a time. Speed might have to be sacrificed for this!
            """
            if taxa_block is not None and len(taxa_block) == 0:
                self._taxa_to_fill = taxa_block
            else:
                self._taxa_to_fill = None
            n, use_ncl = self._get_fp(file_obj)
            if not use_ncl:
                self.purePythonReader.encode_splits = self.encode_splits
                self.purePythonReader.default_rooting = self.default_rooting
                self.purePythonReader.finish_node_func = self.finish_node_func
                for tree in self.purePythonReader.iterate_over_trees(file_obj, taxa_block=taxa_block, dataset=dataset):
                    yield tree
                return                
            if dataset is None:
                dataset = datasets.Dataset()
            if taxa_block is None:
                taxa_block = taxa.TaxaBlock()
            if taxa_block and not (taxa_block in dataset.taxa_blocks):
                dataset.taxa_blocks.append(taxa_block)

            need_tree_event = Event()
            tree_ready_event = Event()
            ntst = NCLTreeStreamThread(n, need_tree_event=need_tree_event, ready_event=tree_ready_event, format=self.format)

            ncl_streamer = ntst.nts
            
            self._register_taxa_context(ntst.reader, dataset.taxa_blocks)
            
            ntst.start()
            need_tree_event.set()
            while True:
                if ntst.done:
                    break
                tree_ready_event.wait()
                tree_ready_event.clear()
                ncl_taxa_block = ncl_streamer.ncl_taxa_block
                self.curr_tree_tokens = ncl_streamer.tree_tokens
                rooted_flag = ncl_streamer.rooted_flag
                ncl_streamer.tree_tokens = None
                need_tree_event.set()
                self.curr_tree = self._ncl_tree_tokens_to_native_tree(ncl_taxa_block, None, self.curr_tree_tokens, rooted_flag=rooted_flag)
                if self.curr_tree:
                    yield self.curr_tree
            del self.curr_tree_tokens
            del self.curr_tree
                    
    

        def _register_taxa_context(self, ncl_reader, incoming_taxa_blocks):
            if not incoming_taxa_blocks:
                return
            num_taxa_blocks = ncl_reader.GetNumTaxaBlocks()
            existing_taxa_blocks = []
            for i in xrange(num_taxa_blocks):
                ncl_tb = ncl_reader.GetTaxaBlock(i)
                labels = list(ncl_tb.GetAllLabels())
                existing_taxa_blocks.append(ncl_tb, labels)
            
            to_add = []
            for tb in incoming_taxa_blocks:
                found = False
                l = [i.label for i in tb]
                for k, v in existing_taxa_blocks:
                    if l == v:
                        found = True
                        self.ncl_taxa_to_native[k.GetInstanceIdentifierString()] = tb
                        break
                if not found:
                    to_add.append(tb)

            for tb in to_add:
                tn = tuple([i.label for i in tb])
                ncl_tb = ncl_reader.RegisterTaxa(tn)
                iid = ncl_tb.GetInstanceIdentifierString()
                self.ncl_taxa_to_native[iid] = tb

        def _ncl_taxa_block_to_native(self, ncl_tb):
            tbiid = ncl_tb.GetInstanceIdentifierString()
            taxa_block = self.ncl_taxa_to_native.get(tbiid)
            if taxa_block is not None:
                return taxa_block
                
            labels = ncl_tb.GetAllLabels()
            if self._taxa_to_fill is None:
                taxa_block =  taxa.TaxaBlock(labels)
            else:
                taxa_block = self._taxa_to_fill
                self._taxa_to_fill = None
                taxa_block.extend(labels)
            self.ncl_taxa_to_native[tbiid] = taxa_block
            return taxa_block
            
        def _ncl_tree_tokens_to_native_tree(self, ncl_tb, taxa_block, tree_tokens, rooted_flag=None):
            if not tree_tokens:
                return None
            if taxa_block is None:
                iid = ncl_tb.GetInstanceIdentifierString()
                taxa_block = self._ncl_taxa_block_to_native(ncl_tb)
                    
            self.taxa_block = taxa_block        
            lti = ListOfTokenIterator(tree_tokens)
            lti.tree_rooted = rooted_flag
            if not self._prev_taxa_block is taxa_block:
                self.tree_translate_dict = {}
                for n, t in enumerate(taxa_block):
                    self.tree_translate_dict[str(n + 1)] = t
                    if self.encode_splits:
                        t.clade_mask = (1 << n)
                self._prev_taxa_block = taxa_block

            return parse_newick_tree_stream(lti, 
                                            taxa_block=taxa_block,
                                            translate_dict=self.tree_translate_dict,
                                            encode_splits=self.encode_splits,
                                            rooted=self.default_rooting,
                                            finish_node_func=self.finish_node_func)
            

        
    NexusReader = NCLBasedReader
    
