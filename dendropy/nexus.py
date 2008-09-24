#! /usr/bin/env python

############################################################################
##  nexus.py
##
##  Part of the DendroPy phylogenetic tree manipulation library.
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
NEXUS format.
"""

import re
import os
import sys
import StringIO

from dendropy import datasets
from dendropy import taxa
from dendropy import trees
from dendropy import characters
from dendropy import newick

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
          result = result + map_to_ambiguity_code(match.group(1))
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

class NexusReader(datasets.Reader):
    """
    Encapsulates loading and parsing of a NEXUS format file.
    """    
    
    ##########################################################
    ## CHARACTER/SYMBOL ANALYSIS AND VALIDATION             ##
    ##########################################################
    
    punctuation = '\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>'
    whitespace = ' \0\t\n\r'
    
    def is_whitespace(char):
        return char in NexusReader.whitespace

    is_whitespace = staticmethod(is_whitespace)
    
    def has_whitespace(s):
        return re.search('['+NexusReader.whitespace+']', s) != None

    has_whitespace = staticmethod(has_whitespace)
    
    def is_punctuation(char):
        if char:
            return char in NexusReader.punctuation

    is_punctuation = staticmethod(is_punctuation)    
    
    def has_punctuation(s):
        return re.search('['+NexusReader.punctuation+']', s) != None

    has_punctuation = staticmethod(has_punctuation)
    
    def is_whitespace_or_punctuation(char):
        return NexusReader.is_whitespace(char) or NexusReader.is_punctuation(char)

    is_whitespace_or_punctuation = staticmethod(is_whitespace_or_punctuation)
    
    def has_whitespace_or_punctuation(s):
        return NexusReader.has_whitespace(s) or NexusReader.has_punctuation(s)

    has_whitespace_or_punctuation = staticmethod(has_whitespace_or_punctuation)

    def validate_identifier(label):
        if NexusReader.has_whitespace_or_punctuation(label):
            if (label[0] == "'" and label[1] == "'") or label[0] != "'":
                label = "'" + label
            if (label[-1] == "'" and label[-2] == "'") or label[-1] != "'":
                label = label + "'"
            return label
        else:
            return label

    validate_identifier = staticmethod(validate_identifier)
    
    ##########################################################
    ## EXCEPTIONS                                           ##
    ##########################################################    
    
    class SyntaxException(Exception):

        def __init__(self, row, column, message):
            self.row = row
            self.column = column
            self.message = message

        def __str__(self):
            return 'ERROR PARSING FILE IN LINE %d: %s' % (self.row, self.message) 

    class NotNexusFileException(SyntaxException):
        def __init__(self, filepath, row, column, message):
            super(NotNexusFileException, self).__init__(filepath, row, column, message)

    ##########################################################
    ## CLASS MEMBER METHODS                                 ##
    ##########################################################

    def __init__(self):
        """
        `tree_factory` is a DendroPy TreeFactory class or derived
        object.
        """
        datasets.Reader.__init__(self)
        self.dataset = None

    ## Implementation of the datasets.Reader interface ##

    def read_dataset(self, fileobj, dataset=None):
        """
        Instantiates and returns a DataSet object based on the
        NEXML-formatted contents read from the file descriptor object
        `fileobj`.
        """
        self.filehandle = fileobj
        return self.parse_nexus_file(dataset)

    ## Class-specific ##
    
    def tree_iter(self, filepath=None, fileobj=None, text=None): 
        """
        Iterates through trees in file, returning them one-by-one instead of
        parsing the entire dataset.
        """
        self.filehandle = datasets.Reader.get_file_handle(filepath=filepath, fileobj=fileobj, text=text)
        self.dataset = datasets.Dataset()
        self.interleave = False
        self.char_block_type = characters.DnaCharactersBlock
        self.gap_char = '-'
        self.missing_char = '?'
        self.match_char = '.'
        self.symbols = None
        self.current_file_char = None
        self.eof = False
        self.current_line_number = 1
        self.current_col_number = 1
        self.previous_file_char = None        
        self.tree_translate_dict = {}
        token = self.read_next_token_ucase()
        if token != "#NEXUS":
            ### Assume nexus tree file ###
            self.filehandle.seek(0)
            statement_block = self.filehandle.read()
            statement_block = statement_block.replace('\n','').replace('\r','')
            trees_block = self.trees_block_factory()
            for statement in statement_block.split(';'):
                statement = statement.strip() + ';'
                newick_parser = newick.NewickTreeParser()
                tree = newick_parser.parse_tree_statement(statement, trees_block)
                trees_block.pop()
                yield tree        
        else:
            while not self.eof:
                token = self.read_next_token_ucase()
                while token != None and token != 'BEGIN' and not self.eof:
                    token = self.read_next_token_ucase()
                token = self.read_next_token_ucase()
                if token == 'TREES':
                    trees_block = self.trees_block_factory()
                    trees_block.taxa_block = self.get_default_taxa_block()
                    self.dataset.add_trees_block(trees_block=trees_block)
                    self.skip_to_semicolon() # move past BEGIN command
                    while not (token == 'END' or token == 'ENDBLOCK') and not self.eof and not token==None:
                        token = self.read_next_token_ucase()
                        if token == 'TRANSLATE':
                            self.parse_translate_statement()                         
                        if token == 'TREE':
                            tree = self.parse_tree_statement(trees_block)  
                            trees_block.pop()
                            yield tree
                    self.skip_to_semicolon() # move past END command    
                else:
                    # unknown block
                    while not (token == 'END' or token == 'ENDBLOCK') and not self.eof and not token==None:
                        #print token
                        self.skip_to_semicolon()
                        token = self.read_next_token_ucase()   
        
    def parse_nexus_file(self, dataset=None):
        """
        Main file parsing driver.
        """
        if dataset is None:
            self.dataset = datasets.Dataset()
        else:
            self.dataset = dataset
        self.interleave = False
        self.char_block_type = characters.StandardCharactersBlock        
        self.symbols = "012"        
        self.gap_char = '-'
        self.missing_char = '?'
        self.match_char = '.'
        self.current_file_char = None
        self.eof = False
        self.current_line_number = 1
        self.current_col_number = 1
        self.previous_file_char = None        
        self.tree_translate_dict = {}
        token = self.read_next_token_ucase()
        if token != "#NEXUS":
            raise self.syntax_exception('Expecting "#NEXUS", but found "%s"' % token)
        else:
            while not self.eof:
                token = self.read_next_token_ucase()
                while token != None and token != 'BEGIN' and not self.eof:
                    token = self.read_next_token_ucase()
                token = self.read_next_token_ucase()
                if token == 'TAXA':
                    self.skip_to_semicolon() # move past BEGIN statement
                    while not (token == 'END' or token == 'ENDBLOCK') and not self.eof and not token==None:
                        token = self.read_next_token_ucase()
                        if token == 'DIMENSIONS':
                            self.parse_dimensions_statement()                        
                        if token == 'TAXLABELS':
                            self.parse_taxlabels_statement()
                    self.skip_to_semicolon() # move past END statement
                elif token == 'CHARACTERS':
                    self.skip_to_semicolon() # move past BEGIN command
                    while not (token == 'END' or token == 'ENDBLOCK') and not self.eof and not token==None:
                        token = self.read_next_token_ucase()
                        if token == 'DIMENSIONS':
                            self.parse_dimensions_statement()    
                        if token == 'FORMAT':
                            self.parse_format_statement()                              
                        if token == 'MATRIX':
                            self.parse_matrix_statement()
                    self.skip_to_semicolon() # move past END command
                elif token == 'DATA':
                    self.skip_to_semicolon() # move past BEGIN command
                    while not (token == 'END' or token == 'ENDBLOCK') and not self.eof and not token==None:
                        token = self.read_next_token_ucase()
                        if token == 'DIMENSIONS':
                            self.parse_dimensions_statement()     
                        if token == 'FORMAT':
                            self.parse_format_statement()                               
                        if token == 'MATRIX':
                            self.parse_matrix_statement()
                    self.skip_to_semicolon() # move past END command
                elif token == 'TREES':
                    trees_block = self.trees_block_factory()
                    trees_block.taxa_block = self.get_default_taxa_block()
                    self.dataset.add_trees_block(trees_block=trees_block)
                    self.skip_to_semicolon() # move past BEGIN command
                    while not (token == 'END' or token == 'ENDBLOCK') and not self.eof and not token==None:
                        token = self.read_next_token_ucase()
                        if token == 'TRANSLATE':
                            self.parse_translate_statement()                         
                        if token == 'TREE':
                            self.parse_tree_statement(trees_block)    
                    self.skip_to_semicolon() # move past END command    
                else:
                    # unknown block
                    while not (token == 'END' or token == 'ENDBLOCK') and not self.eof and not token==None:
                        #print token
                        self.skip_to_semicolon()
                        token = self.read_next_token_ucase()
        return self.dataset
    
    def _get_current_file_char(self):
        """
        Returns the current character from the file stream.
        """
        if self.__current_file_char == None:
            self.__current_file_char = self.read_next_char()
        return self.__current_file_char

    def _set_current_file_char(self, new_char):
        self.__current_file_char = new_char

    current_file_char = property(_get_current_file_char, _set_current_file_char)
    
    def read_next_char(self):
        """
        Advances the file stream cursor to the next character and returns 
        it.
        """
        if self.filehandle:
            read_char = self.filehandle.read(1) # returns empty string if EOF
            if read_char == '':
                self.eof = True
            else:
                if self.previous_file_char == '\n':
                    self.current_line_number = self.current_line_number + 1
                    self.current_col_number = 0
                self.previous_file_char = self.__current_file_char
                self.current_col_number = self.current_col_number + 1
            self.current_file_char = read_char
            return self.__current_file_char
        else:
            return None    
            
    def skip_comment(self):
        """
        Reads characters from the file until current comment block (and 
        any nested comment block terminates. Assumes current cursor 
        position is on first character inside comment block.
        """
        while self.current_file_char != ']' and not self.eof:
            if self.read_next_char() == '[':
                self.skip_comment()
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
        """
        Advances to the first non-whitespace character outside a comment block.
        """
        while (NexusReader.is_whitespace(self.current_file_char) or self.current_file_char=='[') and not self.eof:
            if self.current_file_char=='[':
                self.read_noncomment_character()
            else:
                self.read_next_char()
        return self.current_file_char
        
    def read_next_token(self, ignore_punctuation=None, preserve_quotes=False):
        """
        Reads the next token in the file stream. A token in this context is any word or punctuation character
        outside of a comment block.
        """
        if ignore_punctuation == None:
            ignore_punctuation = []        
        if not self.eof:
            token = ''
            self.skip_to_significant_character()
            if not self.eof:
                if self.current_file_char == "'":
                    self.read_next_char()
                    end_quote = False
                    while not end_quote and not self.eof:  
                        if self.current_file_char == "'":                        
                            self.read_next_char()
                            if self.current_file_char == "'":
                                token = token + "'"
                                self.read_next_char()
                            else:
                                end_quote = True
                        else:
                            token = token + self.current_file_char
                            self.read_next_char()
                    if preserve_quotes:
                        token = "'" + token + "'"                   
                else:
                    if NexusReader.is_punctuation(self.current_file_char) and self.current_file_char not in ignore_punctuation:
                        token = self.current_file_char
                        self.read_next_char()
                    else:
                        while not self.eof and (not NexusReader.is_whitespace_or_punctuation(self.current_file_char) or self.current_file_char in ignore_punctuation): 
                            token = token + self.current_file_char
                            self.read_next_char()
                self.current_token = token
            else:
                self.current_token = None
        else:
            self.current_token = None
        return self.current_token

    def read_next_token_ucase(self, ignore_punctuation=[]):
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
        """
        Advances the file stream cursor to the next semi-colon.
        """
        token = self.read_next_token()
        while token != ';' and not self.eof and token != None: 
            token = self.read_next_token()
            pass
    
    def syntax_exception(self, message):
        """
        Returns an exception object parameterized with line and 
        column number values.
        """
        return NexusReader.SyntaxException(self.current_line_number, self.current_col_number, message)
    
    def parse_format_statement(self):
        """
        Processes a FORMAT command. Assumes that the file reader is 
        positioned right after the "FORMAT" token in a FORMAT command.
        """
        token = self.read_next_token_ucase()
        while token != ';':
            #print token
            if token == 'DATATYPE':
                token = self.read_next_token_ucase()
                if token == '=':
                    token = self.read_next_token_ucase()
                    if token == "DNA" or token == "NUCLEOTIDES":                        
                        self.char_block_type = characters.DnaCharactersBlock
                    elif token == "RNA":
                        self.char_block_type = characters.RnaCharactersBlock
                    elif token == "PROTEIN":
                        self.char_block_type = characters.ProteinCharactersBlock
                    elif token == "STANDARD":
                        self.char_block_type = characters.StandardCharactersBlock
                        self.symbols = "12"
                else:
                    raise self.syntax_exception('Expecting "=" after DATATYPE keyword')
            elif token == 'SYMBOLS':
                token = self.read_next_token_ucase()
                if token == '=':
                    token = self.read_next_token_ucase()
                    if token == '"':
                        self.symbols = ""
                        token = self.read_next_token_ucase()
                        while token != '"':
                            if token not in self.symbols:
                                self.symbols = self.symbols + token
                            token = self.read_next_token_ucase()
                    else:
                        raise self.syntax_exception('Expecting \'"\' before beginning SYMBOLS list')
                else:
                    raise self.syntax_exception('Expecting "=" after SYMBOLS keyword')                    
            elif token == 'GAP':
                token = self.read_next_token_ucase()
                if token == '=':
                    token = self.read_next_token_ucase()
                    self.gap_char = token
                else:
                    raise self.syntax_exception('Expecting "=" after GAP keyword')
            elif token == 'MISSING':
                token = self.read_next_token_ucase()
                if token == '=':
                    token = self.read_next_token_ucase()
                    self.missing_char = token
                else:
                    raise self.syntax_exception('Expecting "=" after MISSING keyword')       
            elif token == 'MATCHCHAR':
                token = self.read_next_token_ucase()
                if token == '=':
                    token = self.read_next_token_ucase()
                    self.match_char = token
                else:
                    raise self.syntax_exception('Expecting "=" after MISSING keyword')                 
            token = self.read_next_token_ucase()   

    def parse_tree_statement(self, trees_block):
        """
        Processes a TREE command. Assumes that the file reader is 
        positioned right after the "TREE" token in a TREE command.
        Calls on the NewickStatementParser of the trees module.
        """                
        token = self.read_next_token()
        if token == '*':
            token = self.read_next_token()        
        tree_name = token
        token = self.read_next_token()
        if token != '=':
            raise self.syntax_exception('Expecting "=" in definition of Tree "%s" but found "%s"' % (tree_name, token))
        else:
            # collect entire tree statement by accumulating tokens until we reach a semi-colon
            statement = []
            token = self.read_next_token(preserve_quotes=True)
            while token and token != ';': 
                statement.append(token)
                token = self.read_next_token(preserve_quotes=True)
            newick_parser = newick.NewickTreeParser()
            tree = newick_parser.parse_tree_statement(tree_statement=''.join(statement), 
                                                      trees_block=trees_block, 
                                                      translate_dict=self.tree_translate_dict)
            tree.label = tree_name 
        if self.current_token != ';':
            self.skip_to_semicolon()        
        return tree
            
    def parse_translate_statement(self):
        """
        Processes a TRANSLATE command. Assumes that the file reader is 
        positioned right after the "TRANSLATE" token in a TRANSLATE command.
        """     
        token = self.current_token
        while token and token != ';':
            translation_token = self.read_next_token()
            translation_label = self.read_next_token()
            token = self.read_next_token() # ","
            #print translation_token, translation_label, token
            if token != ',' and token != ';':
                raise self.syntax_exception('Expecting "," in TRANSLATE statement after definition for %s = "%s", but found "%s" instead' % (translation_token, translation_label, token))
            else:
                self.tree_translate_dict[translation_token] = translation_label

    def parse_dimensions_statement(self):
        """
        Processes a DIMENSIONS command. Assumes that the file reader is 
        positioned right after the "DIMENSIONS" token in a DIMENSIONS command.
        """
        token = self.read_next_token_ucase()
        while token != ';':
            #print token
            if token == 'NTAX':
                token = self.read_next_token_ucase()
                if token == '=':
                    token = self.read_next_token_ucase()
                    if token.isdigit():
                        self.file_specified_ntax = int(token)
                    else:
                        raise self.syntax_exception('Expecting numeric value for NTAX')
                else:
                    raise self.syntax_exception('Expecting "=" after NTAX keyword')
            elif token == 'NCHAR':
                token = self.read_next_token_ucase()
                if token == '=':
                    token = self.read_next_token_ucase()
                    if token.isdigit():
                        self.file_specified_nchar = int(token)
                    else:
                        raise self.syntax_exception('Expecting numeric value for NCHAR')
                else:
                    raise self.syntax_exception('Expecting "=" after NCHAR keyword')
            token = self.read_next_token_ucase()
    
    def get_default_taxa_block(self, taxa_block=None):
        if taxa_block is None:
            if len(self.dataset.taxa_blocks) == 0:
                taxa_block = self.taxa_block_factory()
                self.dataset.add_taxa_block(taxa_block=taxa_block)
            else:
                taxa_block = self.dataset.taxa_blocks[0]
        else:
            if taxa_block not in self.dataset.taxa_blocks:
                self.dataset.taxa_blocks.insert(taxa_block)
        return taxa_block
         
    def parse_taxlabels_statement(self, taxa_block=None):
        """
        Processes a TAXLABELS command. Assumes that the file reader is 
        positioned right after the "TAXLABELS" token in a TAXLABELS command.
        """
        taxa_block = self.get_default_taxa_block(taxa_block)
        token = self.read_next_token()
        while token != ';':
            taxon = taxa.Taxon(label=token)
            taxa_block.append(taxon)
            token = self.read_next_token()            
            
    def build_state_alphabet(self, char_block, symbols):
        sa = characters.StateAlphabet()
        for symbol in symbols:
            sa.append(characters.StateAlphabetElement(symbol=symbol))        
        if self.missing_char:            
            sa.append(characters.StateAlphabetElement(symbol=self.missing_char,
                                           multistate=characters.StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=sa.get_states(symbols=symbols)))            
        char_block.state_alphabets = [sa]
        char_block.default_state_alphabet = char_block.state_alphabets[0]
        return char_block
    
    def parse_matrix_statement(self, taxa_block=None):
        """
        Processes a MATRIX command. Assumes that the file reader 
        is positioned right after the "MATRIX" token in a MATRIX command, 
        and that NTAX and NCHAR have been specified accurately.
        """
        taxa_block = self.get_default_taxa_block(taxa_block)
        if not self.file_specified_ntax:
            raise self.syntax_exception('NTAX must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        elif not self.file_specified_nchar:
            raise self.syntax_exception('NCHAR must be defined by DIMENSIONS command to non-zero value before MATRIX command')
        else:
            char_block = self.char_block_type()
            char_block.taxa_block = taxa_block
            if isinstance(char_block, characters.StandardCharactersBlock):
                char_block = self.build_state_alphabet(char_block, self.symbols)
            self.dataset.add_char_block(char_block=char_block)           
            symbol_state_map = char_block.default_state_alphabet.symbol_state_map()
            if True: # future: trap and handle no labels, transpose etc.
                token = self.read_next_token()
                while token != ';' and not self.eof:
                    taxon = taxa_block.find_taxon(label=token, update=True)
                    if taxon not in char_block:
                        char_block[taxon] = characters.CharacterDataVector(taxon=taxon)
                    if self.interleave:
                        while self.current_file_char != '\n' and self.current_file_char != '\r':
                            if self.current_file_char not in [' ', '\t']:
                                state = symbol_state_map[self.current_file_char]
                                char_block[taxon].append(characters.CharacterDataCell(value=state))
                            self.read_next_char()
                    else:
                        while len(char_block[taxon]) < self.file_specified_nchar and not self.eof:
                            char_group = self.read_next_token()
                            char_group = parse_sequence_iupac_ambiguities(char_group)
                            for char in char_group:
                                state = symbol_state_map[char]
                                char_block[taxon].append(characters.CharacterDataCell(value=state))                   
                    token = self.read_next_token()
            else:
                ## TODO: NO LABELS/TRANSPOSED ##
                pass

class NexusWriter(datasets.Writer):
    """
    Implements the DataWriter interface for handling NEXML files.
    """

    def __init__(self, simple=False):
        """
        Calls the base class constructor.
        """
        datasets.Writer.__init__(self)
        self.simple = simple
        self.comment = []
        
    def compose_taxlabel(label):
        if NexusReader.has_whitespace(label) or NexusReader.has_punctuation(label):
            return "'" + label + "'"
        else:
            return label
    compose_taxlabel = staticmethod(compose_taxlabel)

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
        newick_writer = newick.NewickTreeWriter()
        block.append('BEGIN TREES;')
        for treeidx, tree in enumerate(trees_block):
            if tree.label:
                tree_name = tree.label
            else:
                tree_name = str(treeidx)
            newick_str = newick_writer.compose_node(tree.seed_node)
            block.append('    tree %s = %s;' % (tree_name, newick_str))
        block.append('END;\n\n')
        dest.write('\n'.join(block))
        
    def write_taxa_block(self, taxa_block, dest):
        block = []
        block.append('BEGIN TAXA;')
        block.append('    DIMENSIONS NTAX=%d;' % len(taxa_block))
        block.append('    TAXLABELS')
        for taxon in taxa_block:
            block.append('        %s' % self.compose_taxlabel(taxon.label))
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
        taxlabels = [taxon.label for taxon in char_block]
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
#        nexus.append('    FORMAT DATATYPE=%s GAP=%s MISSING=%s MATCHCHAR=%s;' % (self.char_block_type, self.gap_char, self.missing_char, self.match_char))
        nexus.append('    MATRIX')
        taxa = char_block.keys()
        taxa.sort()
        for taxon in taxa:
            seq_vec = char_block[taxon]
            #print t
            seq = ''.join([str(seq_sym) for seq_sym in seq_vec])
#             seq.replace('~','-')
            seq = seq.ljust(nchar, '-')
            nexus.append('%s    %s' % (self.compose_taxlabel(taxon.label).ljust(max_label_len), seq))            
        nexus.append('    ;')
        nexus.append('END;\n\n')
        dest.write('\n'.join(nexus))        



def _io_test(source):
    nexus_reader = NexusReader()
    dataset = nexus_reader.get_dataset(filepath=source)
    from dendropy import nexml
    nexmlw = nexml.NexmlWriter()
    print nexmlw.compose_dataset(dataset)
    nexus_writer = NexusWriter()
    print nexus_writer.compose_dataset(dataset)
    
if __name__ == "__main__":
    source1 = "/Volumes/KANSAS/jeet/Documents/Codeworks/Portfolios/Python/Projects/Phylogenetics/DendroPy/versions/trunk/dendropy/tests/files/primates.nex"
    source2 = "/Volumes/KANSAS/jeet/Documents/Codeworks/Portfolios/Python/Projects/Phylogenetics/DendroPy/versions/trunk/dendropy/tests/files/nexus_trees.tre"
    #source = "/home/jeet/Documents/Codeworks/Portfolios/Python/Projects/Phylogenetics/DendroPy/versions/trunk/dendropy/tests/files/primate-mtDNA-interleaved.nex"        
    _io_test(source1)
