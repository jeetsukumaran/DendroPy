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
Functions and utilities for handling NEXUS- and NEWICK-formatted data, used
by both `dendropy.newick` and `dendropy.nexus` modules.
"""

import re
from cStringIO import StringIO
from dendropy.utility import containers
from dendropy.utility.errors import DataFormatError
from dendropy import dataobject

###############################################################################
## RootingInterpretation

class RootingInterpretation:
    UNKNOWN_DEF_ROOTED = 0
    UNKNOWN_DEF_UNROOTED = 1
    ROOTED = 2
    UNROOTED = 3

###############################################################################
## StrToTaxon

class StrToTaxon(object):

    def __init__(self, taxon_set, translate_dict=None):
        self.taxon_set = taxon_set
        self.translate = translate_dict or {}

    def get_taxon(self, label):
        return self.translate.get(label)

    def require_taxon(self, label):
        v = self.get_taxon(label)
        if v is not None:
            return v
        t = self.taxon_set.require_taxon(label=label)
        if t is not None:
            self.translate[label] = t #@this could lead to problems when we support multiple taxon blocks, but now it'll speed thing up
        return t

    def index(self, t):
        return self.taxon_set.index(t)

###############################################################################
## parse_tree_from_stream

def parse_tree_from_stream(stream_tokenizer, **kwargs):
    """
    Processes a (SINGLE) TREE statement. Assumes that the input stream is
    located at the beginning of the statement (i.e., the first
    parenthesis that defines the tree).
    """

    translate_dict = kwargs.get("translate_dict", None)
    encode_splits = kwargs.get("encode_splits", False)
    rooted = kwargs.get("rooted", RootingInterpretation.UNKNOWN_DEF_UNROOTED)
    finish_node_func = kwargs.get("finish_node_func", None)
    edge_len_type = kwargs.get("edge_len_type", float)
    taxon_set = kwargs.get("taxon_set", None)
    token = stream_tokenizer.read_next_token()
    if not token:
        return None
    if taxon_set is None:
        taxon_set = dataobject.TaxonSet()
    tree = dataobject.Tree(taxon_set=taxon_set)

    if encode_splits:
        if len(taxon_set) == 0:
            raise Exception("When encoding splits on a tree as it is being parsed, a "
                + "fully pre-populated TaxonSet object must be specified using the 'taxon_set' keyword " \
                + "to avoid taxon/split bitmask values changing as new Taxon objects are created " \
                + "and added to the TaxonSet.")
        if rooted == RootingInterpretation.UNKNOWN_DEF_ROOTED \
                or rooted == RootingInterpretation.UNKNOWN_DEF_UNROOTED:
            r = stream_tokenizer.tree_rooted_comment()
            if r is not None:
                if r:
                    rooted = RootingInterpretation.ROOTED
                else:
                    rooted = RootingInterpretation.UNROOTED

        if rooted == RootingInterpretation.UNKNOWN_DEF_ROOTED \
                or rooted == RootingInterpretation.ROOTED:
            tree.split_edges = {}
        else:
            atb = taxon_set.all_taxa_bitmask()
            d = containers.NormalizedBitmaskDict(mask=atb)
            tree.split_edges = d
        split_map = tree.split_edges

    stt = StrToTaxon(taxon_set, translate_dict)

    tree.seed_node = dataobject.Node()
    curr_node = tree.seed_node
    if encode_splits:
        curr_node.edge.clade_mask = 0L

    while True:
        if not token or token == ';':
            if curr_node is not tree.seed_node:
                raise stream_tokenizer.syntax_exception("Unbalanced parentheses -- not enough ')' characters found in tree description")
            if encode_splits:
                split_map[curr_node.edge.clade_mask] = curr_node.edge
            break
        if token == '(':
            tmp_node = dataobject.Node()
            if encode_splits:
                tmp_node.edge.clade_mask = 0L
            curr_node.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
        elif token == ',':
            tmp_node = dataobject.Node()
            if curr_node.is_leaf() and not curr_node.taxon:
#                 curr_node.taxon = taxon_set.Taxon(oid="UNAMED_" + str(id(curr_node)), label='')
#                 taxon_set.add(curr_node.taxon)
                raise stream_tokenizer.syntax_exception("Missing taxon specifier in a tree -- found either a '(,' or ',,' construct.")
            p = curr_node.parent_node
            if not p:
                raise stream_tokenizer.syntax_exception("Comma found one the 'outside' of a newick tree description")
            if encode_splits:
                tmp_node.edge.clade_mask = 0L
                e = curr_node.edge
                u = e.clade_mask
                split_map[u] = e
                p.edge.clade_mask |= u
            if finish_node_func is not None:
                finish_node_func(curr_node, tree)
            p.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
        else:
            if token == ')':
                if curr_node.is_leaf() and not curr_node.taxon:
#                     curr_node.taxon = dataobject.Taxon(oid="UNAMED_" + str(id(curr_node)), label='')
#                     taxon_set.add(curr_node.taxon)
                    raise stream_tokenizer.syntax_exception("Missing taxon specifier in a tree -- found either a '(,' or ',,' construct.")
                p = curr_node.parent_node
                if not p:
                    raise stream_tokenizer.syntax_exception("Unbalanced parentheses -- too many ')' characters found in tree description")
                if encode_splits:
                    e = curr_node.edge
                    u = e.clade_mask
                    p.edge.clade_mask |= u
                    split_map[u] = curr_node.edge
                if finish_node_func is not None:
                    finish_node_func(curr_node, tree)
                curr_node = p
            else:
                is_leaf = curr_node.is_leaf()
                if is_leaf:
                    t = stt.require_taxon(label=token)
                else:
                    t = stt.get_taxon(label=token)
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
                edge_length_str = stream_tokenizer.read_next_token(ignore_punctuation='-+.')
                try:
                    curr_node.edge.length = edge_len_type(edge_length_str)
                except:
#                     print "FAILED TO COERCE '%s'" % edge_length_str
                    curr_node.edge.length = edge_length_str
                token = stream_tokenizer.read_next_token()
    return tree


###############################################################################
## NexusTokenizer

class NexusTokenizer(object):
    "Encapsulates reading NEXUS/NEWICK tokens from file."

    #######################################################################
    ## STATIC METHODS

    punctuation = '\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>'
    whitespace = ' \0\t\n\r'
    multi_char_needs_quoting = re.compile(' \0\t\n\r\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>_')
    single_char_needs_quoting = re.compile(' \0\t\n\r\[\'_')
    def is_whitespace(char):
        return char in NexusTokenizer.whitespace

    is_whitespace = staticmethod(is_whitespace)

    def has_whitespace(s):
        return re.search('['+NexusTokenizer.whitespace+']', s) != None

    has_whitespace = staticmethod(has_whitespace)

    def is_punctuation(char):
        return char in NexusTokenizer.punctuation

    is_punctuation = staticmethod(is_punctuation)

    def has_punctuation(s):
        return re.search('['+NexusTokenizer.punctuation+']', s) != None

    has_punctuation = staticmethod(has_punctuation)

    def is_whitespace_or_punctuation(char):
        return NexusTokenizer.is_whitespace(char) or NexusTokenizer.is_punctuation(char)

    is_whitespace_or_punctuation = staticmethod(is_whitespace_or_punctuation)

    def has_whitespace_or_punctuation(s):
        return NexusTokenizer.has_whitespace(s) or NexusTokenizer.has_punctuation(s)

    has_whitespace_or_punctuation = staticmethod(has_whitespace_or_punctuation)

    def validate_identifier(label):
        if NexusTokenizer.has_whitespace_or_punctuation(label):
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
        if self._current_file_char == None:
            self._current_file_char = self.read_next_char()
        return self._current_file_char

    def _set_current_file_char(self, new_char):
        self._current_file_char = new_char

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
                self.previous_file_char = self._current_file_char
                self.current_col_number += 1
            self.current_file_char = read_char
            return self._current_file_char
        return None

    def _raw_read_next_char(self):
            read_char = self.stream_handle.read(1) # returns empty string if EOF
            if read_char == '':
                raise StopIteration()
            if self.previous_file_char == '\n':
                self.current_line_number = self.current_line_number + 1
                self.current_col_number = 0
            self.previous_file_char = self._current_file_char
            self.current_col_number += 1
            self.current_file_char = read_char
            return self._current_file_char

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
        while (NexusTokenizer.is_whitespace(self.current_file_char) or self.current_file_char=='[') and not self.eof:
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
            quick_check = NexusTokenizer.is_punctuation
            if quick_check(c) and c not in ignore_punctuation:
                if c == '_':
                    c = ' '
                self.current_token = c
                self.read_next_char()
                return c
            token = StringIO()
            fget = self.stream_handle.read
            quick_check = NexusTokenizer.is_whitespace_or_punctuation

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
            return DataFormatError(row=self.current_line_number,
                                   column=self.current_col_number,
                                   message=message)


