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
from dendropy.utility.error import DataParseError
from dendropy import dataobject
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

###############################################################################
## RootingInterpreter

class RootingInterpreter(object):

    def evaluate_as_rooted_kwargs(kwdict, default):
        if "as_rooted" in kwdict and "as_unrooted" in kwdict \
                and (kwdict["as_rooted"] != (not kwdict["as_unrooted"])):
            raise TypeError("Conflict rooting specification: 'as_rooted'=%s and 'as_unrooted'=%s" \
                    % (kwdict["as_rooted"], kwdict["as_unrooted"]))
        if "as_rooted" in kwdict:
            if kwdict["as_rooted"] is True or kwdict["as_rooted"] is False:
                return kwdict["as_rooted"]
            else:
                raise ValueError("Invalid value for 'as_rooted' (expecting True/False, but received '%s')" \
                        % kwdict["as_rooted"])
        elif "as_unrooted" in kwdict:
            if kwdict["as_unrooted"] is True or kwdict["as_unrooted"] is False:
                return not kwdict["as_unrooted"]
            else:
                raise ValueError("Invalid value for 'as_unrooted' (expecting True/False, but received '%s')" \
                        % kwdict["as_rooted"])
        else:
            return default

    evaluate_as_rooted_kwargs = staticmethod(evaluate_as_rooted_kwargs)

    def evaluate_default_as_rooted_kwargs(kwdict, default=None):
        if "default_as_rooted" in kwdict and "default_as_unrooted" in kwdict \
                and (kwdict["default_as_rooted"] != (not kwdict["default_as_unrooted"])):
            raise TypeError("Conflict rooting specification: 'default_as_rooted'=%s and 'default_as_unrooted'=%s" \
                    % (kwdict["default_as_rooted"], kwdict["default_as_unrooted"]))
        if "default_as_rooted" in kwdict:
            if kwdict["default_as_rooted"] is True or kwdict["default_as_rooted"] is False:
                return kwdict["default_as_rooted"]
            else:
                raise ValueError("Invalid value for 'default_as_rooted' (expecting True/False, but received '%s')" \
                        % kwdict["default_as_rooted"])
        elif "default_as_unrooted" in kwdict:
            if kwdict["default_as_unrooted"] is True or kwdict["default_as_unrooted"] is False:
                return not kwdict["default_as_unrooted"]
            else:
                raise ValueError("Invalid value for 'default_as_unrooted' (expecting True/False, but received '%s')" \
                        % kwdict["default_as_rooted"])
        else:
            return default

    evaluate_default_as_rooted_kwargs = staticmethod(evaluate_default_as_rooted_kwargs)

    def __init__(self, **kwargs):
        self._as_rooted = RootingInterpreter.evaluate_as_rooted_kwargs(kwargs, None)
        self._default_as_rooted = RootingInterpreter.evaluate_default_as_rooted_kwargs(kwargs, False)

    def update(self, **kwargs):
        self._as_rooted = RootingInterpreter.evaluate_as_rooted_kwargs(kwargs, self._as_rooted)
        self._default_as_rooted = RootingInterpreter.evaluate_default_as_rooted_kwargs(kwargs, self._default_as_rooted)

    def _get_as_rooted(self):
        return self._as_rooted

    def _set_as_rooted(self, val):
        self._as_rooted = val

    as_rooted = property(_get_as_rooted, _set_as_rooted)

    def _get_as_unrooted(self):
        return not self._as_rooted

    def _set_as_unrooted(self, val):
        self._as_rooted = not val

    as_unrooted = property(_get_as_unrooted, _set_as_unrooted)

    def _get_default_as_rooted(self):
        return self._default_as_rooted

    def _set_default_as_rooted(self, val):
        self._default_as_rooted = val

    default_as_rooted = property(_get_default_as_rooted, _set_default_as_rooted)

    def _get_default_as_unrooted(self):
        return not self._default_as_rooted

    def _set_default_as_unrooted(self, val):
        self._default_as_rooted = not val

    default_as_unrooted = property(_get_default_as_unrooted, _set_default_as_unrooted)

    def interpret_as_rooted(self, tree_rooting_comment=None, **kwargs):
        if self.as_rooted is not None:
            return self.as_rooted
        elif tree_rooting_comment is not None:
            return tree_rooting_comment.upper() == "&R"
        else:
            return self.default_as_rooted

    def interpret_as_unrooted(self, **kwargs):
        return not self.interpret_as_rooted(**kwargs)

###############################################################################
## StrToTaxon

class StrToTaxon(object):

    class MultipleTaxonUseError(DataParseError):
        def __init__(self, *args, **kwargs):
            DataParseError.__init__(self, *args, **kwargs)

    def __init__(self, taxon_set, translate_dict=None, allow_repeated_use=False):
        """If `allow_repeated_use` is True, then get_taxon and require_taxon
        can be called multiple times to get the same taxon.  If it is false then
        calling the functions with the same label will generate a DataParseError
        indicating that the taxon has been used multiple times."""
        self.taxon_set = taxon_set
        self.translate = translate_dict or {}
        if allow_repeated_use:
            self.returned_taxon_set = None
        else:
            self.returned_taxon_set = set()

    def _returning(self, t, label):
        #_LOG.debug("Checking if we can return %s" % str(t))
        if (self.returned_taxon_set is not None) and (t is not None):
            if t in self.returned_taxon_set:
                raise StrToTaxon.MultipleTaxonUseError(message="Taxon %s used twice (it appears as %s the second time)" % (str(t), label))
            self.returned_taxon_set.add(t)
        return t

    def get_taxon(self, label):
        t = self.translate.get(label)
        if t is None:
            t = self.taxon_set.get_taxon(label=label)
        return self._returning(t, label)

    def require_taxon(self, label):
        v = self.get_taxon(label)
        if v is not None:
            return v
        t = self.taxon_set.require_taxon(label=label)
        return self._returning(t, label)
#        if t is not None:
#            self.translate[label] = t #@this could lead to problems when we support multiple taxon blocks, but now it'll speed thing up

    def index(self, t):
        return self.taxon_set.index(t)

###############################################################################
## parse_tree_from_stream

def parse_tree_from_stream(stream_tokenizer, **kwargs):
    """
    Processes a (SINGLE) TREE statement. Assumes that the input stream is
    located at the beginning of the statement (i.e., the first non-comment
    token should be the opening parenthesis of the tree definition).

    str_to_taxon kwarg (if used) must supply the StrToTaxon interface).
    """
    translate_dict = kwargs.get("translate_dict", None)
    encode_splits = kwargs.get("encode_splits", False)
    rooting_interpreter = kwargs.get("rooting_interpreter", RootingInterpreter(**kwargs))
    finish_node_func = kwargs.get("finish_node_func", None)
    edge_len_type = kwargs.get("edge_len_type", float)
    taxon_set = kwargs.get("taxon_set", None)
    suppress_internal_node_taxa = kwargs.get("suppress_internal_node_taxa", False)
    if taxon_set is None:
        taxon_set = dataobject.TaxonSet()
    tree = dataobject.Tree(taxon_set=taxon_set)

    stream_tokenizer.tree_rooting_comment = None # clear previous comment
    token = stream_tokenizer.read_next_token()
    if not token:
        return None
    tree.is_rooted = rooting_interpreter.interpret_as_rooted(stream_tokenizer.tree_rooting_comment)
    if encode_splits:
        if len(taxon_set) == 0:
            raise Exception("When encoding splits on a tree as it is being parsed, a "
                + "fully pre-populated TaxonSet object must be specified using the 'taxon_set' keyword " \
                + "to avoid taxon/split bitmask values changing as new Taxon objects are created " \
                + "and added to the TaxonSet.")
        if tree.is_rooted:
            tree.split_edges = {}
        else:
            atb = taxon_set.all_taxa_bitmask()
            d = containers.NormalizedBitmaskDict(mask=atb)
            tree.split_edges = d
        split_map = tree.split_edges

    stt = kwargs.get('str_to_taxon')
    if stt is None:
        stt = StrToTaxon(taxon_set, translate_dict, allow_repeated_use=False)

    tree.seed_node = dataobject.Node()
    curr_node = tree.seed_node
    if encode_splits:
        curr_node.edge.split_bitmask = 0L

    ### NHX format support ###
    def store_node_comments(active_node):
        if stream_tokenizer.comments:
            active_node.comments.extend(stream_tokenizer.comments)

    stream_tokenizer.clear_comments()
    while True:
        if not token or token == ';':
            if curr_node is not tree.seed_node:
                raise stream_tokenizer.data_format_error("Unbalanced parentheses -- not enough ')' characters found in tree description")
            if encode_splits:
                split_map[curr_node.edge.split_bitmask] = curr_node.edge
            break
        if token == '(':
            if not curr_node.parent_node:
                if curr_node.child_nodes():
                    raise stream_tokenizer.data_format_error("Unexpected '(' after the tree description.  Expecting a label for the root or a ;")
            tmp_node = dataobject.Node()
            if encode_splits:
                tmp_node.edge.split_bitmask = 0L
            curr_node.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
            store_node_comments(curr_node)
        elif token == ',':
            tmp_node = dataobject.Node()
            if curr_node.is_leaf() and not curr_node.taxon:
#                 curr_node.taxon = taxon_set.Taxon(oid="UNAMED_" + str(id(curr_node)), label='')
#                 taxon_set.add(curr_node.taxon)
                raise stream_tokenizer.data_format_error("Missing taxon specifier in a tree -- found either a '(,' or ',,' construct.")
            p = curr_node.parent_node
            if not p:
                raise stream_tokenizer.data_format_error("Comma found one the 'outside' of a newick tree description")
            if encode_splits:
                tmp_node.edge.split_bitmask = 0L
                e = curr_node.edge
                u = e.split_bitmask
                split_map[u] = e
                p.edge.split_bitmask |= u
            if finish_node_func is not None:
                finish_node_func(curr_node, tree)
            p.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
            store_node_comments(curr_node)
        else:
            if token == ')':
                if curr_node.is_leaf() and not curr_node.taxon:
                    raise stream_tokenizer.data_format_error("Missing taxon specifier in a tree -- found either a '(,' or ',,' construct.")
                p = curr_node.parent_node
                if not p:
                    raise stream_tokenizer.data_format_error("Unbalanced parentheses -- too many ')' characters found in tree description")
                if encode_splits:
                    e = curr_node.edge
                    u = e.split_bitmask
                    p.edge.split_bitmask |= u
                    split_map[u] = curr_node.edge
                if finish_node_func is not None:
                    finish_node_func(curr_node, tree)
                curr_node = p
            else:
                is_leaf = curr_node.is_leaf()
                if is_leaf:
                    if curr_node.taxon:
                        raise stream_tokenizer.data_format_error("Multiple labels found for the same leaf (taxon '%s' and label '%s')" % (str(curr_node.taxon), token))
                    try:
                        t = stt.require_taxon(label=token)
                    except StrToTaxon.MultipleTaxonUseError, e:
                        raise stream_tokenizer.data_format_error(e.msg)
                else:
                    if curr_node.label:
                        raise stream_tokenizer.data_format_error("Multiple labels found for the same leaf (taxon '%s' and label '%s')" % (curr_node.label, token))
                    if suppress_internal_node_taxa:
                        t = None
                    else:
                        try:
                            t = stt.get_taxon(label=token)
                        except StrToTaxon.MultipleTaxonUseError, e:
                            raise stream_tokenizer.data_format_error(e.msg)
                if t is None:
                    curr_node.label = token
                else:
                    curr_node.taxon = t
                    if encode_splits:
                        try:
                            cm = t.split_bitmask
                        except:
                            cm = 1 << (stt.index(t))
                        e = curr_node.edge
                        e.split_bitmask = cm
                        split_map[cm] = e

            token = stream_tokenizer.read_next_token()
            store_node_comments(curr_node)
            if token == ':':
                edge_length_str = stream_tokenizer.read_next_token(ignore_punctuation='-+.')
                store_node_comments(curr_node)
                if not edge_length_str:
                    raise stream_tokenizer.data_format_error("Expecting a branch length after : but encountered the end of the tree description" )
                try:
                    curr_node.edge.length = edge_len_type(edge_length_str)
                except:
                    curr_node.edge.length = edge_length_str
                token = stream_tokenizer.read_next_token()
                store_node_comments(curr_node)
    return tree


###############################################################################
## NexusTokenizer

class NexusTokenizer(object):
    "Encapsulates reading NEXUS/NEWICK tokens from file."

    #######################################################################
    ## FOR COMMUNICATING PARSER POSITION/STATUS

    class BlockTerminatedException(Exception):
        pass

    class EolException(Exception):
        pass

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

    def __init__(self, stream_handle=None, **kwargs):
        self._reset()
        self.preserve_underscores = kwargs.get('preserve_underscores', False)
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
        self.comments = []
        self.tree_rooting_comment = None
        self.last_comment_parsed = None
        self.preserve_underscores = False

    def _get_current_file_char(self):
        "Returns the current character from the file stream."
        if self._current_file_char == None:
            self._current_file_char = self.read_next_char()
        return self._current_file_char

    def _set_current_file_char(self, new_char):
        self._current_file_char = new_char

    current_file_char = property(_get_current_file_char, _set_current_file_char)

    def clear_comments(self):
        self.comments = []

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
        comment = cmt_body.getvalue()
        self.comments.append(comment)
        self.last_comment_parsed = comment
        if comment.strip().upper() == "&R":
            self.tree_rooting_comment = "&R"
        elif comment.strip().upper() == "&U":
            self.tree_rooting_comment = "&U"
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

    def read_statement_tokens_till_eol(self, token_list=None, ignore_punctuation=None):
        """
        Reads all tokens in the file stream until EOL, EOF, or ';'
        is encountered.
        """
        _is_newline = lambda x: x == '\n' or x == '\r'
        self.comments = []
        if ignore_punctuation == None:
            ignore_punctuation = []
        self.current_token = None
        if self.eof:
            return None

        # skip to first significant
        while (NexusTokenizer.is_whitespace(self.current_file_char) \
                or self.current_file_char=='[') \
                and not self.eof:
            if self.current_file_char == '[':
                self.read_noncomment_character()
            elif _is_newline(self.current_file_char):
                return None
            else:
                self.read_next_char()
        if self.eof:
            self.data_format_error("Unexpected end of file before statement completion")
        elif self.current_file_char == ';':
            return None

        c = self.current_file_char
        if token_list is None:
            token_list = []
        try:
            while not self.eof:
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
                        raise self.data_format_error("Unexpected end of file inside quoted token")
                    token_list.append(token.getvalue())
                else:
                    token = StringIO()
                    quick_check = NexusTokenizer.is_whitespace_or_punctuation
                    while True:
                        if quick_check(c):
                            if c == '[':
                                self.skip_comment()
                                c = self.current_file_char
                                continue
                            if quick_check(c) and not c in ignore_punctuation:
                                if not NexusTokenizer.is_whitespace(c):
                                    token.write(c)
                                break
                        if c == '_' and not self.preserve_underscores:
                            c = ' '
                        token.write(c)
                        c = self.read_next_char()
                        if not c:
                            break
                    token = token.getvalue()
                    if token == ";":
                        raise self.BlockTerminatedException()
                    elif token:
                        token_list.append(token)
                    if c == "\r" or c == "\n":
                        raise self.EolException()
                    else:
                        c = self.read_next_char()
        except self.EolException:
            pass
        except self.BlockTerminatedException:
            raise
        return token_list

    def skip_to_significant_character(self):
        "Advances to the first non-whitespace character outside a comment block."
        while (NexusTokenizer.is_whitespace(self.current_file_char) \
                or self.current_file_char=='[') and not self.eof:
            if self.current_file_char=='[':
                self.read_noncomment_character()
            else:
                self.read_next_char()
        return self.current_file_char

    def read_next_token(self, ignore_punctuation=None):
        """
        Reads the next token in the file stream. A token in this context
        is any word or punctuation character outside of a comment block.
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
                raise self.data_format_error("Unexpected end of file inside quoted token")
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
            fget = self.read_next_char
            quick_check = NexusTokenizer.is_whitespace_or_punctuation
            if ignore_punctuation:
                while True:
                    if quick_check(c):
                        if c == '[':
                            self.skip_comment()
                            c = self.current_file_char
                            continue
                        if quick_check(c) and not c in ignore_punctuation:
                            break
                    if c == '_' and not self.preserve_underscores:
                        c = ' '
                    token.write(c)
                    c = fget()
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
                    if c == '_' and not self.preserve_underscores:
                        c = ' '
                    token.write(c)
                    c = fget()
                    if not c:
                        break
            tokenstr = token.getvalue()
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

    def data_format_error(self, message):
            """
            Returns an exception object parameterized with line and
            column number values.
            """
            return DataParseError(message=message,
                                   row=self.current_line_number,
                                   column=self.current_col_number)
