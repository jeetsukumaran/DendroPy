#! /usr/bin/env python

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

import sys
from dendropy.utility import error

##############################################################################
## Tokenizer

class Tokenizer(object):
    """
    Stream tokenizer.
    """

    class TokenizerError(error.DataParseError):

        def __init__(self,
                message=None,
                line_num=None,
                col_num=None,
                stream=None):
            error.DataParseError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class UnterminatedQuoteError(TokenizerError):

        def __init__(self,
                quote_char=None,
                line_num=None,
                col_num=None,
                stream=None):
            Tokenizer.TokenizerError.__init__(self,
                    message="Unterminated quote: {}".format(quote_char),
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    class UnexpectedEndOfStreamError(TokenizerError):

        def __init__(self,
                message=None,
                line_num=None,
                col_num=None,
                stream=None):
            Tokenizer.TokenizerError.__init__(self,
                    message=message,
                    line_num=line_num,
                    col_num=col_num,
                    stream=stream)

    def __init__(self,
            src,                        # source stream
            uncaptured_delimiters,      # delimiters between tokens (not returned)
            captured_delimiters,        # delimiters between tokens (returned as tokens)
            quote_chars,                # characters enclosing literals
            escape_quote_by_doubling,   # should two consecutive quote characters indicate a literal character (rather than a quote)?
            escape_chars,               # characters indicating beginning of escaped character
            comment_begin,              # string indicating beginning of comment
            comment_end,                # string indicating end of comment
            capture_comments,           # are comments to be stored?
            preserve_unquoted_underscores,       # are unquoted underscores to be preserved
            ):
        # Tokenizer behavior customization
        self.uncaptured_delimiters = uncaptured_delimiters
        self.captured_delimiters = captured_delimiters
        self.quote_chars = quote_chars
        self.escape_quote_by_doubling = escape_quote_by_doubling
        self.escape_chars = escape_chars
        self.comment_begin = comment_begin
        self.comment_end = comment_end
        self.capture_comments = capture_comments
        self.preserve_unquoted_underscores = preserve_unquoted_underscores

        # State (internals)
        self.src = src
        self._cur_char = None
        self.current_token = None
        self.is_token_quoted = False

        # Meta-information
        self.captured_comments = []
        self.current_line_num = 1
        self.current_column_num = 0
        self.token_line_num = 0
        self.token_column_num = 0

    def reset(self):
        self.set_stream(src=None)

    def set_stream(self, src=None):
        self.src = src
        self._cur_char = None
        self.current_token = None
        self.is_token_quoted = False
        self.captured_comments = []
        self.current_line_num = 1
        self.current_column_num = 0
        self.token_line_num = 0
        self.token_column_num = 0

    def is_eof(self):
        return self._cur_char == ""

    def has_captured_comments(self):
        return len(self.captured_comments) > 0

    def next_token(self):
        try:
            t = self.__next__()
            return t
        except StopIteration:
            self.current_token = None
            return None

    def require_next_token(self):
        try:
            t = self.__next__()
            return t
        except StopIteration:
            # In Python 3, if you catch an exception and then raise an
            # exception that is not a subclass of the original exception,
            # the original exception is not considered to have been
            # handled.
            # In Python > 3.3, this can be solved by:
            #
            #   raise Tokenizer.UnexpectedEndOfStreamError(
            #                   message="Unexpected end of stream",
            #                   line_num=self.current_line_num,
            #                   col_num=self.current_column_num,
            #                   stream=self.src) from None
            #
            # To accommodate other versions, the following
            # is required:
            exc = Tokenizer.UnexpectedEndOfStreamError(
                            message="Unexpected end of stream",
                            line_num=self.current_line_num,
                            col_num=self.current_column_num,
                            stream=self.src)
            exc.__context__ = None # Python 3.0, 3.1, 3.2
            exc.__cause__ = None # Python 3.3, 3.4
            raise exc

    def clear_captured_comments(self):
        del self.captured_comments[:]

    def pull_captured_comments(self):
        if not self.captured_comments:
            return None
        c = self.captured_comments[:]
        del self.captured_comments[:]
        return c

    def __iter__(self):
        return self

    def __next__(self):
        self.is_token_quoted = False
        if self._cur_char is None:
            self._get_next_char()
        self._skip_to_significant_char()
        if self._cur_char == "":
            raise StopIteration
        if self._cur_char in self.captured_delimiters:
            self.current_token = self._cur_char
            self.token_line_num = self.current_line_num
            self.token_column_num = self.current_column_num
            self._get_next_char()
            return self.current_token
        elif self._cur_char in self.quote_chars:
            self.token_line_num = self.current_line_num
            self.token_column_num = self.current_column_num
            dest = []
            self.is_token_quoted = True
            cur_quote_char = self._cur_char
            self._get_next_char()
            while True:
                if self._cur_char == "":
                    raise Tokenizer.UnterminatedQuoteError(
                            quote_char=cur_quote_char,
                            line_num=self.current_line_num,
                            col_num=self.current_column_num,
                            stream=src)
                if self._cur_char == cur_quote_char:
                    self._get_next_char()
                    if self.escape_quote_by_doubling:
                        if self._cur_char == cur_quote_char:
                            # dest.write(cur_quote_char)
                            dest.append(cur_quote_char)
                            self._get_next_char()
                        else:
                            break
                    else:
                        self._get_next_char()
                        break
                else:
                    # dest.write(self._cur_char)
                    dest.append(self._cur_char)
                    self._get_next_char()
            # self.current_token = dest.getvalue()
            self.current_token = "".join(dest)
            return self.current_token
        else:
            # unquoted
            self.token_line_num = self.current_line_num
            self.token_column_num = self.current_column_num
            dest = []
            self.is_token_quoted = False
            while self._cur_char != "":
                if self._cur_char in self.uncaptured_delimiters:
                    self._get_next_char()
                    break
                elif self._cur_char in self.captured_delimiters:
                    break
                elif self._cur_char in self.comment_begin:
                    self._handle_comment()
                    if self._cur_char == "":
                        break
                else:
                    if self._cur_char == "_" and not self.preserve_unquoted_underscores:
                        self._cur_char = " "
                    dest.append(self._cur_char)
                    self._get_next_char()
            # self.current_token = dest.getvalue()
            self.current_token = "".join(dest)
            if self.current_token == "":
                if self._cur_char != "":
                    self.__next__()
                else:
                    raise StopIteration
            return self.current_token
    next = __next__ # Python 2 legacy support

    def _skip_to_significant_char(self):
        if self._cur_char == "":
            return
        if self._cur_char is None:
            self._get_next_char()
        if self._cur_char not in self.uncaptured_delimiters:
            return
        while self._cur_char != "" and self._cur_char in self.uncaptured_delimiters:
            self._get_next_char()
        return

    def _get_next_char(self):
        self._cur_char = self.src.read(1)
        if self._cur_char != "":
            if self._cur_char == "\n":
                self.current_line_num += 1
                self.current_column_num = 1
            else:
                # print("@@@ {}: {}".format(self.current_column_num, self._cur_char))
                self.current_column_num += 1
        return self._cur_char

    def _handle_comment(self):
        dest = []
        nesting = 0
        comment_complete = False
        while self._cur_char != "":
            if self._cur_char in self.comment_end:
                nesting -= 1
                if nesting <= 0:
                    comment_complete = True
                    self._get_next_char()
                    break
            elif self._cur_char in self.comment_begin:
                nesting += 1
            elif self.capture_comments:
                # dest.write(self._cur_char)
                dest.append(self._cur_char)
            self._get_next_char()
        if self.capture_comments:
            # self.captured_comments.append(dest.getvalue())
            self.captured_comments.append("".join(dest))

