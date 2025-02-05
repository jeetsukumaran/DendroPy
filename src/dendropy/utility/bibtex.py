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
BibTeX interface.
"""

import re
import textwrap
from dendropy.utility.container import OrderedCaselessDict
from dendropy.utility import deprecate

# default order of fields
BIBTEX_FIELDS = [
    'author',
    'year',
    'title',
    'journal',
    'volume',
    'number',
    'editor',
    'booktitle',
    'series',
    'chapter',
    'pages',
    'publisher',
    'institution',
    'address',
    'annote',
    'edition',
    'howpublished',
    'key',
    'month',
    'organization',
    'school',
    'note',
    'abstract',
    'keywords',
    'localfile',
    ]

def _clean_parsed_text(text):
    """
    Strips outer quotes, curly braces, etc.; remove multiple
    consecutive whitespaces, etc.
    """
    if text.startswith('{') and text.endswith('}'):
        text = text[1:-1]
    elif text.startswith('"') and text.endswith('"'):
        text = text[1:-1]
    text = re.sub(r"[\s]+", " ", text).strip()
    return text

def _format_bibtex_value(text, col_start=1, wrap=True, width=78):
    """
    Formats text of a BibTeX field.
    """
    ftext = re.sub(r"[\s]+", " ", text).strip()
    col_indent = " " * col_start
    if not ftext[0].isdigit():
        if wrap:
            initial_indent = '{'
            subsequent_indent = ' '
        else:
            ftext = '{' + ftext
        ftext = ftext + '}'
    else:
        initial_indent = ''
        subsequent_indent = ''
    if wrap:
        wrapped = textwrap.wrap(ftext,
                                width=width,
                                initial_indent=initial_indent,
                                subsequent_indent=subsequent_indent)
        for index, line in enumerate(wrapped[1:]):
            wrapped[index+1] = (" " * col_start) + wrapped[index+1]
        return '\n'.join(wrapped)
    else:
        return ftext

class BibTexEntry(object):
    """
    Tracks a single BibTeX entry.
    """
    decompose_pattern = re.compile(r'^@(\w*)\s*{\s*([\w|\:|\-]*),(.*)}')
    # works, but misses last field
    field_pattern = re.compile(r'\s*([\w|\-]*?)\s*=\s*(.*?),(?=\s*[\w|\-]*\s*\=)')
    # get the last field
    last_field_pattern = re.compile(r'\s*([\w|\-]*?)\s*=\s*(.*?)\s*[,]*\s*$')

    def __init__(self, citation=None):
        """
        Sets up internal dictionary of BibTeX fields, and initializes
        if argument is given.
        """

        deprecate.dendropy_deprecation_warning(
            message=(
                "BibTexEntry is deprecated as of DendroPy 5.0 and will be removed in a future release ."
                "It no longer maintaned, and known to be broken. "
                "If this functionality is needed, please open an issue on GitHub."
            ),
        )
        
        self.bibtype = None
        self.citekey = None
        if isinstance(citation, BibTexEntry):
            self._entry_dict = OrderedCaselessDict(citation._entry_dict)
        elif isinstance(citation, dict):
            self._entry_dict = OrderedCaselessDict()
            for k, v in citation.items():
                self._entry_dict[k.lower()] = v
            self.bibtype = self._entry_dict.get("bibtype", None)
            self.citekey = self._entry_dict.get("citekey", None)
        else:
            self._entry_dict = OrderedCaselessDict()
            self.parse_text(citation)

    def __getattr__(self, name):
        """
        Allows bibtex fields (and any additional ones) to be treated
        like object attributes.
        """
        entry_dict = self._get_entry_dict()
        if name == '_entry_dict' or name == '_BibTexEntry_entry_dict':
            return entry_dict
        elif name == '__dict__':
            return object.__getattribute__(self, '__dict__')
        elif name == 'bibtype' and hasattr(self, 'bibtype'):
            return object.__getattribute__(self, '__dict__')['bibtype']
        elif name == 'citekey' and hasattr(self, 'citekey'):
            return object.__getattribute__(self, '__dict__')['citekey']
        elif name in entry_dict:
            return entry_dict[name]
        elif name in BIBTEX_FIELDS:
            return ""
        else:
            raise AttributeError(name)

    def __setattr__(self, name, value):
        """
        Allows bibtex fields (and any additional ones) to be treated
        like object attributes.
        """
        entry_dict = self._get_entry_dict()
        if name == '_entry_dict' or name == '_BibTexEntry_entry_dict':
            entry_dict = value
        elif name == 'bibtype' or name == 'citekey':
            object.__setattr__(self, name, value)
        else:
            self._entry_dict[name] = value

    def __delattr__(self, name):
        """
        Allows bibtex fields (and any additional ones) to be treated
        like object attributes.
        """
        entry_dict = self._get_entry_dict()
        if name == '_entry_dict' or name == '_BibTexEntry_entry_dict':
            object.__delattr__(self, '_entry_dict')
        elif name in entry_dict:
            del(entry_dict[name])
        elif name in BIBTEX_FIELDS:
            pass
        elif name in object.__getattribute__(self, '__dict__'):
            object.__delattr__(name)
        else:
            raise AttributeError(name)

    def __str__(self):
        """
        String representation of self.
        """
        return self.as_bibtex()

    def __repr__(self):
        """
        Internal representation of self.
        """
        repr_dict = {}
        repr_dict['bibtype'] = self.bibtype
        repr_dict['citekey'] = self.citekey
        repr_dict.update(self.fields_as_dict())
        return repr_dict

    def _get_entry_dict(self):
        """
        Returns the internal field dictionary, creating it first if
        neccessary.
        """
        if not hasattr(self, '_entry_dict'):
            object.__setattr__(self, '_entry_dict', {})
        return object.__getattribute__(self, '_entry_dict')

    def _get_fields(self):
        """
        Returns list of populated fields in order (does not include
        bibtype and citekey).
        """
        fields = []
        for field in BIBTEX_FIELDS:
            if field in self._entry_dict:
                fields.append(field)
        for key in self._entry_dict:
            if key not in fields:
                fields.append(key)
        return fields

    fields = property(_get_fields)

    def parse_text(self, text):
        """
        Parses a BibTeX text entry.
        """
        text = text.replace("\n", "")
        self.bibtype = None
        self.citekey = None
        text = text.strip()
        decompose_match = self.decompose_pattern.match(text)
        try:
            self.bibtype = decompose_match.group(1)
        except AttributeError as exception:
            raise ValueError("Failed to parse bibtype: {}".format(text))
        try:
            self.citekey = decompose_match.group(2)
        except AttributeError as exception:
            raise ValueError("Failed to parse citekey: {}".format(text))
        remaining = decompose_match.group(3)
        field_match = self.field_pattern.match(remaining)
        while field_match:
            field_match = self.field_pattern.match(remaining)
            if field_match:
                field_name = field_match.group(1).lower()
                field_value = _clean_parsed_text(field_match.group(2))
                self._entry_dict[field_name] = field_value
                remaining = remaining.replace(field_match.group(), '')
        if remaining:
            last_field_match = self.last_field_pattern.match(remaining)
        if last_field_match:
            field_name = last_field_match.group(1).lower()
            field_value = _clean_parsed_text(last_field_match.group(2))
            self._entry_dict[field_name] = field_value

    def fields_as_dict(self):
        """
        Returns the fields (i.e., all public attributes except for
        bibtype and citekey as a dictionary).
        """
        return dict(self._entry_dict)

    def as_bibtex(self, wrap_width=78):
        """
        Composes entry in BibTex format.
        """
        entry = []
        sep = "  =  "
        entry.append('@{}{{},'.format((self.bibtype, self.citekey)))
        fields = self.fields
#         maxlen = max([len(field) for field in fields])
        maxlen = max([len(field) for field in BIBTEX_FIELDS])
        for field in fields:
            if field != 'url':
                wrap = True
            else:
                wrap = False
            field_header = field.ljust(maxlen)
            field_value = _format_bibtex_value(self._entry_dict[field],
                                      wrap=wrap,
                                      width = wrap_width - maxlen - len(sep) + 2,
                                      col_start = maxlen + len(sep) + 2 )
            entry.append("  {}{}{},".format((field_header, sep, field_value)))
        entry.append('}')
        return '\n'.join(entry)

    def as_compact_bibtex(self):
        """
        Composes entry in BibTex format.
        """
        entry = []
        entry.append('@{}{{{},'.format((self.bibtype, self.citekey)))
        fields = self.fields
        for field in fields:
            field_value = _format_bibtex_value(self._entry_dict[field],
                                      wrap=False,
                                      width=None,
                                      col_start=1)
            entry.append("{}={},".format((field, field_value)))
        entry.append('}')
        return ''.join(entry)

class BibTexDb(object):
    """
    A BibTeX database.
    """


    def __init__(self, bibfilepath=None):
        """
        Initializes database, optionally from file.
        """

        deprecate.dendropy_deprecation_warning(
            message=(
                "BibTexEntry is deprecated as of DendroPy 5.0 and will be removed in a future release ."
                "It no longer maintaned, and known to be broken. "
                "If this functionality is needed, please open an issue on GitHub."
            ),
        )
        
        self.entries = []
        if bibfilepath:
           self.add_from_file(bibfilepath)

    def add_from_file(self, filepath):
        """
        Reads and loads a BibTeX database file.
        """
        bibfile = open(filepath, 'r')
        self.add_from_text(bibfile.read())
        bibfile.close()

    def add_from_text(self, text):
        """
        Loads from text.
        """
        text = text.replace('\n','')
        entry_pattern = re.compile(r'@\w*{([^{}]+{[^{}]*({[^{}]*}[^{}]*)*})+}')
        for match in entry_pattern.finditer(text):
            entry = BibTexEntry(match.group())
            self.entries.append(entry)
