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
Specialized tokenizer for processing NEXUS/Newick streams.
"""

import re
import itertools
import decimal
from dendropy.dataio.tokenizer import Tokenizer
from dendropy.utility import textprocessing
from dendropy.utility import container
from dendropy.datamodel import basemodel

##############################################################################
## NexusTokenizer

class NexusTokenizer(Tokenizer):

    def __init__(self, src,
            preserve_unquoted_underscores=False):
        Tokenizer.__init__(self,
            src=src,
            uncaptured_delimiters=set(" \t\n\r"),
            captured_delimiters=set(r'{}(),;:=\"'),
            quote_chars=set("'"),
            escape_quote_by_doubling=True,
            escape_chars=set(""),
            comment_begin=set("["),
            comment_end=set("]"),
            capture_comments=True,
            preserve_unquoted_underscores=preserve_unquoted_underscores)
        # self.preserve_unquoted_underscores = preserve_unquoted_underscores

    # def __next__(self):
        # Tokenizer.__next__(self)
        # if (not self.preserve_unquoted_underscores
        #         and not self.is_token_quoted):
        #     self.current_token = self.current_token.replace("_", " ")
        # return self.current_token

    def set_capture_eol(self, capture_eol):
        if capture_eol:
            try:
                self.uncaptured_delimiters.discard("\n")
            except ValueError:
                pass
            try:
                self.uncaptured_delimiters.discard("\r")
            except ValueError:
                pass
            if "\n" not in self.captured_delimiters:
                self.captured_delimiters.add("\n")
            if "\r" not in self.captured_delimiters:
                self.captured_delimiters.add("\r")
        else:
            try:
                self.captured_delimiters.discard("\n")
            except ValueError:
                pass
            try:
                self.captured_delimiters.discard("\r")
            except ValueError:
                pass
            if "\n" not in self.uncaptured_delimiters:
                self.uncaptured_delimiters.add("\n")
            if "\r" not in self.uncaptured_delimiters:
                self.uncaptured_delimiters.add("\r")

    def set_hyphens_as_captured_delimiters(self, hyphens_as_captured_delimiters):
        if hyphens_as_captured_delimiters:
            if "-" not in self.captured_delimiters:
                self.captured_delimiters.add("-")
        else:
            try:
                self.captured_delimiters.discard("-")
            except ValueError:
                pass

    def require_next_token_ucase(self):
        t = self.require_next_token()
        t = t.upper()
        self.current_token = t
        return t

    def next_token_ucase(self):
        try:
            t = self.__next__()
            t = t.upper()
            self.current_token = t
            return t
        except StopIteration:
            self.current_token = None
            return None

    def cast_current_token_to_ucase(self):
        if self.current_token:
            self.current_token = self.current_token.upper()
        return self.current_token

    def process_and_clear_comments_for_item(self,
            item,
            extract_comment_metadata):
        process_comments_for_item(item,
                self.captured_comments,
                extract_comment_metadata)
        del self.captured_comments[:]

    def skip_to_semicolon(self):
        token = self.next_token()
        while token != ';' and not self._cur_char == "" and token != None:
            token = self.next_token()

###############################################################################
## Taxon Handling

class NexusTaxonSymbolMapper(object):
    """
    Manages |TaxonNamespace| and |Taxon| object look-ups when
    parsing NEXUS and NEWICK formatted data.

    Operational taxonomic unit concepts in NEXUS files can be referenced using
    one of three types of symbols:

        - the "TRANSLATE" block token
        - the taxon label
        - the taxon number

    In the event of redundant over overdetermined symbols, the resolution order
    is as given above.

    This class encapsulates creating, looking-up and retrieving |Taxon|
    objects corresponding to operation taxonomic unit concept references
    encountered when reading NEXUS or NEWICK data sources from the
    |TaxonNamespace| that it wraps and manages. It keeps track of
    "TRANSLATE" block tokens, operational taxonomic unit labels, and
    operational taxonomic unit indexes in mapping containers that allow for
    quick retrieval of corresponding |Taxon| objects. The symbol look-up
    is case-insensitive, as per NEXUS/NEWICK convention.

    If a |Taxon| object is not found for a particular symbol, it will
    create a new |Taxon| object with that symbol for its label, and
    register it in all the other supplemental mappings appropriately.

    Note that the |TaxonNamespace| object passed to this class and the
    member |Taxon| objects should not be modified during the lifespan of
    this class or, at least, the tenure of the management of
    |TaxonNamespace| and member |Taxon| objects by this class.
    This is to ensure that the various supplementatl mappings (in particular,
    the label mapping and the taxon number mapping) are synchronized.
    To this end, the of the |TaxonNamespace| object is locked, and all
    |Taxon| object creation should be through this class's native
    methods.
    """

    def __init__(self,
            taxon_namespace,
            enable_lookup_by_taxon_number=True,
            case_sensitive=False):
        self._taxon_namespace = None
        self.taxon_namespace_original_mutability_state = None
        self.case_sensitive = case_sensitive
        if not self.case_sensitive:
            if taxon_namespace.is_case_sensitive:
                raise ValueError("Attempting case insensitive read with case sensitive TaxonNamespace")
            self.token_taxon_map = container.CaseInsensitiveDict()
            self.label_taxon_map = container.CaseInsensitiveDict()
        else:
            if not taxon_namespace.is_case_sensitive:
                raise ValueError("Attempting case sensitive read with case insensitive TaxonNamespace")
            self.token_taxon_map = {}
            self.label_taxon_map = {}
        self.number_taxon_map = {}
        self.number_taxon_label_map = {}
        self.enable_lookup_by_taxon_number = enable_lookup_by_taxon_number
        self._set_taxon_namespace(taxon_namespace)

    def restore_taxon_namespace_mutability(self):
        if self._taxon_namespace is not None:
            if self.taxon_namespace_original_mutability_state is not None:
                self._taxon_namespace.is_mutable = self.taxon_namespace_original_mutability_state
        self.taxon_namespace_original_mutability_state = None

    def __del__(self):
        self.restore_taxon_namespace_mutability()

    def _get_taxon_namespace(self):
        return self._taxon_namespace
    def _set_taxon_namespace(self, taxon_namespace):
        if self._taxon_namespace is not None:
            self.restore_taxon_namespace_mutability()
        self._taxon_namespace = taxon_namespace
        self.taxon_namespace_original_mutability_state = self._taxon_namespace.is_mutable
        self._taxon_namespace.is_mutable = False
        self.reset_supplemental_mappings()
    taxon_namespace = property(_get_taxon_namespace, _set_taxon_namespace)

    def reset_supplemental_mappings(self):
        self.token_taxon_map.clear()
        if not self.case_sensitive:
            self.label_taxon_map = container.CaseInsensitiveDict(self._taxon_namespace.label_taxon_map())
        else:
            self.label_taxon_map = self._taxon_namespace.label_taxon_map()
        self.number_taxon_map.clear()
        for idx, taxon in enumerate(self._taxon_namespace):
            s = str(idx+1)
            self.number_taxon_map[s] = taxon
            self.number_taxon_label_map[s] = taxon.label

    def add_translate_token(self, token, taxon):
        if not textprocessing.is_str_type(token):
            token = str(token)
        self.token_taxon_map[token] = taxon

    def lookup_taxon_symbol(self, symbol, create_taxon_if_not_found=True):
        # if symbol in self.token_taxon_map:
        #     return self.token_taxon_map[symbol]
        # if symbol in self.label_taxon_map:
        #     return self.label_taxon_map[symbol]
        # if symbol in self.number_taxon_map:
        #     return self.number_taxon_map[symbol]
        if not textprocessing.is_str_type(symbol):
            symbol = str(symbol)
        try:
            return self.token_taxon_map[symbol]
        except KeyError:
            pass
        try:
            return self.label_taxon_map[symbol]
        except KeyError:
            pass
        if self.enable_lookup_by_taxon_number:
            try:
                return self.number_taxon_map[symbol]
            except KeyError:
                pass
        if create_taxon_if_not_found:
            return self.new_taxon(symbol)
        return None

    def require_taxon_for_symbol(self, symbol):
        return self.lookup_taxon_symbol(symbol=symbol, create_taxon_if_not_found=True)

    def new_taxon(self, label):
        self._taxon_namespace.is_mutable = self.taxon_namespace_original_mutability_state
        t = self._taxon_namespace.new_taxon(label)
        self._taxon_namespace.is_mutable = False
        self.label_taxon_map[label] = t
        taxon_number = str(len(self._taxon_namespace))
        self.number_taxon_map[taxon_number] = t
        return t

    def add_taxon(self, taxon):
        self._taxon_namespace.is_mutable = self.taxon_namespace_original_mutability_state
        self._taxon_namespace.add_taxon(taxon)
        self._taxon_namespace.is_mutable = False
        self.label_taxon_map[taxon.label] = taxon
        taxon_number = str(len(self._taxon_namespace))
        self.number_taxon_map[taxon_number] = taxon
        return taxon

###############################################################################
## Metadata

FIGTREE_COMMENT_FIELD_PATTERN = re.compile(r'(.+?)=({.+?,.+?}|.+?)(,|$)')
NHX_COMMENT_FIELD_PATTERN = re.compile(r'(.+?)=({.+?,.+?}|.+?)(:|$)')

def parse_comment_metadata_to_annotations(
        comment,
        annotations=None,
        field_name_map=None,
        field_value_types=None,
        strip_leading_trailing_spaces=True):
    """
    Returns set of |Annotation| objects corresponding to metadata
    given in comments.

    Parameters
    ----------
    ``comment`` : string
        A comment token.
    ``annotations`` : |AnnotationSet| or ``set``
        Set of |Annotation| objects to which to add this annotation.
    ``field_name_map`` : dict
        A dictionary mapping field names (as given in the comment string)
        to strings that should be used to represent the field in the
        metadata dictionary; if not given, no mapping is done (i.e., the
        comment string field name is used directly).
    ``field_value_types`` : dict
        A dictionary mapping field names (as given in the comment
        string) to the value type (e.g. {"node-age" : float}.
    ``strip_leading_trailing_spaces`` : boolean
        Remove whitespace from comments.

    Returns
    -------
    metadata : :py:``set`` [|Annotation|]
        Set of |Annotation| objects corresponding to metadata
        parsed.
    """
    if annotations is None:
        annotations = set()
    if field_name_map is None:
        field_name_map = {}
    if field_value_types is None:
        field_value_types = {}
    if comment.startswith("&&NHX:"):
        pattern = NHX_COMMENT_FIELD_PATTERN
        comment = comment[6:]
    elif comment.startswith("&&"):
        pattern = NHX_COMMENT_FIELD_PATTERN
        comment = comment[2:]
    elif comment.startswith("&"):
        pattern = FIGTREE_COMMENT_FIELD_PATTERN
        comment = comment[1:]
    else:
        # unrecognized metadata pattern
        return annotations
    for match_group in pattern.findall(comment):
        key, val = match_group[:2]
        if strip_leading_trailing_spaces:
            key = key.strip()
            val = val.strip()
        if key in field_value_types:
            value_type = field_value_types[key]
        else:
            value_type = None
        if val.startswith('{'):
            if value_type is not None:
                val = [value_type(v) for v in val[1:-1].split(',')]
            else:
                val = val[1:-1].split(',')
        elif val.startswith('"') and val.endswith('"'):
            val = val[1:-1]
        elif val.lower() == "false":
            val = False
        elif val.lower() == "true":
            val = True
        else:
            if value_type is not None:
                val = value_type(val)
        if key in field_name_map:
            key = field_name_map[key]
        annote = basemodel.Annotation(
                name=key,
                value=val,
                # datatype_hint=datatype_hint,
                # name_prefix=name_prefix,
                # namespace=namespace,
                # name_is_prefixed=name_is_prefixed,
                # is_attribute=False,
                # annotate_as_reference=annotate_as_reference,
                # is_hidden=is_hidden,
                )
        annotations.add(annote)
    return annotations

def process_comments_for_item(item,
        item_comments,
        extract_comment_metadata):
    if not item_comments or item is None:
        return
    for comment in item_comments:
        if extract_comment_metadata and comment.startswith("&"):
            annotations = parse_comment_metadata_to_annotations(comment)
            if annotations:
                item.annotations.update(annotations)
            else:
                item.comments.append(comment)
        else:
            item.comments.append(comment)

###############################################################################
## NEWICK/NEXUS formatting support.

def format_annotated_value(
        value,
        annotated_real_value_format_specifier,
        real_value_format_specifier=None,
        override_annotation_format_specifier=False):
    if isinstance(value, float) or isinstance(value, decimal.Decimal):
        if override_annotation_format_specifier and real_value_format_specifier is not None:
            fmtspec = real_value_format_specifier
        elif annotated_real_value_format_specifier is not None:
            fmtspec = annotated_real_value_format_specifier
        elif real_value_format_specifier is not None:
            fmtspec = real_value_format_specifier
        else:
            fmtspec = ""
    else:
        if annotated_real_value_format_specifier is not None:
            fmtspec = annotated_real_value_format_specifier
        else:
            fmtspec = ""
    return "{:{fmtspec}}".format(value, fmtspec=fmtspec)

def format_item_annotations_as_comments(
        annotated,
        nhx=False,
        real_value_format_specifier=None,
        override_annotation_format_specifier=False):
    """
    ``annotated`` - Annotated object
    ``nhx``       - render as NHX '[&& ...]'? Otherwise as '[& ...]'
    ``real_value_format_specifier`` - Format specification for real/float values.
                   The format specifier should be given in Python's string
                   format specification mini-language. E.g. ".8f", ".4E",
                   "8.4f". If the annotation has its ``format_specifier``
                   attribute set, then this argument is ignored for rendering
                   that particular annotation unless
                   ``override_annotation_format_specifier`` is |True|. Defaults to "".
    ``override_annotation_format_specifier``
                    If the annotation has its ``format_specifier`` attribute set,
                    then this it will be used in preference to the
                    ``real_value_format_specifier`` above unless this argument is
                    |True|. Defaults to |False|.
    """
    if not annotated.annotations:
        return ""
    parts = []
    for annote in annotated.annotations:
        if annote.is_hidden:
            continue
        key = annote.name
        value = annote.value
        if isinstance(value, list) or isinstance(value, tuple):
            items = []
            for item in value:
                items.append(format_annotated_value(
                    value=item,
                    annotated_real_value_format_specifier=annote.real_value_format_specifier,
                    real_value_format_specifier=real_value_format_specifier,
                    override_annotation_format_specifier=override_annotation_format_specifier))
            items = ",".join(items)
            parts.append("%s={%s}" % (key, items))
        elif isinstance(value, dict):
            raise TypeError("Dictionary types not supported for rendering as NEWICK-formatted metadata")
        else:
            x = format_annotated_value(
                    value=value,
                    annotated_real_value_format_specifier=annote.real_value_format_specifier,
                    real_value_format_specifier=real_value_format_specifier,
                    override_annotation_format_specifier=override_annotation_format_specifier)
            parts.append("{key}={value}".format(key=key, value=x))
    if nhx:
        prefix = "[&&NHX:"
        separator = ":"
        suffix = "]"
    else:
        prefix = "[&"
        separator = ","
        suffix = "]"
    body = separator.join(parts)
    return prefix + body + suffix

def escape_nexus_token(
    label,
    preserve_spaces=False,
    quote_underscores=True,
    protect_regex=r'''[()[\]{}\\\/,;:=*'"`+\-<>\0\t\n]''',
):
    """
    Properly protects a NEXUS token.

    Kwarg protect_regex allows less eager quoting when working with non-Nexus
    Newick strings.
    """
    if label is None:
        return ""
    if not preserve_spaces \
            and "_" not in label \
            and not re.search(protect_regex, label):
        label = label.replace(' ', '_').replace('\t', '_')
    elif re.search(protect_regex, label) \
        or " " in label \
        or (quote_underscores and "_" in label):
        s = label.split("'")
        if len(s) == 1:
            return "'" + label + "'"
        return "'{}'".format("''".join(s))
    return label

def bitmask_as_newick_string(split, taxon_set, preserve_spaces=False, quote_underscores=True):
    """
    Represents a split as a newick string.
    """
    taxlabels = [escape_nexus_token(label, preserve_spaces=preserve_spaces, quote_underscores=quote_underscores) for label in taxon_set.labels()]

    # do not do the root
    if split == 0 or (split == taxon_set.all_taxa_bitmask()):
        return "({});".format(",".join(taxlabels))

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
    return "(({}), ({}));".format(", ".join(left), ", ".join(right))

def group_ranges(L):
    """
    Collapses a list of integers into a list of the start and end of
    consecutive runs of numbers. Returns a generator of generators.

    >>> [list(x) for x in group_ranges([1, 2, 3, 5, 6, 8])]
    [[1, 3], [5, 6], [8]]
    """
    for w, z in itertools.groupby(sorted(L), lambda x, y=itertools.count(): next(y)-x):
        grouped = list(z)
        yield (x for x in [grouped[0], grouped[-1]][:len(grouped)])

def get_rooting_argument(**kwargs):
    if "is_rooted" in kwargs and "is_unrooted" in kwargs:
        raise ValueError("Must specify exactly one of 'is_rooted' or 'is_unrooted'")
    elif "is_rooted":
        is_rooted = kwargs["is_rooted"]
    elif "is_unrooted":
        is_rooted = not kwargs["is_unrooted"]
    if is_rooted is True:
        return "force-rooted"
    elif is_rooted is False:
        return "force-unrooted"
    else:
        return "default-unrooted"
