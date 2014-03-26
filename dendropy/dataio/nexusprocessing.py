#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Specialized tokenizer for processing NEXUS/Newick streams.
"""

import re
import io
from dendropy.dataio.tokenizer import Tokenizer
from dendropy.utility import container

##############################################################################
## NexusTokenizer

class NexusTokenizer(Tokenizer):

    def __init__(self, src):
        Tokenizer.__init__(self,
            src=src,
            uncaptured_delimiters=" \t\n\r",
            captured_delimiters="(),;:",
            quote_chars="\"'",
            escape_quote_by_doubling=True,
            escape_chars="",
            comment_begin="[",
            comment_end="]",
            capture_comments=True)

###############################################################################
## Taxon Handling

class NexusTaxonSymbolMapper(object):
    """
    Manages :class:`TaxonNamespace` and :class:`Taxon` object look-ups when
    parsing NEXUS and NEWICK formatted data.

    Operational taxonomic unit concepts in NEXUS files can be referenced using
    one of three types of symbols:

        - the "TRANSLATE" block token
        - the taxon label
        - the taxon number

    In the event of redundant over overdetermined symbols, the resolution order
    is as given above.

    This class encapsulates creating, looking-up and retrieving :class:`Taxon`
    objects corresponding to operation taxonomic unit concept references
    encountered when reading NEXUS or NEWICK data sources from the
    :class:`TaxonNamespace` that it wraps and manages. It keeps track of
    "TRANSLATE" block tokens, operational taxonomic unit labels, and
    operational taxonomic unit indexes in mapping containers that allow for
    quick retrieval of corresponding :class:`Taxon` objects. The symbol look-up
    is case-insensitive, as per NEXUS/NEWICK convention.

    If a :class:`Taxon` object is not found for a particular symbol, it will
    create a new :class:`Taxon` object with that symbol for its label, and
    register it in all the other supplemental mappings appropriately.

    Note that the :class:`TaxonNamespace` object passed to this class and the
    member :class:`Taxon` objects should not be modified during the lifespan of
    this class or, at least, the tenure of the management of
    :class:`TaxonNamespace` and member :class:`Taxon` objects by this class.
    This is to ensure that the various supplementatl mappings (in particular,
    the label mapping and the taxon number mapping) are synchronized.
    To this end, the of the :class:`TaxonNamespace` object is locked, and all
    :class:`Taxon` object creation should be through this class's native
    methods.
    """

    def __init__(self, taxon_namespace, case_insensitive=True):
        self._taxon_namespace = None
        self.taxon_namespace_original_mutability_state = None
        self.case_insensitive = case_insensitive
        if self.case_insensitive:
            self.token_taxon_map = container.CaseInsensitiveDict()
            self.label_taxon_map = container.CaseInsensitiveDict()
        else:
            self.token_taxon_map = {}
            self.label_taxon_map = {}
        self.number_taxon_map = {}
        self.number_taxon_label_map = {}
        self._set_taxon_namespace(taxon_namespace)

    def restore_taxon_namespace_mutability(self):
        if self._taxon_namespace is not None:
            if self.taxon_namespace_original_mutability_state is not None:
                self._taxon_namespace.is_mutable = self.taxon_namespace_original_mutability_state
        self.taxon_namespace_original_mutability_state = None

    def __del__(self):
        self.restore_taxon_namespace_mutability()

    def _get_taxon_namespace(self):
        return self.taxon_namespace
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
        if self.case_insensitive:
            self.label_taxon_map = container.CaseInsensitiveDict(self._taxon_namespace.label_taxon_map())
        else:
            self.label_taxon_map = self._taxon_namespace.label_taxon_map()
        self.number_taxon_map.clear()
        for idx, taxon in enumerate(self._taxon_namespace):
            s = str(idx+1)
            self.number_taxon_map[s] = taxon
            self.number_taxon_label_map[s] = taxon.label

    def add_translate_token(self, token, taxon):
        if not isinstance(token, str):
            token = str(token)
        self.token_taxon_map[token] = taxon

    def lookup_taxon_symbol(self, symbol, create_taxon_if_not_found=True):
        # if symbol in self.token_taxon_map:
        #     return self.token_taxon_map[symbol]
        # if symbol in self.label_taxon_map:
        #     return self.label_taxon_map[symbol]
        # if symbol in self.number_taxon_map:
        #     return self.number_taxon_map[symbol]
        if not isinstance(symbol, str):
            symbol = str(symbol)
        try:
            return self.token_taxon_map[symbol]
        except KeyError:
            pass
        try:
            return self.label_taxon_map[symbol]
        except KeyError:
            pass
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
        self._taxon_namespace.is_mutable = True
        t = self._taxon_namespace.new_taxon(label)
        self._taxon_namespace.is_mutable = False
        self.label_taxon_map[label] = t
        taxon_number = str(len(self._taxon_namespace))
        self.number_taxon_map[taxon_number] = t
        return t

    def add_taxon(self, taxon):
        self._taxon_namespace.is_mutable = True
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
    Returns set of :class:`Annotation` objects corresponding to metadata
    given in comments.

    Parameters
    ----------
    `comment` : string
        A comment token.
    `annotations` : :class:`AnnotationSet` or `set`
        Set of :class:`Annotation` objects to which to add this annotation.
    `field_name_map` : dict
        A dictionary mapping field names (as given in the comment string)
        to strings that should be used to represent the field in the
        metadata dictionary; if not given, no mapping is done (i.e., the
        comment string field name is used directly).
    `field_value_types` : dict
        A dictionary mapping field names (as given in the comment
        string) to the value type (e.g. {"node-age" : float}.
    `strip_leading_trailing_spaces` : boolean
        Remove whitespace from comments.

    Returns
    -------
    metadata : :py:class::`set` [:class:`Annotation`]
        Set of :class:`Annotation` objects corresponding to metadata
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
        else:
            if value_type is not None:
                val = value_type(val)
        if key in field_name_map:
            key = field_name_map[key]
        annote = base.Annotation(
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

###############################################################################
## NEWICK/NEXUS formatting support.

def escape_nexus_token(label, preserve_spaces=False, quote_underscores=True):
    """
    Properly protects a NEXUS token.
    """
    if label is None:
        return ""
    if not preserve_spaces \
            and "_" not in label \
            and not re.search('[\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>\0\t\n]', label):
        label = label.replace(' ', '_').replace('\t', '_')
    elif re.search('[\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>\0\t\n\r ]', label) \
        or quote_underscores and "_" in label:
        s = label.split("'")
        if len(s) == 1:
            return "'" + label + "'"
        return "'{}'".format("''".join(s))
    return label

def split_as_newick_string(split, taxon_set, preserve_spaces=False, quote_underscores=True):
    """
    Represents a split as a newick string.
    """
    taxlabels = [escape_nexus_token(label, preserve_spaces=preserve_spaces, quote_underscores=quote_underscores) for label in taxon_set.labels()]

    # do not do the root
    if split == 0 or (split == taxon_set.all_taxa_bitmask()):
        return "({})".format(",".join(taxlabels))

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
    return "(({}), ({}))".format(", ".join(left), ", ".join(right))

def group_ranges(L):
    """
    Collapses a list of integers into a list of the start and end of
    consecutive runs of numbers. Returns a generator of generators.

    >>> [list(x) for x in group_ranges([1, 2, 3, 5, 6, 8])]
    [[1, 3], [5, 6], [8]]
    """
    for w, z in itertools.groupby(L, lambda x, y=itertools.count(): next(y)-x):
        grouped = list(z)
        yield (x for x in [grouped[0], grouped[-1]][:len(grouped)])


