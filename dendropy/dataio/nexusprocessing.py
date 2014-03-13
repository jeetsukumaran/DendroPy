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

    def __init__(self, taxon_namespace, case_insensitive=True):
        self._taxon_namespace = None
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

    def _get_taxon_namespace(self):
        return self.taxon_namespace
    def _set_taxon_namespace(self, taxon_namespace):
        self._taxon_namespace = taxon_namespace
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
        self.token_taxon_map[token] = taxon

    def lookup_taxon_symbol(self, symbol, create_taxon_if_not_found=True):
        # if symbol in self.token_taxon_map:
        #     return self.token_taxon_map[symbol]
        # if symbol in self.label_taxon_map:
        #     return self.label_taxon_map[symbol]
        # if symbol in self.number_taxon_map:
        #     return self.number_taxon_map[symbol]
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
            t = self._taxon_namespace.new_taxon(symbol)
            self.label_taxon_map[symbol] = t
            return t
        return None

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


