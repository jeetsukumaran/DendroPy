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
Implementation of PHYLIP-format data writer.
"""

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.utility import textutils


STRICT_MODE_MAX_LABEL_LENGTH = 10

class PhylipWriter(ioservice.DataWriter):
    "Implements the DataWriter interface for writing PHYLIP files."

    def __init__(self, **kwargs):
        """
        __init__ recognizes the following keywords (in addition to those of `DataWriter.__init__`):

            - `strict` (boolean)
            - `spaces_to_underscores` (boolean)
            - `force_unique_taxon_labels` (boolean)
        """
        ioservice.DataWriter.__init__(self, **kwargs)
        self.strict = kwargs.get("strict", False)
        self.spaces_to_underscores = kwargs.get("spaces_to_underscores", False)
        self.force_unique_taxon_labels = kwargs.get("force_unique_taxon_labels", False)

    def get_taxon_label_map(self, taxon_set):
        taxon_label_map = {}
        if self.strict:
            max_label_len = STRICT_MODE_MAX_LABEL_LENGTH
        else:
            max_label_len = 0
        for taxon in taxon_set:
            label = taxon.label
            if self.spaces_to_underscores:
                label = label.replace(' ', '_')
            if self.strict:
                label = label[:max_label_len]
            taxon_label_map[taxon] = label
        taxon_label_map = textutils.unique_taxon_label_map(taxon_set, taxon_label_map, max_label_len, _LOG)
        if self.strict:
            for t in taxon_label_map:
                label = taxon_label_map[t]
                if len(label) < STRICT_MODE_MAX_LABEL_LENGTH:
                    taxon_label_map[t] = label.ljust(STRICT_MODE_MAX_LABEL_LENGTH)
        return taxon_label_map

    def write(self, stream):
        "Writes dataset to a full PHYLIP document."

        if self.exclude_chars:
            return self.dataset

        assert self.dataset is not None, \
            "PhylipWriter instance is not attached to a DataSet: no source of data"

        char_matrix = None
        if len(self.dataset.char_matrices) == 0:
            raise ValueError("No character data in DataSet")
        if self.attached_taxon_set is not None:
            taxon_set_matrices = [cmat for cmat in self.dataset.char_matrices if cmat.taxon_set is self.attached_taxon_set]
            if len(taxon_set_matrices) == 0:
                raise ValueError("No character matrix associated with attached TaxonSet '%s'" % (repr(self.attached_taxon_set)))
            if len(taxon_set_matrices) > 1:
                raise ValueError("Multiple character matrices associated with attached TaxonSet '%s'" % (repr(self.attached_taxon_set)))
            char_matrix = taxon_set_matrices_map[self.attached_taxon_set]
        else:
            if len(self.dataset.char_matrices) > 1:
                raise ValueError("Multiple character matrices found")
            char_matrix = self.dataset.char_matrices[0]

        assert char_matrix is not None, \
            "Failed to identify suitable CharacterMatrix"

        if self.strict or self.force_unique_taxon_labels:
            taxon_label_map = self.get_taxon_label_map(char_matrix.taxon_set)
            if not self.strict:
                spacer = "  "
            else:
                spacer = ""
        else:
            taxon_label_map = {}
            for taxon in char_matrix.taxon_set:
                label = taxon.label
                if self.spaces_to_underscores:
                    label = label.replace(' ', '_')
                taxon_label_map[taxon] = label
            spacer = "  "
        maxlen = max([len(str(label)) for label in taxon_label_map.values()])

        n_seqs = len(char_matrix)
        n_sites = len(char_matrix.values()[0])
        stream.write("%d %d\n" % (n_seqs, n_sites))

        for taxon in char_matrix.taxon_set:
            label = taxon_label_map[taxon]
            try:
                seq_vec = char_matrix[taxon]
            except KeyError:
                continue
            stream.write("%s%s%s\n" % ( label.ljust(maxlen), spacer, str(seq_vec.symbols_as_string()).replace(' ', '')))
