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
from dendropy.dataio import ioservice
from dendropy.utility import text

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

    def _write(self,
            stream,
            taxon_namespaces=None,
            tree_lists=None,
            char_matrices=None,
            global_annotations_target=None):
        for char_matrix in char_matrices:
            if (self.attached_taxon_namespace is not None
                    and char_matrix.taxon_namespace is not self.attached_taxon_namespace):
                continue
            self._write_char_matrix(stream, char_matrix)

    def _write_char_matrix(self, stream, char_matrix):
        "Writes dataset to a full PHYLIP document."

        if self.strict or self.force_unique_taxon_labels:
            taxon_label_map = self.get_taxon_label_map(char_matrix.taxon_namespace)
            if not self.strict:
                spacer = "  "
            else:
                spacer = ""
        else:
            taxon_label_map = {}
            for taxon in char_matrix.taxon_namespace:
                label = taxon.label
                if self.spaces_to_underscores:
                    label = label.replace(' ', '_')
                taxon_label_map[taxon] = label
            spacer = "  "
        maxlen = max([len(str(label)) for label in taxon_label_map.values()])
        n_seqs = len(char_matrix)
        n_sites = len(char_matrix.values()[0])
        stream.write("%d %d\n" % (n_seqs, n_sites))
        for taxon in char_matrix.taxon_namespace:
            label = taxon_label_map[taxon]
            try:
                seq_vec = char_matrix[taxon]
            except KeyError:
                continue
            stream.write("%s%s%s\n" % ( label.ljust(maxlen), spacer, str(seq_vec.symbols_as_string())))

    def get_taxon_label_map(self, taxon_namespace):
        taxon_label_map = {}
        if self.strict:
            max_label_len = STRICT_MODE_MAX_LABEL_LENGTH
        else:
            max_label_len = 0
        for taxon in taxon_namespace:
            label = taxon.label
            if self.spaces_to_underscores:
                label = label.replace(' ', '_')
            if self.strict:
                label = label[:max_label_len]
            taxon_label_map[taxon] = label
        taxon_label_map = text.unique_taxon_label_map(taxon_namespace, taxon_label_map, max_label_len, _LOG)
        if self.strict:
            for t in taxon_label_map:
                label = taxon_label_map[t]
                if len(label) < STRICT_MODE_MAX_LABEL_LENGTH:
                    taxon_label_map[t] = label.ljust(STRICT_MODE_MAX_LABEL_LENGTH)
        return taxon_label_map

