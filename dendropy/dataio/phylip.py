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
Implementation of PHYLIP-schema i/o client(s).
"""

from dendropy.utility import iosys
from dendropy.utility import texttools
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

class PhylipWriter(iosys.DataWriter):
    "Implements the DataWriter interface for handling PHYLIP files."

    def __init__(self, **kwargs):
        "Calls the base class constructor."
        iosys.DataWriter.__init__(self, **kwargs)
        self.strict = kwargs.get("strict", False)
        self.spaces_to_underscores = kwargs.get("spaces_to_underscores", False)
        self.taxon_label_map = {}

    def get_taxon_label_map(self, dataset):
        if self.attached_taxon_set is not None:
            taxon_sets = [self.attached_taxon_set]
        else:
            taxon_sets = dataset.taxon_sets
        taxon_label_map = {}
        taxa = []
        if self.strict:
            max_label_len = 10
        else:
            max_label_len = 0
        for taxon_set in taxon_sets:
            for taxon in taxon_set:
                label = taxon.label
                if self.spaces_to_underscores:
                    label = label.replace(' ', '')
                if self.strict:
                    label = label[:max_label_len]
                taxa.append(taxon)
                taxon_label_map[taxon] = label
        taxon_label_map = texttools.unique_taxon_label_map(taxa, taxon_label_map, max_label_len, _LOG)
        return taxon_label_map


    def write(self, stream, **kwargs):
        "Writes dataset to a full PHYLIP document."

        # update directives
        self.dataset = kwargs.get("dataset", self.dataset)
        self.attached_taxon_set = kwargs.get("taxon_set", self.attached_taxon_set)
        self.exclude_trees = kwargs.get("exclude_trees", self.exclude_trees)
        self.exclude_chars = kwargs.get("exclude_chars", self.exclude_chars)
        self.strict = kwargs.get("strict", False)
        self.spaces_to_underscores = kwargs.get("spaces_to_underscores", False)
        self.taxon_label_map = {}

        if self.exclude_chars:
            return self.dataset

        if self.strict:
            self.taxon_label_map = self.get_taxon_label_map(self.dataset)

        assert self.dataset is not None, \
            "PhylipWriter instance is not attached to a DataSet: no source of data"

        char_matrix = self.dataset.char_matrices[0]
        n_seqs = len(char_matrix)
        n_sites = len(char_matrix.values()[0])
        stream.write("%d %d\n" % (n_seqs, n_sites))
        maxlen = max([len(str(taxon)) for taxon in char_matrix])
        for taxon in char_matrix.taxon_set:
            stream.write("%s        %s\n" % ( str(taxon).ljust(maxlen), str(char_matrix[taxon].symbols_as_string()).replace(' ', '')))
