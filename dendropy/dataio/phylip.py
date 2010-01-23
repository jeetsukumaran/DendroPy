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

from cStringIO import StringIO
from dendropy.utility import iosys
from dendropy.utility import texttools
from dendropy.utility import error
from dendropy.utility.messaging import get_logger
from dendropy import dataobject
_LOG = get_logger(__name__)

STRICT_MODE_MAX_LABEL_LENGTH = 10

class PhylipReader(iosys.DataReader):
    "Implements the DataReader interface for parsing PHYLIP files."

    supported_data_types = ['dna', 'rna', 'protein', 'standard', 'restriction', 'infinite']
    supported_matrix_types = [dataobject.DnaCharacterMatrix,
                              dataobject.RnaCharacterMatrix,
                              dataobject.ProteinCharacterMatrix,
                              dataobject.StandardCharacterMatrix,
                              dataobject.RestrictionSitesCharacterMatrix,
                              dataobject.InfiniteSitesCharacterMatrix]

    def __init__(self, **kwargs):
        iosys.DataReader.__init__(self, **kwargs)
        if "char_matrix_type" in kwargs and "data_type" in kwargs:
            raise ValueError("Cannot specify both 'data_type' and 'char_matrix_type'")
        if "data_type" in kwargs:
            data_type = kwargs["data_type"].lower()
            if data_type not in PhylipReader.supported_data_types:
                raise ValueError("'%s' is not a valid data type specification; must be one of: %s" \
                    % (", ".join([("'" + d + "'") for d in PhylipReader.supported_data_types])))
            else:
                self.char_matrix_type = dataobject.character_data_type_label_map[data_type]
#        elif "char_matrix_type" in kwargs:
#            self.char_matrix_type = kwargs.get("char_matrix_type")
#        else:
#            raise ValueError("Must specify at 'data_type' or 'char_matrix_type'")
        else:
            self.char_matrix_type = kwargs.get("char_matrix_type", dataobject.DnaCharacterMatrix)
        if self.char_matrix_type not in Phylip.supported_matrix_types:
            raise ValueError("'%s' is not a supported data type for PhylipReader" % self.char_matrix_type.__name__)
        self.strict = kwargs.get("strict", False)
        self.underscores_to_spaces = kwargs.get("underscores_to_spaces", False)

    def read(self, stream, **kwargs):
        self.dataset = kwargs.get("dataset", self.dataset)
        self.attached_taxon_set = kwargs.get("taxon_set", self.attached_taxon_set)
        self.exclude_trees = kwargs.get("exclude_trees", self.exclude_trees)
        self.exclude_chars = kwargs.get("exclude_chars", self.exclude_chars)
        self.strict = kwargs.get("strict", False)
        self.underscores_to_spaces = kwargs.get("underscores_to_spaces", False)
        if self.exclude_chars:
            return self.dataset
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        if self.attached_taxon_set is not None \
                and self.attached_taxon_set is not in self.dataset.taxon_sets:
            self.dataset.add(self.attached_taxon_set)
        else:
            self.attached_taxon_set = self.dataset.new_taxon_set()
        if self.strict:
            self._parse_strict(stream)
        else:
            self._parse_relaxed(stream)
        return self.dataset

    def _parse_strict(self, stream):
        raise NotImplementedError()

    def _parse_relaxed(self, stream):
        raise NotImplementedError()

class PhylipWriter(iosys.DataWriter):
    "Implements the DataWriter interface for writing PHYLIP files."

    def __init__(self, **kwargs):
        "Calls the base class constructor."
        iosys.DataWriter.__init__(self, **kwargs)
        self.strict = kwargs.get("strict", False)
        self.spaces_to_underscores = kwargs.get("spaces_to_underscores", False)
        self.force_unique_taxon_labels = kwargs.get("force_unique_taxon_labels", False)

    def get_taxon_label_map(self, dataset):
        if self.attached_taxon_set is not None:
            taxon_sets = [self.attached_taxon_set]
        else:
            taxon_sets = dataset.taxon_sets
        taxon_label_map = {}
        taxa = []
        if self.strict:
            max_label_len = STRICT_MODE_MAX_LABEL_LENGTH
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
        if self.strict:
            for t in taxon_label_map:
                label = taxon_label_map[t]
                if len(label) < STRICT_MODE_MAX_LABEL_LENGTH:
                    taxon_label_map[t] = label.ljust(STRICT_MODE_MAX_LABEL_LENGTH)
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
        self.force_unique_taxon_labels = kwargs.get("force_unique_taxon_labels", False)

        if self.exclude_chars:
            return self.dataset

        assert self.dataset is not None, \
            "PhylipWriter instance is not attached to a DataSet: no source of data"

        char_matrix = None
        if len(self.dataset.char_matrices) == 0:
            raise ValueError("No character data in DataSet")
        if self.attached_taxon_set is not None:
            for c in self.dataset.char_matrices:
                if c.taxon_set is self.attached_taxon_set:
                    char_matrix = c
                    break
            else:
                raise ValueError("No character matrix found associated with TaxonSet %s" % repr(self.attached_taxon_set))
        else:
            if len(char_matrix) > 1:
                raise ValueError('Multiple character matrices in DataSet, and Writer is not bound to a TaxonSet')
            char_matrix = self.dataset.char_matrices[0]
        assert char_matrix is not None

        if self.strict or self.force_unque_taxon_labels:
            taxon_label_map = self.get_taxon_label_map(self.dataset)
        else:
            taxon_label_map = {}
            for taxon_set in taxon_sets:
                for taxon in taxon_set:
                    taxon_label_map[taxon] = taxon.label
        maxlen = max([len(str(label)) for label in taxon_label_map.values()])

        char_matrix = self.dataset.char_matrices[0]
        n_seqs = len(char_matrix)
        n_sites = len(char_matrix.values()[0])
        stream.write("%d %d\n" % (n_seqs, n_sites))

        for taxon in char_matrix.taxon_set:
            stream.write("%s        %s\n" % ( str(taxon).ljust(maxlen), str(char_matrix[taxon].symbols_as_string()).replace(' ', '')))
