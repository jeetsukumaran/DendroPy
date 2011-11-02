#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
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
Implementation of PHYLIP-schema i/o client(s).
"""

from cStringIO import StringIO
import re
from dendropy.utility import iosys
from dendropy.utility import textutils
from dendropy.utility import fileutils
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

    class PhylipStrictSequentialError(error.DataParseError):
        def __init__(self, *args, **kwargs):
            error.DataParseError.__init__(self, *args, **kwargs)

    class PhylipStrictInterleavedError(error.DataParseError):
        def __init__(self, *args, **kwargs):
            error.DataParseError.__init__(self, *args, **kwargs)

    class PhylipRelaxedSequentialError(error.DataParseError):
        def __init__(self, *args, **kwargs):
            error.DataParseError.__init__(self, *args, **kwargs)

    class PhylipRelaxedInterleavedError(error.DataParseError):
        def __init__(self, *args, **kwargs):
            error.DataParseError.__init__(self, *args, **kwargs)

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
        elif "char_matrix_type" in kwargs:
            self.char_matrix_type = kwargs.get("char_matrix_type")
        else:
            raise ValueError("Must specify 'data_type' for PHYLIP format, one of: %s" % (PhylipReader.supported_data_types))
        if self.char_matrix_type not in PhylipReader.supported_matrix_types:
            raise ValueError("'%s' is not a supported data type for PhylipReader" % self.char_matrix_type.__name__)
        self.strict = kwargs.get("strict", False)
        self.interleaved = kwargs.get("interleaved", False)
        self.multispace_delimiter = kwargs.get("multispace_delimiter", False)
        self.underscores_to_spaces = kwargs.get("underscores_to_spaces", False)
        self.ignore_invalid_chars = kwargs.get("ignore_invalid_chars", False)
        self.ntax = None
        self.nchar = None

    def describe_mode(self):
        parts = []
        if self.strict:
            parts.append("strict")
        else:
            parts.append("relaxed")
        if self.interleaved:
            parts.append("interleaved")
        else:
            parts.append("sequential")
        return ", ".join(parts)

    def read(self, stream):
        if self.exclude_chars:
            return self.dataset
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        taxon_set = self.get_default_taxon_set()
        self.char_matrix = self.dataset.new_char_matrix(char_matrix_type=self.char_matrix_type,
                taxon_set=taxon_set)
        if isinstance(self.char_matrix, dataobject.StandardCharacterMatrix) \
            and len(self.char_matrix.state_alphabets) == 0:
                self.char_matrix.state_alphabets.append(dataobject.get_state_alphabet_from_symbols("0123456789"))
                self.char_matrix.default_state_alphabet = self.char_matrix.state_alphabets[0]
        if self.char_matrix.default_state_alphabet is not None:
            self.symbol_state_map = self.char_matrix.default_state_alphabet.symbol_state_map()
        elif len(self.char_matrix.state_alphabets) == 0:
            raise ValueError("No state alphabets defined")
        elif len(self.char_matrix.state_alphabets) > 1:
            raise NotImplementedError("Mixed state-alphabet matrices not supported")
        else:
            self.symbol_state_map = self.char_matrix.state_alphabets[0]

        self.stream = stream
        lines = fileutils.get_lines(self.stream)
        if len(lines) == 0:
            raise error.DataSourceError("No data in source", stream=self.stream)
        elif len(lines) <= 2:
            raise error.DataParseError("Expecting at least 2 lines in PHYLIP format data source", stream=self.stream)

        desc_line = lines[0]
        lines = lines[1:]
        m = re.match('\s*(\d+)\s+(\d+)\s*$', desc_line)
        if m is None:
            raise self._data_parse_error("Invalid data description line: '%s'" % desc_line)
        self.ntax = int(m.groups()[0])
        self.nchar = int(m.groups()[1])
        if self.ntax == 0 or self.nchar == 0:
            raise error.DataSourceError("No data in source", stream=self.stream)
        if self.interleaved:
            self._parse_interleaved(lines)
        else:
            self._parse_sequential(lines)
        self.stream = None
        return self.dataset

    def _data_parse_error(self, message, line_index=None):
        if line_index is None:
            row = None
        else:
            row = line_index + 2
        if self.strict and self.interleaved:
            error_type = PhylipReader.PhylipStrictInterleavedError
        elif self.strict:
            error_type = PhylipReader.PhylipStrictSequentialError
        elif self.interleaved:
            error_type = PhylipReader.PhylipRelaxedInterleavedError
        else:
            error_type = PhylipReader.PhylipStrictSequentialError
        return error_type(message, row=row, stream=self.stream)

    def _parse_taxon_from_line(self, line, line_index):
        if self.strict:
            seq_label = line[:10].strip()
            line = line[10:]
        else:
            if self.multispace_delimiter:
                parts = re.split('[ \t]{2,}', line, maxsplit=1)
            else:
                parts = re.split('[ \t]{1,}', line, maxsplit=1)
            seq_label = parts[0]
            if len(parts) < 2:
                line = ''
            else:
                line = parts[1]
        seq_label = seq_label.strip()
        if not seq_label:
            raise self._data_parse_error("Expecting taxon label", line_index=line_index)
        if self.underscores_to_spaces:
            seq_label = seq_label.replace('_', ' ')
        current_taxon = self.char_matrix.taxon_set.require_taxon(label=seq_label)
        if current_taxon not in self.char_matrix:
            self.char_matrix[current_taxon] = dataobject.CharacterDataVector(taxon=current_taxon)
        else:
            if len(self.char_matrix[current_taxon]) >= self.nchar:
                raise self._data_parse_error("Cannot add characters to sequence for taxon '%s': already has declared number of characters (%d)" \
                        % (current_taxon.label, self.char_matrix[current_taxon]), line_index=line_index)
        return current_taxon, line

    def _parse_sequence_from_line(self, current_taxon, line, line_index):
        for c in line:
            if c in [' ', '\t']:
                continue
            try:
                state = self.symbol_state_map[c.upper()]
            except KeyError:
                if not self.ignore_invalid_chars:
                    raise self._data_parse_error("Invalid state symbol for taxon '%s': '%s'" % (current_taxon.label, c),
                            line_index=line_index)
            else:
                self.char_matrix[current_taxon].append(dataobject.CharacterDataCell(value=state))

    def _parse_sequential(self, lines, line_num_start=1):
        seq_labels = []
        current_taxon = None
        for line_index, line in enumerate(lines):
            line = line.rstrip()
            if line == '':
                continue
            if current_taxon is None:
                seq_label = None
                if line.strip() and len(self.char_matrix.taxon_set) >= self.ntax:
                    raise self._data_parse_error("Cannot add new sequence '%s': declared number of sequences (%d) already defined" \
                                % (line.strip(), len(self.char_matrix.taxon_set)), line_index=line_index)
                current_taxon, line = self._parse_taxon_from_line(line, line_index)
            self._parse_sequence_from_line(current_taxon, line, line_index)
            if len(self.char_matrix[current_taxon]) >= self.nchar:
                current_taxon = None

    def _parse_interleaved(self, lines, line_num_start=1):
        seq_labels = []
        current_taxon = None
        paged = False
        paged_row = -1
        for line_index, line in enumerate(lines):
            current_taxon = None
            line = line.rstrip()
            if line == '':
                continue
            paged_row += 1
            if paged_row >= self.ntax:
                paged_row = 0
            if paged:
                current_taxon = self.char_matrix.taxon_set[paged_row]
            else:
                current_taxon, line = self._parse_taxon_from_line(line, line_index)
                if len(self.char_matrix.taxon_set) == self.ntax:
                    paged = True
                    paged_row = -1
            self._parse_sequence_from_line(current_taxon, line, line_index)


class PhylipWriter(iosys.DataWriter):
    "Implements the DataWriter interface for writing PHYLIP files."

    def __init__(self, **kwargs):
        """
        __init__ recognizes the following keywords (in addition to those of `DataWriter.__init__`):

            - `strict` (boolean)
            - `spaces_to_underscores` (boolean)
            - `force_unique_taxon_labels` (boolean)
        """
        iosys.DataWriter.__init__(self, **kwargs)
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
            stream.write("%s%s%s\n" % ( label.ljust(maxlen), spacer, str(char_matrix[taxon].symbols_as_string()).replace(' ', '')))
