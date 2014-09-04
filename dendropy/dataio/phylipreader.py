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
Implementation of PHYLIP-format data reader.
"""


import re
from dendropy.dataio import ioservice
from dendropy.utility import filesys
from dendropy.utility import error

class PhylipReader(ioservice.DataReader):
    "Implements the DataReader interface for parsing PHYLIP files."

    # supported_data_types = ['dna', 'rna', 'protein', 'standard', 'restriction', 'infinite']
    # supported_matrix_types = [dataobject.DnaCharacterMatrix,
    #                           dataobject.RnaCharacterMatrix,
    #                           dataobject.ProteinCharacterMatrix,
    #                           dataobject.StandardCharacterMatrix,
    #                           dataobject.RestrictionSitesCharacterMatrix,
    #                           dataobject.InfiniteSitesCharacterMatrix]

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
        ioservice.DataReader.__init__(self)
        self.datatype_name = kwargs.pop("datatype_name", None)
        # if "char_matrix_type" in kwargs and "data_type" in kwargs:
        #     raise ValueError("Cannot specify both 'data_type' and 'char_matrix_type'")
        # if "data_type" in kwargs:
        #     data_type = kwargs["data_type"].lower()
        #     if data_type not in PhylipReader.supported_data_types:
        #         raise ValueError("'%s' is not a valid data type specification; must be one of: %s" \
        #             % (", ".join([("'" + d + "'") for d in PhylipReader.supported_data_types])))
        #     else:
        #         self.char_matrix_type = dataobject.character_data_type_label_map[data_type]
        # elif "char_matrix_type" in kwargs:
        #     self.char_matrix_type = kwargs.pop("char_matrix_type")
        # else:
        #     raise ValueError("Must specify 'data_type' for PHYLIP format, one of: %s" % (PhylipReader.supported_data_types))
        # if self.char_matrix_type not in PhylipReader.supported_matrix_types:
        #     raise ValueError("'%s' is not a supported data type for PhylipReader" % self.char_matrix_type.__name__)
        self.strict = kwargs.pop("strict", False)
        self.interleaved = kwargs.pop("interleaved", False)
        self.multispace_delimiter = kwargs.pop("multispace_delimiter", False)
        self.underscores_to_spaces = kwargs.pop("underscores_to_spaces", False)
        self.ignore_invalid_chars = kwargs.pop("ignore_invalid_chars", False)
        self.check_for_unused_keyword_arguments(kwargs)
        self.ntax = None
        self.nchar = None
        self.char_matrix = None
        self.taxon_namespace = None

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

    def reset(self):
        self.ntax = None
        self.nchar = None
        self.char_matrix = None
        self.taxon_namespace = None

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            state_alphabet_factory=None,
            global_annotations_target=None):
        self.reset()
        self.taxon_namespace = taxon_namespace_factory(label=None)
        if self.datatype_name is None:
            raise TypeError("Data type must be specified for this schema")
        self.char_matrix = char_matrix_factory(
                self.datatype_name,
                label=None,
                taxon_namespace=self.taxon_namespace)
        if self.datatype_name == "standard":
            state_alphabet = state_alphabet_factory(
                fundamental_states="0123456789",
                no_data_symbol="?",
                gap_symbol="-",
                case_sensitive=False)
            self.char_matrix.state_alphabets.append(state_alphabet)
        lines = filesys.get_lines(stream)
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
        product = self.Product(
                taxon_namespaces=None,
                tree_lists=None,
                char_matrices=[self.char_matrix])
        return product

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
        current_taxon = self.char_matrix.taxon_namespace.require_taxon(label=seq_label)
        if current_taxon not in self.char_matrix:
            self.char_matrix[current_taxon] = self.char_matrix.new_sequence(taxon=current_taxon)
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
                state = self.char_matrix.default_state_alphabet[c]
            except KeyError:
                if not self.ignore_invalid_chars:
                    raise self._data_parse_error("Invalid state symbol for taxon '%s': '%s'" % (current_taxon.label, c),
                            line_index=line_index)
            else:
                self.char_matrix[current_taxon].append(state)

    def _parse_sequential(self, lines, line_num_start=1):
        seq_labels = []
        current_taxon = None
        for line_index, line in enumerate(lines):
            line = line.rstrip()
            if line == '':
                continue
            if current_taxon is None:
                seq_label = None
                if line.strip() and len(self.char_matrix.taxon_namespace) >= self.ntax:
                    raise self._data_parse_error("Cannot add new sequence '%s': declared number of sequences (%d) already defined" \
                                % (line.strip(), len(self.char_matrix.taxon_namespace)), line_index=line_index)
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
                current_taxon = self.char_matrix.taxon_namespace[paged_row]
            else:
                current_taxon, line = self._parse_taxon_from_line(line, line_index)
                if len(self.char_matrix.taxon_namespace) == self.ntax:
                    paged = True
                    paged_row = -1
            self._parse_sequence_from_line(current_taxon, line, line_index)

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

