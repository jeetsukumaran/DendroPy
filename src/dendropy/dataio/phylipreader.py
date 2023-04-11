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
        """
        Keyword Arguments
        -----------------
        data_type: str
            When reading into a |DataSet| object, the type of data must be
            specified: "dna", "rna", "protein", "restriction", "infinite",
            "standard", or "continuous".
        default_state_alphabet: |StateAlphabet| instance
            A |StateAlphabet| object to be used to manage the alphabet of the
            characters (|StandardCharacterMatrix| **only**).
        strict : bool
            If |True|, then data is given in 'strict' format, where first 10
            characters are the taxon label and remaining characters are the sequence.
            Default is |False|: relaxed format, where taxon labels are of
            arbitrary length and separation of sequences are is by one or more (if
            ``multispace_delimiter`` is |False|) or two or more (if
            ``multispace_delimiter`` is |True|) spaces.
        interleaved : bool
            If |True|, then data is in interleaved format.
            Default is |False|: data is non-interleaved.
        multispace_delimiter: bool
            If |True| (and ``strict`` is |False|), then at least two spaces are
            required to delimit taxon label and associated sequence. Default is
            |False|: one or more spaces delimit taxon label and associated
            sequence.
        underscore_to_spaces: bool
            If |True|, then underscores in taxon labels are converted to
            spaces. Default is |False|: underscores are not converted.
        ignore_invalid_chars : bool
            If |True| then any invalid characters in sequences will be ignored.
            Default is |False|: invalid characters result in errors.
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.
        """
        ioservice.DataReader.__init__(self)
        self.data_type = kwargs.pop("data_type", None)
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
        self.default_state_alphabet = kwargs.pop("default_state_alphabet", None)
        if self.default_state_alphabet is not None:
            if self.data_type is None:
                self.data_type = "standard"
            elif self.data_type != "standard":
                raise ValueError("Cannot specify 'default_state_alphabet' with data type of '{}'".format(self.data_type))
        self.check_for_unused_keyword_arguments(kwargs)
        self.ntax = None
        self.taxa_processed = None
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
        self.stream = None
        self.taxa_processed = set()

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            state_alphabet_factory=None,
            global_annotations_target=None):
        self.reset()
        self.stream = stream
        self.taxon_namespace = taxon_namespace_factory(label=None)
        if self.data_type is None:
            raise TypeError("Data type must be specified for this schema")
        if self.data_type == "standard" and self.default_state_alphabet is not None:
            self.char_matrix = char_matrix_factory(
                    self.data_type,
                    label=None,
                    taxon_namespace=self.taxon_namespace,
                    default_state_alphabet=self.default_state_alphabet,
                    )
        else:
            self.char_matrix = char_matrix_factory(
                    self.data_type,
                    label=None,
                    taxon_namespace=self.taxon_namespace)
            if self.data_type == "standard":
                state_alphabet = state_alphabet_factory(
                    fundamental_states="0123456789",
                    no_data_symbol="?",
                    gap_symbol="-",
                    case_sensitive=False)
                self.char_matrix.state_alphabets.append(state_alphabet)
        lines = filesys.get_lines(stream)
        if len(lines) == 0:
            raise error.DataParseError("No data in source", stream=self.stream)
        elif len(lines) <= 2:
            raise error.DataParseError("Expecting at least 2 lines in PHYLIP format data source", stream=self.stream)
        desc_line = lines[0]
        lines = lines[1:]
        m = re.match(r'\s*(\d+)\s+(\d+)\s*$', desc_line)
        if m is None:
            raise self._data_parse_error("Invalid data description line: '%s'" % desc_line)
        self.ntax = int(m.groups()[0])
        self.nchar = int(m.groups()[1])
        if self.ntax == 0 or self.nchar == 0:
            raise error.DataParseError("No data in source", stream=self.stream)
        if self.interleaved:
            self._parse_interleaved(lines)
        else:
            self._parse_sequential(lines)
        if len(self.taxa_processed) != self.ntax:
            self._taxon_error(num_expected=self.ntax, found=self.taxa_processed)
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
        self.taxa_processed.add(current_taxon)
        if len(self.taxa_processed) > self.ntax:
            self._taxon_error(num_expected=self.ntax, found=self.taxa_processed)
        return current_taxon, line

    def _parse_sequence_from_line(self, current_taxon, line, line_index):
        if self.data_type == "continuous":
            for c in line.split():
                if not c:
                    continue
                try:
                    state = float(c)
                except ValueError:
                    if not self.ignore_invalid_chars:
                        raise self._data_parse_error("Invalid state for taxon '%s': '%s'" % (current_taxon.label, c),
                                line_index=line_index)
                else:
                    self.char_matrix[current_taxon].append(state)
        else:
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
                current_taxon, line = self._parse_taxon_from_line(line, line_index)
                # if current_taxon not in self.char_matrix and len(self.char_matrix.taxon_namespace) >= self.ntax:
                #     raise self._data_parse_error("Cannot add new sequence %s: declared number of sequences (%d) already defined" \
                #                 % (current_taxon, len(self.char_matrix.taxon_namespace)), line_index=line_index)
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
        return error_type(message, line_num=row, stream=self.stream)

    def _taxon_error(self, num_expected, found):
        if num_expected == 1:
            n1 = "taxon"
        else:
            n1 = "taxa"
        if len(found) == 1:
            n2 = "taxon"
        else:
            n2 = "taxa"
        if num_expected > len(found):
            a = "only "
        else:
            a = ""
        raise error.DataParseError("{} {} expected but {}{} {} found: {}".format(
            num_expected,
            n1,
            a,
            len(found),
            n2,
            ", ".join("{}".format(t) for t in found)))
