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
Implementation of FASTA-format data reader and writer.
"""

from cStringIO import StringIO
import re
import textwrap

from dendropy import dataobject
from dendropy.utility import iosys
from dendropy.utility.error import DataParseError

class FastaReader(iosys.DataReader):
    "Encapsulates loading and parsing of a FASTA format file."

    supported_data_types = ['dna', 'rna', 'protein', 'standard', 'restriction', 'infinite']
    supported_matrix_types = [dataobject.DnaCharacterMatrix,
                              dataobject.RnaCharacterMatrix,
                              dataobject.ProteinCharacterMatrix,
                              dataobject.StandardCharacterMatrix,
                              dataobject.RestrictionSitesCharacterMatrix,
                              dataobject.InfiniteSitesCharacterMatrix]

    def __init__(self, **kwargs):
        """
        __init__ requires either `char_matrix_type` or `data_type` kwargs.

        The supported kwargs are:

            - `row_type` can be RICH or STR,
            - `char_matrix_type` should be one of the `CharacterMatrix` types.
            - `data_type` (should be in FastaReader.supported_data_types)
        """
        iosys.DataReader.__init__(self, **kwargs)
        self.simple_rows = kwargs.get('row_type', 'rich').upper() == 'STR'

        self.char_matrix_type = kwargs.get("char_matrix_type")
        data_type = kwargs.get("data_type", '').lower()
        if data_type:
            if self.char_matrix_type is not None:
                raise ValueError("Cannot specify both 'data_type' and 'char_matrix_type'")
            if data_type not in FastaReader.supported_data_types:
                raise ValueError("'%s' is not a valid data type specification; must be one of: %s" \
                    % (", ".join([("'" + d + "'") for d in FastaReader.supported_data_types])))
            else:
                self.char_matrix_type = dataobject.character_data_type_label_map[data_type]
        elif self.char_matrix_type is None:
            raise ValueError("Must specify 'data_type' for FASTA format, one of: %s" % (FastaReader.supported_data_types))
        if self.char_matrix_type not in FastaReader.supported_matrix_types:
            raise ValueError("'%s' is not a supported data type for FastaReader" % self.char_matrix_type.__name__)

    def read(self, stream):
        """
        Main file parsing driver.
        """

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

        curr_vec = None
        curr_taxon = None

        if self.simple_rows:
            legal_chars = self.char_matrix.default_state_alphabet.get_legal_symbols_as_str()

        for line_index, line in enumerate(stream):
            s = line.strip()
            if not s:
                continue
            if s.startswith('>'):
                if self.simple_rows and curr_taxon and curr_vec:
                    self.char_matrix[curr_taxon] = "".join(curr_vec)
                name = s[1:].strip()
                curr_taxon = taxon_set.require_taxon(label=name)
                if curr_taxon in self.char_matrix:
                    raise DataParseError(message="Fasta error: Repeated sequence name (%s) found" % name, row=line_index + 1, stream=stream)
                if curr_vec is not None and len(curr_vec) == 0:
                    raise DataParseError(message="Fasta error: Expected sequence, but found another sequence name (%s)" % name, row=line_index + 1, stream=stream)
                if self.simple_rows:
                    curr_vec = []
                else:
                    curr_vec = dataobject.CharacterDataVector(taxon=curr_taxon)
                    self.char_matrix[curr_taxon] = curr_vec
            elif curr_vec is None:
                raise DataParseError(message="Fasta error: Expecting a lines starting with > before sequences", row=line_index + 1, stream=stream)
            else:
                if self.simple_rows:
                    for col_ind, c in enumerate(s):
                        c = c.strip()
                        if not c:
                            continue
                        if c not in legal_chars:
                            DataParseError(message='Unrecognized sequence symbol "%s"' % c, row=line_index + 1, column=col_ind + 1, stream=stream)
                        curr_vec.append(c)
                else:
                    for col_ind, c in enumerate(s):
                        c = c.strip()
                        if not c:
                            continue
                        try:
                            state = self.symbol_state_map[c]
                            curr_vec.append(dataobject.CharacterDataCell(value=state))
                        except:
                            raise DataParseError(message='Unrecognized sequence symbol "%s"' % c, row=line_index + 1, column=col_ind + 1, stream=stream)
        if self.simple_rows and curr_taxon and curr_vec:
            self.char_matrix[curr_taxon] = "".join(curr_vec)
        return self.dataset

class DNAFastaReader(FastaReader):
    def __init__(self, **kwargs):
        FastaReader.__init__(self, char_matrix_type=dataobject.DnaCharacterMatrix, **kwargs)

class RNAFastaReader(FastaReader):
    def __init__(self, **kwargs):
        FastaReader.__init__(self, char_matrix_type=dataobject.RnaCharacterMatrix, **kwargs)

class ProteinFastaReader(FastaReader):
    def __init__(self, **kwargs):
        FastaReader.__init__(self, char_matrix_type=dataobject.ProteinCharacterMatrix, **kwargs)

class FastaWriter(iosys.DataWriter):
    """
    Implements the DataWriter interface for handling FASTA files.
    Additional keyword arguments:

        -`wrap`: if True, wraps text; slows things down
          quite a bit, so default is no wrap.
    """

    def __init__(self, **kwargs):
        iosys.DataWriter.__init__(self, **kwargs)
        self.wrap = kwargs.get("wrap", False)
        self.wrap_width = kwargs.get("wrap_width", 70)

    def write(self, stream):
        """
        Writes attached DataSet in FASTA format to a
        file-like object `stream`.
        """
        assert self.dataset is not None, \
            "FastaWriter instance is not attached to a DataSet: no source of data"
        if self.exclude_chars:
            return

        try:
            tw = textwrap.TextWrapper(width=self.wrap_width,
                    break_long_words=True,
                    break_on_hyphens=False)
        except TypeError:
            # Python 2.5 does not support break_on_hyphens
            tw = textwrap.TextWrapper(width=self.wrap_width,
                    break_long_words=True)

        for char_matrix in self.dataset.char_matrices:
            if self.attached_taxon_set is not None \
                    and char_matrix.taxon_set is not self.attached_taxon_set:
                continue
            for taxon in char_matrix.taxon_set:
                try:
                    seqs = char_matrix[taxon]
                    stream.write(">%s\n" % taxon.label)
                    if isinstance(seqs, dataobject.CharacterDataVector):
                        seqs = seqs.symbols_as_string()
                    if self.wrap:
                        seqs = tw.fill(seqs)
                    stream.write("%s\n\n" % seqs)
                except KeyError:
                    # if no sequences found associated with this taxa: skip
                    pass


