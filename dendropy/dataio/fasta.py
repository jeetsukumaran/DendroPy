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
Implementation of FASTA-format data reader and writer.
"""

from cStringIO import StringIO
import re
import textwrap

from dendropy import dataobject
from dendropy.utility import iosys
from dendropy.utility.error import DataFormatError

class FastaReader(iosys.DataReader):
    "Encapsulates loading and parsing of a FASTA format file."

    def __init__(self, **kwargs):
        """
        Keywords `row_type` kwarg can be RICH or STR, `char_matrix_type`
        is one of the `CharacterMatrix` types.
        """
        iosys.DataReader.__init__(self, **kwargs)
        self.char_matrix_type = kwargs.get("char_matrix_type", dataobject.DnaCharacterMatrix)

    def read(self, stream, **kwargs):
        """
        Main file parsing driver.
        """
        simple_rows = kwargs.get('row_type', 'rich').upper() == 'STR'
        if self.dataset is None:
            self.dataset = dataobject.DataSet()
        taxon_set = self.get_default_taxon_set(**kwargs)
        char_matrix = self.dataset.new_char_matrix(char_matrix_type=self.char_matrix_type, taxon_set=taxon_set)
        char_matrix.taxon_set = taxon_set
        symbol_state_map = char_matrix.default_state_alphabet.symbol_state_map()

        curr_vec = None
        curr_taxon = None

        if simple_rows:
            legal_chars = char_matrix.default_state_alphabet.get_legal_symbols_as_str()

        for line_index, line in enumerate(stream):
            s = line.strip()
            if not s:
                continue
            if s.startswith('>'):
                if simple_rows and curr_taxon and curr_vec:
                    char_matrix[curr_taxon] = "".join(curr_vec)
                name = s[1:].strip()
                curr_taxon = taxon_set.require_taxon(label=name)
                if curr_taxon in char_matrix:
                    raise DataFormatError(row=line_index + 1, message="Fasta error: Repeated sequence name (%s) found" % name)
                if curr_vec is not None and len(curr_vec) == 0:
                    raise DataFormatError(row=line_index + 1, message="Fasta error: Expected sequence, but found another sequence name (%s)" % name)
                if simple_rows:
                    curr_vec = []
                else:
                    curr_vec = dataobject.CharacterDataVector(taxon=curr_taxon)
                    char_matrix[curr_taxon] = curr_vec
            elif curr_vec is None:
                raise DataFormatError(row=line_index + 1, message="Fasta error: Expecting a lines starting with > before sequences")
            else:
                if simple_rows:
                    for col_ind, c in enumerate(s):
                        c = c.strip()
                        if not c:
                            continue
                        if c not in legal_chars:
                            DataFormatError(row=line_index + 1, column=col_ind + 1, message='Unrecognized sequence symbol "%s"' % c)
                        curr_vec.append(c)
                else:
                    for col_ind, c in enumerate(s):
                        c = c.strip()
                        if not c:
                            continue
                        try:
                            state = symbol_state_map[c]
                            curr_vec.append(dataobject.CharacterDataCell(value=state))
                        except:
                            raise DataFormatError(row=line_index + 1, column=col_ind + 1, message='Unrecognized sequence symbol "%s"' % c)
        if simple_rows and curr_taxon and curr_vec:
            char_matrix[curr_taxon] = "".join(curr_vec)
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

        -`wrap`: if non-zero, wraps text to this width; slows things down
            quite a bit, so default is no wrap (=0).
    """

    def __init__(self, **kwargs):
        "Calls the base class constructor."
        iosys.DataWriter.__init__(self, **kwargs)
        self.wrap = kwargs.get("wrap", 0)

    def write(self, stream, **kwargs):
        """
        Writes attached `DataSource` or `TaxonDomain` in FASTA format to a
        file-like object `stream`.
        """

        assert self.dataset is not None, \
            "FastaWriter instance is not attached to a DataSet: no source of data"
        if self.exclude_chars:
            return

        for char_matrix in self.dataset.char_matrices:
            for taxon in char_matrix.taxon_set:
                stream.write(">%s\n" % taxon.label)
                seqs = char_matrix[taxon]
                if isinstance(seqs, dataobject.CharacterDataVector):
                    seqs = seqs.symbols_as_string()
                if self.wrap > 0:
                    seqs = textwrap.fill(seqs, width=self.wrap, break_long_words=True)
                stream.write("%s\n\n" % seqs)


