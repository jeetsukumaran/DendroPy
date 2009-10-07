#! /usr/bin/env python

############################################################################
##  fasta.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
############################################################################

"""
This module wraps routines needed for reading and writing data in
FASTA format.

*** VERY CRUDE, AND CURRENTLY ONLY SUPPORTS WRITING  ***
"""
import textwrap
from dendropy import get_logger, SyntaxException
from dendropy import datasets
from dendropy import characters
from dendropy.taxa import TaxaBlock

class FastaWriter(datasets.Writer):
    "Implements the DataWriter interface for handling FASTA files."
    
    def __init__(self):
        "Calls the base class constructor."
        datasets.Writer.__init__(self)
        
    def write_dataset(self, dataset, dest, store_chars=True, store_trees=True):
        "Writes dataset to a full PHYLIP document."
        if not store_chars:
            return
        for char_block in dataset.char_blocks:
            for taxon in char_block.matrix:
                dest.write(">%s\n" % str(taxon))
                seqs = str(char_block.matrix[taxon])
                seqsf = textwrap.fill(seqs, width=79, break_long_words=True, break_on_hyphens=False)
                dest.write("%s\n\n" % seqsf)

class FastaReader(datasets.Reader):
    "Encapsulates loading and parsing of a FASTA format file."

    def __init__(self, character_block_type):
        datasets.Reader.__init__(self)
        self.char_block_type = character_block_type

    def read_dataset(self, file_obj, dataset=None, **kwargs):
        """Main file parsing driver.
        row_type kwarg can be RICH or STR
        """
        if dataset is None:
            dataset = datasets.Dataset()

        simple_rows = kwargs.get('row_type', 'rich').upper() == 'STR'
        
        ntb = len(dataset.taxa_blocks)
        if ntb > 1:
            raise ValueError('datasets with multiple taxa blocks are not supported by FASTA')
        if ntb == 0:
            taxa_block = TaxaBlock()
            dataset.taxa_blocks.append(taxa_block)
        else:
            taxa_block = dataset.taxa_block[0]

        char_block = self.char_block_type()
        char_block.taxa_block = taxa_block
        dataset.add_char_block(char_block=char_block)

        symbol_state_map = char_block.default_state_alphabet.symbol_state_map()

        curr_vec = None
        curr_taxon = None
        
        if simple_rows:
            legal_chars = char_block.default_state_alphabet.get_legal_symbols_as_str()

        for line_index, line in enumerate(file_obj):
            s = line.strip()
            if not s:
                continue
            if s.startswith('>'):
                if simple_rows and curr_taxon and curr_vec:
                    char_block[curr_taxon] = "".join(curr_vec)
                name = s[1:].strip()
                curr_taxon = taxa_block.get_taxon(label=name)
                if curr_taxon in char_block:
                    raise SyntaxException(row=line_index + 1, message="Fasta error: Repeated sequence name (%s) found" % name)
                if curr_vec is not None and len(curr_vec) == 0:
                    raise SyntaxException(row=line_index + 1, message="Fasta error: Expected sequence, but found another sequence name (%s)" % name)
                if simple_rows:
                    curr_vec = []
                else:
                    curr_vec = characters.CharacterDataVector(taxon=curr_taxon)
                    char_block[curr_taxon] = curr_vec
            elif curr_vec is None:
                raise SyntaxException(row=line_index + 1, message="Fasta error: Expecting a lines starting with > before sequences")
            else:
                if simple_rows:
                    for col_ind, c in enumerate(s):
                        c = c.strip()
                        if not c:
                            continue
                        if c not in legal_chars:
                            SyntaxException(row=line_index + 1, column=col_ind + 1, message='Unrecognized sequence symbol "%s"' % c)
                        curr_vec.append(c)
                else:
                    for col_ind, c in enumerate(s):
                        c = c.strip()
                        if not c:
                            continue
                        try:
                            state = symbol_state_map[c]
                            curr_vec.append(characters.CharacterDataCell(value=state))
                        except:
                            raise SyntaxException(row=line_index + 1, column=col_ind + 1, message='Unrecognized sequence symbol "%s"' % c)
        if simple_rows and curr_taxon and curr_vec:
            char_block[curr_taxon] = "".join(curr_vec)
                

class DNAFastaReader(FastaReader):
    def __init__(self):
        FastaReader.__init__(self, character_block_type=characters.DnaCharactersBlock)

class RNAFastaReader(FastaReader):
    def __init__(self):
        FastaReader.__init__(self, character_block_type=characters.RnaCharactersBlock)

class ProteinFastaReader(FastaReader):
    def __init__(self):
        FastaReader.__init__(self, character_block_type=characters.ProteinCharactersBlock)
