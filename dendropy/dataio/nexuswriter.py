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
Implementation of NEXUS-schema data writer.
"""

from cStringIO import StringIO
import re

from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)
from dendropy import dataobject
from dendropy.utility import textutils
from dendropy.utility import iosys
from dendropy.dataio import newick
from dendropy.dataio import nexustokenizer

###############################################################################
## NexusWriter

class NexusWriter(iosys.DataWriter):
    "Implements the DataWriter interface for handling NEXML files."

    def __init__(self, **kwargs):
        """
        __init__ recognizes the following keywords (in addition to those of `DataWriter.__init__`):

           - `simple` : if True, write in simple NEXUS format, i.e. in a
             single "DATA" block, instead of separate "TAXA" and "CHARACTER"
             blocks. [False]
           - `taxa_block` : if False, do not write a "TAXA" block [True]
           - `write_rooting` : if False, do not write a comment before each
             tree indicating its rooting state [True]
           - `store_tree_weights` : tree weights are store
           - `edge_lengths` : if False, edges will not write edge lengths [True]
           - `internal_labels` : if False, internal labels will not be written [True]
           - `annotations_as_comments` : if True, will write annotations as comments
           - `annotations_as_nhx` : if True, will write annotation as NHX statements
           - `write_item_comments` : if True, comments associated with (tree) nodes will be written [False]
           - `comment` : list of lines of text to be added as comments to the
             file
           - `supplemental_blocks` : list of strings to be written after data
           	 (e.g., PAUP blocks, MrBayes blocks etc.)
        """
        iosys.DataWriter.__init__(self, **kwargs)
        self.simple = kwargs.get("simple", False)
        self.exclude_taxa = kwargs.get("exclude_taxa", False)
        self.is_write_rooting = kwargs.get("write_rooting", True)
        self.store_tree_weights = kwargs.get('store_tree_weights', False)
        self.is_write_edge_lengths = kwargs.get("edge_lengths", True)
        self.is_write_internal_labels = kwargs.get("internal_labels", True)
        self.is_write_block_titles = kwargs.get("block_titles", None)
        self.preserve_spaces = kwargs.get("preserve_spaces", False)
        self.quote_underscores = kwargs.get('quote_underscores', True)
        self.annotations_as_comments = kwargs.get("annotations_as_comments", True)
        self.annotations_as_nhx = kwargs.get("annotations_as_nhx", False)
        self.nhx_key_to_func = kwargs.get("nhx_key_to_func_dict")
        self.is_write_item_comments = kwargs.get("write_item_comments", not self.simple)
        self.comment = kwargs.get("comment", [])
        self.supplemental_blocks = kwargs.get("supplemental_blocks", [])

    def write(self, stream):
        """
        Writes attached `DataSource` or `TaxonDomain` to the file-like object
        `stream`.
        """

        assert self.dataset is not None, \
            "NexusWriter instance is not attached to a DataSet: no source of data"

        if self.is_write_block_titles == False \
                and self.attached_taxon_set is None \
                and (len(self.dataset.taxon_sets) > 1):
            _LOG.warn("Multiple taxon sets in data, but directed not to write block titles: data file may not be interpretable")

        stream.write('#NEXUS\n\n')
        if self.comment is not None:
            if isinstance(self.comment, list):
                for line in self.comment:
                    if line.strip().replace("\n", "").replace("\r", ""):
                        stream.write("[ %s ]\n" % line)
                    else:
                        stream.write("\n")
                stream.write("\n")
            else:
                stream.write("[ %s ]\n\n" % self.comment)
        if (( (not self.exclude_chars) and self.dataset.char_matrices) \
                or ( (not self.exclude_trees) and self.dataset.tree_lists)) \
                and (not self.simple) \
                and (not self.exclude_taxa):
            for taxon_set in self.dataset.taxon_sets:
                if self.attached_taxon_set is None or taxon_set is self.attached_taxon_set:
                    self.write_taxa_block(taxon_set, stream=stream)
        if not self.exclude_chars:
            for char_matrix in self.dataset.char_matrices:
                if self.attached_taxon_set is None or char_matrix.taxon_set is self.attached_taxon_set:
                    self.write_char_block(char_matrix=char_matrix, stream=stream)
        if not self.exclude_trees:
            for tree_list in self.dataset.tree_lists:
                if self.attached_taxon_set is None or tree_list.taxon_set is self.attached_taxon_set:
                    self.write_trees_block(tree_list=tree_list, stream=stream)
        if self.supplemental_blocks:
            for block in self.supplemental_blocks:
                stream.write(block)
                stream.write("\n")

    def _link_blocks(self):
        """
        If only one taxon set in dataset, or in attached taxon set mode, then
        unless the 'block_titles' directive has been explicitly set to True
        by the user, block titles and links will not be written.
        """
        if self.is_write_block_titles is None:
            if self.attached_taxon_set is None and len(self.dataset.taxon_sets) > 1:
                return True
            else:
                return False
        else:
            return self.is_write_block_titles

    def compose_block_title(self, block):
        # if self.is_write_block_titles is False then no block titles;
        # if only one taxon set, or attached taxon set mode, unless self.is_write_block_titles
        # is explicitly True, then again, we do not write block titles
        if not self._link_blocks():
            return ""
        if not block.label:
            block.label = block.oid
        if block.label:
            return "TITLE %s" % textutils.escape_nexus_token(block.label, preserve_spaces=self.preserve_spaces, quote_underscores=self.quote_underscores)
        else:
            return ""

    def write_taxa_block(self, taxon_set, stream):
        block = []
        block.append('BEGIN TAXA;')
        if self._link_blocks():
            title = self.compose_block_title(taxon_set)
            if title:
                block.append('    %s;' % title)
        block.append('    DIMENSIONS NTAX=%d;' % len(taxon_set))
        block.append('    TAXLABELS')
        for taxon in taxon_set:
            block.append('        %s' % textutils.escape_nexus_token(taxon.label, preserve_spaces=self.preserve_spaces, quote_underscores=self.quote_underscores))
        block.append('  ;')
        block.append('END;\n\n')
        stream.write('\n'.join(block))

    def write_trees_block(self, tree_list, stream):
        block = []
        newick_writer = newick.NewickWriter(
                edge_lengths=self.is_write_edge_lengths,
                internal_labels=self.is_write_internal_labels,
                preserve_spaces=self.preserve_spaces,
                quote_underscores=self.quote_underscores,
                annotations_as_comments=self.annotations_as_comments,
                annotations_as_nhx=self.annotations_as_nhx,
                nhx_key_to_func_dict=self.nhx_key_to_func,
                write_item_comments=self.is_write_item_comments)
        block.append('BEGIN TREES;')
        if self._link_blocks():
            title = self.compose_block_title(tree_list)
            if title:
                block.append('    %s;' % title)
            if tree_list.taxon_set.labels:
                block.append('    LINK TAXA = %s;' % textutils.escape_nexus_token(tree_list.taxon_set.label, preserve_spaces=self.preserve_spaces, quote_underscores=self.quote_underscores))
        for treeidx, tree in enumerate(tree_list):
            if tree.label:
                tree_name = tree.label
            else:
                tree_name = str(treeidx)
            newick_str = newick_writer.compose_tree(tree)
            block.append('    TREE %s = %s' % (textutils.escape_nexus_token(tree_name, preserve_spaces=self.preserve_spaces, quote_underscores=self.quote_underscores),
                newick_str))
        block.append('END;\n\n')
        stream.write('\n'.join(block))

    def write_char_block(self, char_matrix, stream):
        nexus = []
        taxlabels = [textutils.escape_nexus_token(taxon.label, preserve_spaces=self.preserve_spaces, quote_underscores=self.quote_underscores) for taxon in char_matrix.taxon_set]
        max_label_len = max([len(label) for label in taxlabels])
        nchar = max([len(seq) for seq in char_matrix.values()])
        if self.simple:
            nexus.append('BEGIN DATA;')
            ntaxstr = "NTAX=%d" % len(taxlabels)
        else:
            nexus.append('BEGIN CHARACTERS;')
            ntaxstr = ""
        if self._link_blocks():
            title = self.compose_block_title(char_matrix)
            if title:
                nexus.append('    %s;' % title)
            if char_matrix.taxon_set.label:
                nexus.append('    LINK TAXA = %s;' % textutils.escape_nexus_token(char_matrix.taxon_set.label, preserve_spaces=self.preserve_spaces, quote_underscores=self.quote_underscores))
        nexus.append('    DIMENSIONS %s NCHAR=%d;' % (ntaxstr, nchar))
        nexus.append('    FORMAT %s;' % self.compose_format_terms(char_matrix))
        nexus.append('    MATRIX')
        state_string_map = {}
        if isinstance(char_matrix, dataobject.ContinuousCharacterMatrix):
            for taxon in char_matrix.taxon_set:
                seq = " ".join([str(v) for v in char_matrix[taxon]])
                nexus.append('%s    %s' % (textutils.escape_nexus_token(taxon.label, preserve_spaces=self.preserve_spaces, quote_underscores=self.quote_underscores).ljust(max_label_len), seq))
        else:
            for taxon in char_matrix.taxon_set:
                seq_vec = char_matrix[taxon]
                seq = StringIO()
                for cell in seq_vec:
                    state = cell.value
                    assert state is not None, "Undefined state encountered in character sequence."
                    try:
                        seq.write(state_string_map[state])
                    except:
                        if state.symbol is not None:
                            state_string_map[state] = state.symbol
                        elif state.multistate == dataobject.StateAlphabetElement.AMBIGUOUS_STATE:
                            state_string_map[state] = "{%s}" % ("".join(state.fundamental_symbols))
                        elif state.multistate == dataobject.StateAlphabetElement.POLYMORPHIC_STATE:
                            state_string_map[state] = "(%s)" % ("".join(state.fundamental_symbols))
                        else:
                            raise Exception("Could not match character state to symbol: '%s'." % state)
                        seq.write(state_string_map[state])
                nexus.append('%s    %s' % (textutils.escape_nexus_token(taxon.label, preserve_spaces=self.preserve_spaces, quote_underscores=self.quote_underscores).ljust(max_label_len), seq.getvalue()))
        nexus.append('    ;')
        nexus.append('END;\n\n')
        stream.write('\n'.join(nexus))

    def compose_format_terms(self, char_matrix):
        format = []
        if isinstance(char_matrix, dataobject.DnaCharacterMatrix):
            format.append("DATATYPE=DNA")
            format.append("GAP=- MISSING=? MATCHCHAR=.")
        elif isinstance(char_matrix, dataobject.RnaCharacterMatrix):
            format.append("DATATYPE=RNA")
            format.append("GAP=- MISSING=? MATCHCHAR=.")
        elif isinstance(char_matrix, dataobject.NucleotideCharacterMatrix):
            format.append("DATATYPE=NUCLEOTIDE")
            format.append("GAP=- MISSING=? MATCHCHAR=.")
        elif isinstance(char_matrix, dataobject.ProteinCharacterMatrix):
            format.append("DATATYPE=PROTEIN")
            format.append("GAP=- MISSING=? MATCHCHAR=.")
        elif isinstance(char_matrix, dataobject.ContinuousCharacterMatrix):
            format.append("DATATYPE=CONTINUOUS ITEMS=(STATES)")
        else:
            format.append("DATATYPE=STANDARD")

            fundamental_symbols = set()
            for state_alphabet in char_matrix.state_alphabets:
                for s in state_alphabet.fundamental_states():
                    if s.symbol is not None:
                        fundamental_symbols.add(s.symbol)
                    else:
                        raise Exception("Could not match character state to symbol: '%s'." % s)
            format.append('SYMBOLS="%s"' % "".join(fundamental_symbols))

            equates = set()
            for state_alphabet in char_matrix.state_alphabets:
                for a in state_alphabet.ambiguous_states():
                    if a.symbol == "?":
                        format.append("MISSING=?")
                    elif a.symbol == "-":
                        format.append("GAP=-")
                    else:
                        if a.symbol is not None:
                            equates.add("%s={%s}" % (a.symbol, "".join(a.fundamental_symbols)))

            for state_alphabet in char_matrix.state_alphabets:
                for p in state_alphabet.polymorphic_states():
                    if p.symbol is not None:
                        equates.add("%s=(%s)" % (p.symbol, "".join(p.fundamental_symbols)))

            if equates:
                format.append('EQUATE="%s"' % equates)

        return ' '.join(format)

