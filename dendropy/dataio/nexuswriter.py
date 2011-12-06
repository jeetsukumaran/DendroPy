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

            `dataset`
                Data to be written.
            `simple`
                If True, write in simple NEXUS format, i.e. in a single "DATA"
                block, instead of separate "TAXA" and "CHARACTER" blocks.
                Default is False.
            `suppress_taxa_block`
                If True, do not write a "TAXA" block. Default is False.
            `file_comments`
                List of lines of text to be added as comments to the file.
            `preamble_blocks`
                List of strings to be written before data (e.g., PAUP blocks
                suppressing warnings etc.). Default is [].
            `supplemental_blocks`
                List of strings to be written after data (e.g., PAUP blocks,
                MrBayes blocks etc.). Default is [].
            `suppress_leaf_taxon_labels`
                If True, then taxon labels will not be printed for leaves.
                Default is False.
            `suppress_leaf_node_labels`
                If False, then node labels (if available) will be printed
                for leaves. Defaults to True. Note that DendroPy distinguishes
                between *taxon* labels and *node* labels. In a typical NEWICK
                string, taxon labels are printed for leaf nodes, while leaf
                node labels are ignored (hence the default 'True' setting to
                suppress leaf node labels).
            `suppress_internal_taxon_labels`
                If True, then taxon labels will not be printed for internal
                nodes.  Default is False.
                NOTE: this replaces the `internal_labels` argument which has
                been deprecated.
            `suppress_internal_node_labels`
                If True, internal node labels will not be written. Default is
                False.
                NOTE: this replaces the `internal_labels` argument which has
                been deprecated.
            `suppress_rooting`
                If True, will not write rooting statement. Default is False.
                NOTE: this replaces the `write_rooting` argument which has been
                deprecated.
            `suppress_edge_lengths`
                If True, will not write edge lengths. Default is False.
                NOTE: this replaces the `edge_lengths` argument which has been
                deprecated.
            `unquoted_underscores`
                If True, labels with underscores will not be quoted, which will
                mean that they will be interpreted as spaces if read again
                ("soft" underscores).  If False, then labels with underscores
                will be quoted, resulting in "hard" underscores.  Default is
                False.
                NOTE: this replaces the `quote_underscores` argument which has
                been deprecated.
            `preserve_spaces`
                If True, spaces not mapped to underscores in labels (which
                means any labels containing spaces will have to be
                quoted). Default is False.
                False.
            `store_tree_weights`
                If True, tree weights are written. Default is False.
            `suppress_annotations`
                If False, will *not* write annotations as comments. Default is
                False if ``simple`` is False; True otherwise (note that, in
                contrast to this, the default for NEWICK formats, which
                is False).
            `annotations_as_nhx`
                If True, and if `suppress_annotations` is False, will write
                annotation as NHX statements. Default is False.
            `suppress_item_comments`
                If True, will write any additional comments. Default is
                False if ``simple`` is False; True otherwise (note that, in
                contrast to this, the default for NEWICK formats, which
                is False).
            `node_label_element_separator`
                If both `suppress_leaf_taxon_labels` and
                `suppress_leaf_node_labels` are False, then this will be the
                string used to join them. Defaults to ' '.
            `node_label_compose_func`
                If not None, should be a function that takes a Node object as
                an argument and returns the string to be used to represent the
                node in the tree statement. The return value
                from this function is used unconditionally to print a node
                representation in a tree statement, by-passing the default
                labelling function (and thus ignoring
                `suppress_leaf_taxon_labels`, `suppress_leaf_node_labels=True`,
                `suppress_internal_taxon_labels`, `suppress_internal_node_labels`,
                etc.). Defaults to None.
            `edge_label_compose_func`
                If not None, should be a function that takes an Edge object as
                an argument, and returns the string to be used to represent the
                edge length in the tree statement.

        Typically, these keywords would be passed to the `write_to_path()`,
        `write_to_stream` or `as_string` arguments, when 'nexus' is used as
        the schema::

            d.write_to_path('data.nex', 'nexus',
                    simple=False,
                    suppress_taxa_block=True,
                    file_comments=None,
                    preamble_blocks=[],
                    supplemental_blocks=[],
                    suppress_leaf_taxon_labels=False,
                    suppress_leaf_node_labels=True,
                    suppress_internal_taxon_labels=False,
                    suppress_internal_node_labels=False,
                    suppress_rooting=False,
                    suppress_edge_lengths=False,
                    unquoted_underscores=False,
                    preserve_spaces=False,
                    store_tree_weights=False,
                    suppress_annotations=False,
                    annotations_as_nhx=False,
                    suppress_item_comments=False,
                    node_label_element_separator=' ',
                    node_label_compose_func=None,
                    edge_label_compose_func=None)

        """
        iosys.DataWriter.__init__(self, **kwargs)
        self.simple = kwargs.get("simple", False)
        self.suppress_taxa_block = kwargs.get("suppress_taxa_block", False)
        self.suppress_taxa_block = kwargs.get("exclude_taxa", self.suppress_taxa_block) # legacy
        self.suppress_taxa_block = not kwargs.get("taxa_block", not self.suppress_taxa_block) # legacy
        self.is_write_block_titles = kwargs.get("block_titles", None)
        self.file_comments = kwargs.get("file_comments", [])
        self.file_comments = kwargs.get("comment", self.file_comments) # legacy
        self.preamble_blocks = kwargs.get("preamble_blocks", [])
        self.supplemental_blocks = kwargs.get("supplemental_blocks", [])

        self.suppress_leaf_taxon_labels = kwargs.get("suppress_leaf_taxon_labels", False)
        self.suppress_leaf_node_labels = kwargs.get("suppress_leaf_node_labels", True)
        self.suppress_internal_taxon_labels = kwargs.get("suppress_internal_taxon_labels", False)
        self.suppress_internal_taxon_labels = not kwargs.get("internal_labels", not self.suppress_internal_taxon_labels) # legacy
        self.suppress_internal_node_labels = kwargs.get("suppress_internal_node_labels", False)
        self.suppress_internal_node_labels = not kwargs.get("internal_labels", not self.suppress_internal_node_labels) # legacy

        self.suppress_rooting = kwargs.get("suppress_rooting", False)
        self.suppress_rooting = not kwargs.get("write_rooting", not self.suppress_rooting) # legacy

        self.suppress_edge_lengths = kwargs.get("suppress_edge_lengths", False)
        self.suppress_edge_lengths = not kwargs.get("edge_lengths", not self.suppress_edge_lengths) # legacy

        self.unquoted_underscores = kwargs.get('unquoted_underscores', False)
        self.unquoted_underscores = not kwargs.get('quote_underscores', not self.unquoted_underscores) # legacy

        self.preserve_spaces = kwargs.get("preserve_spaces", False)
        self.store_tree_weights = kwargs.get("store_tree_weights", False)
        self.suppress_annotations = kwargs.get("suppress_annotations", self.simple)
        self.suppress_annotations = not kwargs.get("annotations_as_comments", not self.suppress_annotations) # legacy
        self.annotations_as_nhx = kwargs.get("annotations_as_nhx", False)
        self.suppress_item_comments = kwargs.get("suppress_item_comments", self.simple)
        self.suppress_item_comments = not kwargs.get("write_item_comments", not self.suppress_item_comments)

        self.node_label_element_separator = kwargs.get("node_label_element_separator", ' ')
        self.node_label_compose_func = kwargs.get("node_label_compose_func", None)
        self.edge_label_compose_func = kwargs.get("edge_label_compose_func", None)

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
        if self.file_comments is not None:
            if isinstance(self.file_comments, list):
                for line in self.file_comments:
                    if line.strip().replace("\n", "").replace("\r", ""):
                        stream.write("[ %s ]\n" % line)
                    else:
                        stream.write("\n")
                stream.write("\n")
            else:
                stream.write("[ %s ]\n\n" % self.file_comments)
        if self.preamble_blocks:
            for block in self.preamble_blocks:
                stream.write(block)
                stream.write("\n")
            stream.write("\n")
        if (( (not self.exclude_chars) and self.dataset.char_matrices) \
                or ( (not self.exclude_trees) and self.dataset.tree_lists)) \
                and (not self.simple) \
                and (not self.suppress_taxa_block):
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
            return "TITLE %s" % textutils.escape_nexus_token(block.label, preserve_spaces=self.preserve_spaces, quote_underscores=not self.unquoted_underscores)
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
            block.append('        %s' % textutils.escape_nexus_token(taxon.label, preserve_spaces=self.preserve_spaces, quote_underscores=not self.unquoted_underscores))
        block.append('  ;')
        block.append('END;\n\n')
        stream.write('\n'.join(block))

    def write_trees_block(self, tree_list, stream):
        block = []
        newick_writer = newick.NewickWriter(
                suppress_rooting=self.suppress_rooting,
                suppress_edge_lengths=self.suppress_edge_lengths,
                unquoted_underscores=self.unquoted_underscores,
                preserve_spaces=self.preserve_spaces,
                store_tree_weights=self.store_tree_weights,
                suppress_annotations=self.suppress_annotations,
                annotations_as_nhx=self.annotations_as_nhx,
                suppress_item_comments=self.suppress_item_comments,
                suppress_leaf_taxon_labels=self.suppress_leaf_taxon_labels,
                suppress_leaf_node_labels=self.suppress_leaf_node_labels,
                suppress_internal_taxon_labels=self.suppress_internal_taxon_labels,
                suppress_internal_node_labels=self.suppress_internal_node_labels,
                node_label_element_separator=self.node_label_element_separator,
                node_label_compose_func=self.node_label_compose_func,
                edge_label_compose_func=self.edge_label_compose_func,
                )
        block.append('BEGIN TREES;')
        if self._link_blocks():
            title = self.compose_block_title(tree_list)
            if title:
                block.append('    %s;' % title)
            if tree_list.taxon_set.labels:
                block.append('    LINK TAXA = %s;' % textutils.escape_nexus_token(tree_list.taxon_set.label, preserve_spaces=self.preserve_spaces, quote_underscores=not self.unquoted_underscores))
        for treeidx, tree in enumerate(tree_list):
            if tree.label:
                tree_name = tree.label
            else:
                tree_name = str(treeidx)
            newick_str = newick_writer.compose_tree(tree)
            block.append('    TREE %s = %s' % (textutils.escape_nexus_token(tree_name, preserve_spaces=self.preserve_spaces, quote_underscores=not self.unquoted_underscores),
                newick_str))
        block.append('END;\n\n')
        stream.write('\n'.join(block))

    def write_char_block(self, char_matrix, stream):
        nexus = []
        taxlabels = [textutils.escape_nexus_token(taxon.label, preserve_spaces=self.preserve_spaces, quote_underscores=not self.unquoted_underscores) for taxon in char_matrix.taxon_set]
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
                nexus.append('    LINK TAXA = %s;' % textutils.escape_nexus_token(char_matrix.taxon_set.label, preserve_spaces=self.preserve_spaces, quote_underscores=not self.unquoted_underscores))
        nexus.append('    DIMENSIONS %s NCHAR=%d;' % (ntaxstr, nchar))
        nexus.append('    FORMAT %s;' % self.compose_format_terms(char_matrix))
        nexus.append('    MATRIX')
        state_string_map = {}
        if isinstance(char_matrix, dataobject.ContinuousCharacterMatrix):
            for taxon in char_matrix.taxon_set:
                seq = " ".join([str(v) for v in char_matrix[taxon]])
                nexus.append('%s    %s' % (textutils.escape_nexus_token(taxon.label, preserve_spaces=self.preserve_spaces, quote_underscores=not self.unquoted_underscores).ljust(max_label_len), seq))
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
                nexus.append('%s    %s' % (textutils.escape_nexus_token(taxon.label, preserve_spaces=self.preserve_spaces, quote_underscores=not self.unquoted_underscores).ljust(max_label_len), seq.getvalue()))
        nexus.append('    ;')
        nexus.append('END;\n\n')
        if hasattr(char_matrix, "character_subsets"):
            nexus.append('BEGIN SETS;')
            for label, char_set in char_matrix.character_subsets.items():
                label = textutils.escape_nexus_token(char_set.label,
                        preserve_spaces=self.preserve_spaces,
                        quote_underscores=not self.unquoted_underscores)
                pos = " ".join(str(c+1) for c in char_set.character_indices)
                nexus.append('    charset %s = %s;\n' % (label, pos))
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

