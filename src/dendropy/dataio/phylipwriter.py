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
Implementation of PHYLIP-format data writer.
"""

from dendropy.dataio import ioservice
from dendropy.utility import textprocessing

STRICT_MODE_MAX_LABEL_LENGTH = 10

class PhylipWriter(ioservice.DataWriter):
    "Implements the DataWriter interface for writing PHYLIP files."

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------

        strict : bool
            If |True|, use 'strict' format, i.e., taxon labels given in
            first 10 characters, followed by sequence starting at character 11.
            Default is |False|: use 'relaxed' format, with arbitrary-length
            taxon labels separated from sequences by two or more spaces.
        spaces_to_underscores : bool
            If |True|, all spaces will be converted to underscores. Default is
            |False|: spaces will be preserved.
        taxon_label_fn: function object
            If specified, then this function will be called everytime taxon
            label is required. It will be passed a |Taxon| object as an
            argument and should represent the string or string-like object that
            should serve as the label.
        force_unique_taxon_labels : bool
            If |True|, then taxon labels will be modified to avoid duplicate
            labels. Default is |False|: taxon labels will not be modified.
        suppress_missing_taxa : bool
            If |True|, then taxa with zero characters will not be printed
            Default is |False|: all taxa will be printed
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.
        """
        ioservice.DataWriter.__init__(self, **kwargs)
        self.strict = kwargs.pop("strict", False)
        self.spaces_to_underscores = kwargs.pop("spaces_to_underscores", False)
        self.force_unique_taxon_labels = kwargs.pop("force_unique_taxon_labels", False)
        self.suppress_missing_taxa = kwargs.pop("suppress_missing_taxa", False)
        self.taxon_label_fn = kwargs.pop("taxon_label_fn", None)
        if self.taxon_label_fn is None:
            self.taxon_label_fn = lambda taxon: taxon.label
        self.check_for_unused_keyword_arguments(kwargs)

    def _write(self,
            stream,
            taxon_namespaces=None,
            tree_lists=None,
            char_matrices=None,
            global_annotations_target=None):
        for char_matrix in char_matrices:
            if (self.attached_taxon_namespace is not None
                    and char_matrix.taxon_namespace is not self.attached_taxon_namespace):
                continue
            self._write_char_matrix(stream, char_matrix)

    def _write_char_matrix(self, stream, char_matrix):
        "Writes dataset to a full PHYLIP document."

        if self.strict or self.force_unique_taxon_labels:
            taxon_label_map = self.get_taxon_label_map(char_matrix.taxon_namespace)
            if not self.strict:
                spacer = "  "
            else:
                spacer = ""
        else:
            taxon_label_map = {}
            for taxon in char_matrix.taxon_namespace:
                label = self.taxon_label_fn(taxon)
                if self.spaces_to_underscores:
                    label = label.replace(' ', '_')
                taxon_label_map[taxon] = label
            spacer = "  "
        maxlen = max([len(str(label)) for label in taxon_label_map.values()])
        n_seqs = len(char_matrix)
        n_sites = char_matrix.max_sequence_size
        stream.write("%d %d\n" % (n_seqs, n_sites))
        for taxon in char_matrix.taxon_namespace:
            label = taxon_label_map[taxon]
            if taxon in char_matrix:
                seq_vec = char_matrix[taxon].symbols_as_string()
            else:
                seq_vec = ""
            if len(seq_vec) or (not self.suppress_missing_taxa):
                stream.write("%s%s%s\n" % ( label.ljust(maxlen), spacer, str(seq_vec)))

    def get_taxon_label_map(self, taxon_namespace):
        taxon_label_map = {}
        if self.strict:
            max_label_len = STRICT_MODE_MAX_LABEL_LENGTH
        else:
            max_label_len = 0
        for taxon in taxon_namespace:
            label = self.taxon_label_fn(taxon)
            if self.spaces_to_underscores:
                label = label.replace(' ', '_')
            if self.strict:
                label = label[:max_label_len]
            taxon_label_map[taxon] = label
        taxon_label_map = textprocessing.unique_taxon_label_map(taxon_namespace, taxon_label_map, max_label_len)
        if self.strict:
            for t in taxon_label_map:
                label = taxon_label_map[t]
                if len(label) < STRICT_MODE_MAX_LABEL_LENGTH:
                    taxon_label_map[t] = label.ljust(STRICT_MODE_MAX_LABEL_LENGTH)
        return taxon_label_map


