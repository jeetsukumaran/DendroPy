# !/usr/bin/env python

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
Implementation of FASTA-format data writer.
"""

from dendropy.dataio import ioservice

class FastaWriter(ioservice.DataWriter):
    """
    Formatter for FASTA writer
    """

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------

        wrap: boolean, default: |True|
            If |False|, then sequences are written out as single, unbroken lines.
            Defaults to |True|: wraps sequences at 70 colums.
        """
        ioservice.DataWriter.__init__(self)
        self.wrap = kwargs.get("wrap", True)
        self.wrap_width = kwargs.get("wrap_width", 70)

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
        for taxon in char_matrix:
            stream.write(">{}\n".format(taxon.label))
            seq = char_matrix[taxon]
            if self.wrap:
                col_count = 0
                for c in seq:
                    if col_count == self.wrap_width:
                        stream.write("\n")
                        col_count = 0
                    stream.write(str(c))
                    col_count += 1
            else:
                s = "".join("{}".format(c) for c in seq)
                stream.write("{}\n".format(s))
            stream.write("\n\n")



