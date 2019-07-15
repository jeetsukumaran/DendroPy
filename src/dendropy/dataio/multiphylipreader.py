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
Implementation of multiple-alignment PHYLIP-format data reader.
"""

import re
from dendropy.dataio import ioservice
from dendropy.dataio import phylipreader
from dendropy.utility import filesys
from dendropy.utility import error
from dendropy.utility.textprocessing import StringIO

class MultiPhylipReader(ioservice.DataReader):

    data_block_start_pattern = re.compile(r'\n*(?=\s*\d+\s+\d+\s*\n)', re.MULTILINE)

    def __init__(self, **kwargs):
        """
        See ``dendropy.dataio.phylipreader.PhylipReader`` for keyword arguments.
        """
        self._phylip_reader = phylipreader.PhylipReader(**kwargs)

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            state_alphabet_factory=None,
            global_annotations_target=None):
        lines = stream.read()
        blocks = MultiPhylipReader.data_block_start_pattern.split(lines)
        if len(blocks) == 1:
            raise error.DataParseError("No PHYLIP data blocks found in source", stream=stream)
        char_matrices = []
        for block in blocks:
            # block = block.strip()
            if not block:
                continue
            src = StringIO(block)
            subproduct = self._phylip_reader._read(
                    stream=src,
                    taxon_namespace_factory=taxon_namespace_factory,
                    tree_list_factory=tree_list_factory,
                    char_matrix_factory=char_matrix_factory,
                    state_alphabet_factory=state_alphabet_factory,
                    global_annotations_target=global_annotations_target,
                    )
            char_matrices.extend(subproduct.char_matrices)
        product = self.Product(
                taxon_namespaces=None,
                tree_lists=None,
                char_matrices=char_matrices)
        return product

