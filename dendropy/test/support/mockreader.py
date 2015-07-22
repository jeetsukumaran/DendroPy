#! /usr/bin/env python

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

from dendropy.dataio import ioservice
from dendropy import dataio

class MockReader(ioservice.DataReader):

    def __init__(self):
        self.stream = None
        self.taxon_namespace_factory = None
        self.tree_list_factory = None
        self.char_matrix_factory = None
        self.state_alphabet_factory = None
        self.global_annotations_target = None

    def _read(self,
            stream,
            taxon_namespace_factory,
            tree_list_factory,
            char_matrix_factory,
            state_alphabet_factory,
            global_annotations_target):
        self.stream = stream
        self.taxon_namespace_factory = taxon_namespace_factory
        self.tree_list_factory = tree_list_factory
        self.char_matrix_factory = char_matrix_factory
        self.state_alphabet_factory = state_alphabet_factory
        self.global_annotations_target = global_annotations_target
        return self.process_read_call()

    def process_read_call(self):
        return None

