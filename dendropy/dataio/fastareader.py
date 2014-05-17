#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Implementation of FASTA-format data reader.
"""

from dendropy.dataio import ioservice
from dendropy.utility.error import DataParseError

class FastaReader(ioservice.DataReader):
    "Encapsulates loading and parsing of a FASTA format file."

    def __init__(self, **kwargs):
        self.data_type_name = kwargs.pop("data_type_name", None)
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            global_annotations_target=None):
        taxon_namespace = taxon_namespace_factory(label=None)
        char_matrix = char_matrix_factory(
                data_type_name=self.data_type_name,
                label=None,
                taxon_namespace=taxon_namespace)
        symbol_state_map = char_matrix.default_state_alphabet.full_symbol_state_map
        curr_vec = None
        curr_taxon = None
        for line_index, line in enumerate(stream):
            s = line.strip()
            if not s:
                continue
            if s.startswith('>'):
                name = s[1:].strip()
                curr_taxon = taxon_namespace.require_taxon(label=name)
                if curr_taxon in char_matrix:
                    raise DataParseError(message="FASTA error: Repeated sequence name ('{}') found".format(name), line_num=line_index + 1, stream=stream)
                if curr_vec is not None and len(curr_vec) == 0:
                    raise DataParseError(message="FASTA error: Expected sequence, but found another sequence name ('{}')".format(name), line_num=line_index + 1, stream=stream)
                curr_vec = char_matrix[curr_taxon]
            elif curr_vec is None:
                raise DataParseError(message="FASTA error: Expecting a lines starting with > before sequences", line_num=line_index + 1, stream=stream)
            else:
                for col_ind, c in enumerate(s):
                    c = c.strip()
                    if not c:
                        continue
                    try:
                        state = symbol_state_map[c]
                    except KeyError:
                        raise DataParseError(message="Unrecognized sequence symbol '{}'".format(c), line_num=line_index + 1, col_num=col_ind + 1, stream=stream)
                    curr_vec.append(state)
        product = self.Product(
                taxon_namespaces=None,
                tree_lists=None,
                char_matrices=[char_matrix])
        return product
