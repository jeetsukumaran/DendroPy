#! /usr/bin/env python

"""
Implementation of FASTA-format data reader.
"""

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.dataio import ioservice
from dendropy.utility.error import DataParseError

class FastaReader(ioservice.DataReader):
    "Encapsulates loading and parsing of a FASTA format file."

    supported_data_types = ['dna', 'rna', 'protein', 'standard', 'restriction', 'infinite']

    def __init__(self, **kwargs):
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            global_annotations_target=None):
        taxon_namespace = taxon_namespace_factory(label=None)
        char_matrix = char_matrix_factory(label=None, taxon_namespace=taxon_namespace)
        symbol_state_map = char_matrix.default_symbol_state_map

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
                curr_vec = char_matrix.new_character_data_vector(taxon=curr_taxon)
                char_matrix[curr_taxon] = curr_vec
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
