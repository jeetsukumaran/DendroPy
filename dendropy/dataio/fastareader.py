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

"""
Implementation of FASTA-format data reader.
"""

from dendropy.dataio import ioservice
from dendropy.utility.error import DataParseError
from dendropy.utility import deprecate

class FastaReader(ioservice.DataReader):
    "Encapsulates loading and parsing of a FASTA format file."

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------
        data_type: str
            When reading into a |DataSet| object, the type of data must be
            specified: "dna", "rna", "protein", "restriction", "infinite",
            "standard", or "continuous".
        default_state_alphabet: |StateAlphabet| instance
            A |StateAlphabet| object to be used to manage the alphabet of the
            characters (|StandardCharacterMatrix| **only**).
        """
        ioservice.DataReader.__init__(self)
        self.data_type = kwargs.pop("data_type", None)
        self.default_state_alphabet = kwargs.pop("default_state_alphabet", None)
        if self.default_state_alphabet is not None:
            if self.data_type is None:
                self.data_type = "standard"
            elif self.data_type != "standard":
                raise ValueError("Cannot specify 'default_state_alphabet' with data type of '{}'".format(self.data_type))
        self.check_for_unused_keyword_arguments(kwargs)

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            state_alphabet_factory=None,
            global_annotations_target=None):
        taxon_namespace = taxon_namespace_factory(label=None)
        if self.data_type is None:
            raise TypeError("Data type must be specified for this schema")
        if self.data_type == "standard" and self.default_state_alphabet is not None:
            char_matrix = char_matrix_factory(
                    self.data_type,
                    label=None,
                    taxon_namespace=taxon_namespace,
                    default_state_alphabet=self.default_state_alphabet,
                    )
        else:
            char_matrix = char_matrix_factory(
                    self.data_type,
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
                states = []
                for col_ind, c in enumerate(s):
                    c = c.strip()
                    if not c:
                        continue
                    try:
                        state = symbol_state_map[c]
                    except KeyError:
                        raise DataParseError(message="Unrecognized sequence symbol '{}'".format(c), line_num=line_index + 1, col_num=col_ind + 1, stream=stream)
                    states.append(state)
                curr_vec.extend(states)
        product = self.Product(
                taxon_namespaces=None,
                tree_lists=None,
                char_matrices=[char_matrix])
        return product


class DnaFastaReader(FastaReader):

    def __init__(self, **kwargs):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4:",
                old_construct="d = dendropy.CharacterMatrix.get_from_path(schema='dnafasta', ...)\nd = dendropy.DataSet.get_from_path(schema='dnafasta', ...)",
                new_construct="d = dendropy.DnaCharacterMatrix.get(path=..., schema='fasta', ...)\nd = dendropy.DataSet.get(path=..., schema='fasta', data_type='dna', ...)",
                stacklevel=7)
        # raise TypeError("'dnafasta' is no longer a supported schema: use 'schema=\"fasta\"' with the 'DnaCharacterMatrix.get()' method instead or 'schema=\"fasta\"' and 'data_type=\"dna\" with the 'DataSet.get()' or 'DataSet.read()' methods")
        kwargs["data_type"] = "dna"
        FastaReader.__init__(self, **kwargs)

class RnaFastaReader(FastaReader):

    def __init__(self, **kwargs):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4:",
                old_construct="d = dendropy.CharacterMatrix.get_from_path(schema='rnafasta', ...)\nd = dendropy.DataSet.get_from_path(schema='rnafasta', ...)",
                new_construct="d = dendropy.RnaCharacterMatrix.get(path=..., schema='fasta', ...)\nd = dendropy.DataSet.get(path=..., schema='fasta', data_type='rna', ...)",
                stacklevel=7)
        # raise TypeError("'rnafasta' is no longer a supported schema: use 'schema=\"fasta\"' with the 'RnaCharacterMatrix.get()' method instead or 'schema=\"fasta\"' and 'data_type=\"dna\" with the 'DataSet.get()' or 'DataSet.read()' methods")
        kwargs["data_type"] = "rna"
        FastaReader.__init__(self, **kwargs)

class ProteinFastaReader(FastaReader):

    def __init__(self, **kwargs):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4:",
                old_construct="d = dendropy.CharacterMatrix.get_from_path(schema='proteinfasta', ...)\nd = dendropy.DataSet.get_from_path(schema='proteinfasta', ...)",
                new_construct="d = dendropy.ProteinCharacterMatrix.get(path=..., schema='fasta', ...)\nd = dendropy.DataSet.get(path=..., schema='fasta', data_type='protein', ...)",
                stacklevel=7)
        # raise TypeError("'proteinfasta' is no longer a supported schema: use 'schema=\"fasta\"' with the 'ProteinCharacterMatrix.get()' method instead or 'schema=\"fasta\"' and 'data_type=\"dna\" with the 'DataSet.get()' or 'DataSet.read()' methods")
        kwargs["data_type"] = "protein"
        FastaReader.__init__(self, **kwargs)

