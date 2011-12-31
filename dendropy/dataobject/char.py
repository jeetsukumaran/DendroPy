#! /usr/bin/env pthon

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
This module handles the core definition of phylogenetic character data.
"""

import copy
from cStringIO import StringIO
from dendropy.utility import error
from dendropy.utility import iosys
from dendropy.utility import containers
from dendropy.dataobject.base import IdTagged, Annotated
from dendropy.dataobject.taxon import TaxonLinked, TaxonSetLinked, TaxonSet

class ContinuousCharElement(IdTagged):
    def __init__(self, value, column_def,  **kwargs):
        IdTagged.__init__(self, **kwargs)
        self.value = value
        self.column_def = column_def

###############################################################################
## State Alphabet Infrastructure

class StateAlphabetElement(IdTagged):
    """
    A character state definition, which can either be a fundamental state or
    a mapping to a set of other character states (for polymorphic or ambiguous
    characters).
    """

    # multistate enums
    SINGLE_STATE = 0
    AMBIGUOUS_STATE = 1
    POLYMORPHIC_STATE = 2

    def __init__(self,
                 oid=None,
                 label=None,
                 symbol=None,
                 token=None,
                 multistate=SINGLE_STATE,
                 member_states=None):
        IdTagged.__init__(self, oid=oid, label=label)
        self.symbol = symbol
        self.token = token
        self.multistate = multistate
        self.member_states = member_states

    def __str__(self):
        return str(self.symbol)

    def _is_single_state(self):
        return self.multistate == StateAlphabetElement.SINGLE_STATE
    is_single_state = property(_is_single_state)

    def _get_fundamental_states(self):
        """
        Returns value of self in terms of a set of _get_fundamental states (i.e.,
        set of single states) that correspond to this state.
        """
        if self.member_states is None:
            return set([self])
        else:
            states = set()
            for state in self.member_states:
                states.update(state._get_fundamental_states())
            return states

    fundamental_states = property(_get_fundamental_states)

    def _get_fundamental_ids(self):
        "Returns set of id's of all _get_fundamental states to which this state maps."
        return set([state.oid for state in self._get_fundamental_states()])

    fundamental_ids = property(_get_fundamental_ids)

    def _get_fundamental_symbols(self):
        "Returns set of symbols of all _get_fundamental states to which this state maps."
        return set([state.symbol for state in self._get_fundamental_states()])

    fundamental_symbols = property(_get_fundamental_symbols)

    def _get_fundamental_tokens(self):
        "Returns set of tokens of all _get_fundamental states to which this state maps."
        return set([state.token for state in self._get_fundamental_states()])

    fundamental_tokens = property(_get_fundamental_tokens)

class StateAlphabet(IdTagged, list):

    "A list of states available for a particular character type/format."

    def __init__(self, *args, **kwargs):
        IdTagged.__init__(self, *args, **kwargs)
        list.__init__(self, *args)
        self.missing = None
        self.symbol_synonyms = {}
        self.case_sensitive = kwargs.get('case_sensitive', False)

    def get_state(self, attr_name, value):
        "Returns state in self in which attr_name equals value."
        for state in self:
            if getattr(state, attr_name) == value:
                return state
        if attr_name == "symbol":
            if value in self.symbol_synonyms:
                return self.symbol_synonyms[value]
            if value.islower() and not self.case_sensitive:
                return self.get_state('symbol', value.upper())
        raise Exception("State with %s of '%s' not defined" % (attr_name, str(value)))

    def state_index_for_symbol(self, symbol):
        """
        Returns index of the StateAlphabetElement object corresponding to
        the given symbol.
        """
        for idx, state in enumerate(self):
            if state.symbol == symbol:
                return idx
        if value in self.symbol_synonyms:
            return self.index(self.symbol_synonyms[value])
        raise Exception("State with symbol of '%s' not defined" % symbol)

    def state_for_symbol(self, symbol):
        "Returns a StateAlphabetElement object corresponding to given symbol."
        return self.get_state('symbol', symbol)

    def symbol_state_map(self):
        """
        Returns dictionary with symbols as keys and StateAlphabetElement objects
        as values.
        """
        map = {}
        for state in self:
            map[state.symbol] = state
        map.update(self.symbol_synonyms)
        if not self.case_sensitive:
            for state in self:
                if state.symbol.islower():
                    map[state.symbol.upper()] = state
                else:
                    map[state.symbol.lower()] = state
            for symbol, state in self.symbol_synonyms.items():
                if symbol.islower():
                    map[symbol.upper()] = state
                else:
                    map[symbol.lower()] = state
        return map

    def get_legal_symbols_as_str(self):
        m = self.symbol_state_map()
        keys = m.keys()
        for k in keys:
            if len(k) > 1:
                raise ValueError('get_legal_symbols can only be called with alphabets in which all symbols are single characters')
        return "".join(keys)

    def get_states(self, oids=None, symbols=None, tokens=None):
        """
        Returns list of states with ids/symbols/tokens equal to values
        given in a list of ids/symbols/tokens (exact matches, one-to-one
        correspondence between state and attribute value in list).
        """
        if oids is not None:
            attr_name = 'oid'
            values = oids
        elif symbols is not None:
            attr_name = 'symbol'
            values = symbols
        elif tokens is not None:
            attr_name = 'token'
            values = tokens
        else:
            raise Exception("Must specify ids, symbols or tokens")
        return [self.get_state(attr_name=attr_name, value=i) for i in values]

    def match_state(self, oids=None, symbols=None, tokens=None):
        "Returns SINGLE state that has ids/symbols/tokens as member states."
        if oids is not None:
            attr_name = 'fundamental_ids'
            values = oids
        elif symbols is not None:
            attr_name = 'fundamental_symbols'
            values = symbols
        elif tokens is not None:
            attr_name = 'fundamental_tokens'
            values = tokens
        else:
            raise Exception("Must specify ids, symbols or tokens")
        if isinstance(values, list):
            values = set(values)
        elif isinstance(values, str):
            values = set([ch for ch in values])

        for state in self:
            if getattr(state, attr_name) == values:
                return state
        return None

    def id_state_map(self):
        "Returns dictionary of element id's to state objects."
        map = {}
        for state in self:
            map[state.oid] = state
        return map

    def fundamental_states(self):
        "Returns list of fundamental states of this alphabet"
        return [s for s in self if s.multistate == StateAlphabetElement.SINGLE_STATE]

    def multi_states(self):
        "Returns list of multistate states of this alphabet"
        return [s for s in self if s.multistate != StateAlphabetElement.SINGLE_STATE]

    def ambiguous_states(self):
        "Returns list of ambiguous states of this alphabet"
        return [s for s in self if s.multistate == StateAlphabetElement.AMBIGUOUS_STATE]

    def polymorphic_states(self):
        "Returns list of ambiguous states of this alphabet"
        return [s for s in self if s.multistate == StateAlphabetElement.POLYMORPHIC_STATE]

    def get_states_as_cells(self, oids=None, symbols=None, tokens=None):
        """
        Returns (plain) list of CharacterDataCell objects with values set to
        states corresponding to symbols given by `symbols`.
        """
        return [CharacterDataCell(value=s) for \
            s in self.get_states(oids=oids, symbols=symbols, tokens=tokens)]

    def get_states_as_vector(self, oids=None, symbols=None, tokens=None, **kwargs):
        """
        Returns CharacterDataVector object, with member CharacterDataCell objects
        with values set to states corresponding to symbols given by `symbols`.
        If `taxon` is given in keyword arguments, its value will be assigned
        to the `taxon` property of the CharacterDataVector.
        """
        return CharacterDataVector(self.get_states_as_cells(oids=oids, symbols=symbols, tokens=tokens), **kwargs)
    def is_gap_state(self, el):
        """
        Returns True if the Alphabet has an element designated as the gap "state"
            and `el` is this element.
        """
        try:
            return el is self.gap
        except:
            return False

###############################################################################
## Pre-defined State Alphabets

class FixedStateAlphabet(StateAlphabet):

    def __init__(self, *args, **kwargs):
        StateAlphabet.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo):
        o = self
        memo[id(self)] = o
        return o

def _add_iupac(alphabet, states, ambig):
    for sym in states:
        sae = StateAlphabetElement(symbol=sym)
        alphabet.append(sae)
        if sym == '-':
            alphabet.gap = sae
        else:
            setattr(alphabet, sym, sae)

    for a in ambig:
        k, v = a[0], a[1]
        sae = StateAlphabetElement(symbol=k,
                                   multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                   member_states=alphabet.get_states(symbols=v))
        alphabet.append(sae)
        if k == '?':
            alphabet.missing = sae
        elif k == '*':
            alphabet.stop = sae
        else:
            setattr(alphabet, k, sae)

class DnaStateAlphabet(FixedStateAlphabet):
    _states = "ACGT-"
    _ambig = (("?",('A', 'C', 'G', 'T', '-')),
              ("N",('A', 'C', 'G', 'T')),
              ("M", ('A', 'C')),
              ("R", ('A', 'G')),
              ("W",('A', 'T')),
              ("S", ('C', 'G')),
              ("Y", ('C', 'T')),
              ("K", ('G', 'T')),
              ("V", ('A', 'C', 'G')),
              ("H",('A', 'C', 'T')),
              ("D", ('A', 'G', 'T')),
              ("B", ('C', 'G', 'T')),
             )
    unknown_state_symbol = 'N'

    def __init__(self, *args, **kwargs):
        FixedStateAlphabet.__init__(self, *args, **kwargs)
        _add_iupac(self, DnaStateAlphabet._states, DnaStateAlphabet._ambig)
        self.any_residue = self.N
        self.symbol_synonyms['X'] = self.missing

class RnaStateAlphabet(FixedStateAlphabet):
    _states = "ACGU-"
    _ambig = (("?",('A', 'C', 'G', 'U', '-')),
              ("N",('A', 'C', 'G', 'U')),
              ("M", ('A', 'C')),
              ("R", ('A', 'G')),
              ("W",('A', 'U')),
              ("S", ('C', 'G')),
              ("Y", ('C', 'U')),
              ("K", ('G', 'U')),
              ("V", ('A', 'C', 'G')),
              ("H",('A', 'C', 'U')),
              ("D", ('A', 'G', 'U')),
              ("B", ('C', 'G', 'U')),
             )
    unknown_state_symbol = 'N'

    def __init__(self, *args, **kwargs):
        FixedStateAlphabet.__init__(self, *args, **kwargs)
        _add_iupac(self, RnaStateAlphabet._states, RnaStateAlphabet._ambig)
        self.any_residue = self.N
        self.symbol_synonyms['X'] = self.missing

class NucleotideStateAlphabet(DnaStateAlphabet):

    def __init__(self, *args, **kwargs):
        DnaStateAlphabet.__init__(self, *args, **kwargs)
        self.symbol_synonyms['U'] = self.state_for_symbol('T')

class ProteinStateAlphabet(FixedStateAlphabet):
    _states = "ACDEFGHIKLMNPQRSTUVWY*-"
    _ambig = (('B', ('D', 'N')),
               ('Z', ('E', 'Q')),
               ('X', ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '*')),
               ("?", ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '*', '-')),
              )
    unknown_state_symbol = 'X'
    def __init__(self, *args, **kwargs):
        FixedStateAlphabet.__init__(self, *args, **kwargs)
        _add_iupac(self, ProteinStateAlphabet._states, ProteinStateAlphabet._ambig)
        self.any_residue = self.X

class BinaryStateAlphabet(FixedStateAlphabet):

    def __init__(self, *args, **kwargs):
        FixedStateAlphabet.__init__(self, *args, **kwargs)
        self.append(StateAlphabetElement(symbol="0"))
        self.append(StateAlphabetElement(symbol="1"))
        if kwargs.get("allow_gaps", False):
            self.append(StateAlphabetElement(symbol="-"))
            self.gap = self[-1]
            if kwargs.get("allow_missing", False):
                self.missing = StateAlphabetElement(symbol="?",
                                                   multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                                   member_states=self.get_states(symbols=['0', '1', '-']))
                self.append(self.missing)
        elif kwargs.get("allow_missing", False):
            self.missing = StateAlphabetElement(symbol="?",
                                               multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                               member_states=self.get_states(symbols=['0', '1']))
            self.append(self.missing)

class RestrictionSitesStateAlphabet(BinaryStateAlphabet):

    def __init__(self, *args, **kwargs):
        BinaryStateAlphabet.__init__(self, *args, **kwargs)

class InfiniteSitesStateAlphabet(BinaryStateAlphabet):

    def __init__(self, *args, **kwargs):
        BinaryStateAlphabet.__init__(self, *args, **kwargs)


###############################################################################
## GLOBAL STATE ALPHABETS

DNA_STATE_ALPHABET = DnaStateAlphabet()
RNA_STATE_ALPHABET = RnaStateAlphabet()
NUCLEOTIDE_STATE_ALPHABET = NucleotideStateAlphabet()
PROTEIN_STATE_ALPHABET = ProteinStateAlphabet()
RESTRICTION_SITES_STATE_ALPHABET = RestrictionSitesStateAlphabet()
INFINITE_SITES_STATE_ALPHABET = InfiniteSitesStateAlphabet()


###############################################################################
## Data Containers

class CharacterType(IdTagged):
    """
    A character format or type of a particular column: i.e., maps
    a particular set of character state definitions to a column in a character matrix.
    """

    def __init__(self, *args, **kwargs):
        IdTagged.__init__(self, *args, **kwargs)
        self._state_alphabet = None
        self.id_state_map = None
        self.state_alphabet = kwargs.get("state_alphabet", None)

    def _set_state_alphabet(self, value):
        self._state_alphabet = value
        if self._state_alphabet is not None:
            self.id_state_map = self._state_alphabet.id_state_map()
        else:
            self.id_state_map = None

    def _get_state_alphabet(self):
        return self._state_alphabet

    state_alphabet = property(_get_state_alphabet, _set_state_alphabet)

class CharacterDataCell(Annotated):
    """
    A container for that holds the value for a particular cell in a matrix.

    The attributes of CharacterDataCell are:
        'value' = an instnance of a StateAlphabetElement
        'character_type' isa CharacterType or None
    """

    def __init__(self, value=None, character_type=None):
        Annotated.__init__(self)
        self.value = value
        self.character_type = character_type

    def __str__(self):
        return str(self.value)

    def __eq__(self, other):
        if isinstance(other, CharacterDataCell):
            return self.value == other.value
        else:
            return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return NotImplemented
        return not result

class CharacterDataVector(list, TaxonLinked):
    """A list of character data values for a taxon -- a row of a Character Matrix.

    The CharacterDataVector typically contains elements that are instances of
    CharacterDataCell
    """

    def __init__(self, *args, **kwargs):
        if len(args) > 2:
            raise Exception("CharacterDataVector takes at most 1 non-keyword argument, but %d given" % len(args))
        list.__init__(self, *args)
        TaxonLinked.__init__(self, **kwargs)
        self.string_sep = ''

    def set_cell_by_index(self, column_index, cell):
        """
        Sets the cell of a cell at a particular position.
        """
        while len(self) <= column_index:
            self.append(None)
        self[column_index] = cell

    def values(self):
        return [cell.value for cell in self]

    def symbols_as_list(self):
        return [str(cell.value) for cell in self]

    def symbols_as_string(self, sep=""):
        return sep.join(self.symbols_as_list())

    def __str__(self):
        return str(self.symbols_as_string())

class CharacterDataMap(dict, Annotated):
    """
    An annotable dictionary with Taxon objects as keys and
    CharacterDataVectors objects as values.
    """

    def __init__(self):
        dict.__init__(self)
        Annotated.__init__(self)

    def _get_vector_size(self):
        """
        Returns number of characters in *first* sequence.
        """
        if len(self):
            return len(self.values()[0])
        else:
            return 0
    vector_size = property(_get_vector_size, None, None, "Returns number of characters in *first* sequence")

    def __setitem__(self, key, value):
        """
        Synchronizes taxon association.
        """
        if not isinstance(value, str):
            if not isinstance(value, CharacterDataVector):
                value = CharacterDataVector(value, taxon=key)
            else:
                value.taxon = key
        dict.__setitem__(self, key, value)

    def extend_characters(self, other_map):
        """
        Extends this matrix by adding characters from sequences of taxa
        in given matrix to sequences of taxa with correspond labels in
        this one. Taxa in the second matrix that do not exist in the
        current one are ignored.
        """
        label_taxon_map = dict([(taxon.label, taxon) for taxon in other_map])
        for taxon in self:
            if taxon.label in label_taxon_map:
                self[taxon].extend(other_map[label_taxon_map[taxon.label]])

    def extend(self,
            other_map,
            overwrite_existing=False,
            extend_existing=False):
        """
        Extends this matrix by adding taxa and characters from the given
        matrix to this one.  If `overwrite_existing` is True and a taxon
        in the other matrix is already present in the current one, then
        the sequence associated with the taxon in the second matrix
        replaces the sequence in the current one. If `extend_existing`
        is True and a taxon in the other matrix is already present in
        the current one, then the squence associated with the taxon in
        the second matrix will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True,  and a taxon in the other matrix is already present in
        the current one, then the sequence is ignored.
        Note that the containing CharacterMatrix taxa has to be normalized
        after this operation.
        """
        if overwrite_existing and extend_existing:
            raise Exception("Can only specify to overwrite or append, not both")
        label_taxon_map = dict([(taxon.label, taxon) for taxon in self])
        for other_taxon in other_map:
            if other_taxon.label in label_taxon_map:
                this_taxon = label_taxon_map[other_taxon.label]
                if overwrite_existing:
                    self[this_taxon] = other_map[other_taxon]
                elif extend_existing:
                    self[this_taxon].extend(other_map[other_taxon])
            else:
                self[other_taxon] = other_map[other_taxon]

###############################################################################
## Subset of Character (Columns)

class CharacterSubset(IdTagged):
    """
    Tracks definition of a subset of characters.
    """

    def __init__(self, *args, **kwargs):
        """
        Keyword arguments:

            - `label`: name of this subset
            - `character_indices`: list of 0-based (integer) indices
               of column positions that constitute this subset.

        """
        IdTagged.__init__(self, *args, **kwargs)
        self.character_indices = set(kwargs.get("character_indices", []))

    def __len__(self):
        return len(self.character_indices)

    def __iter__(self):
        return iter(self.character_indices)

###############################################################################
## Base Character Matrix

class CharacterMatrix(TaxonSetLinked, iosys.Readable, iosys.Writeable):
    "Character data container/manager manager."

    def concatenate(cls, char_matrices):
        """
        Creates and returns a single character matrix from multiple
        CharacterMatrix objects specified as a list, 'char_matrices'.
        All the CharacterMatrix objects in the list must be of the
        same type, and share the same TaxonSet reference. All taxa
        must be present in all alignments, all all alignments must
        be of the same length. Component parts will be recorded as
        character subsets.
        """
        taxon_set = char_matrices[0].taxon_set
        nseqs = len(char_matrices[0])
        concatenated_chars = cls(taxon_set=taxon_set)
        pos_start = 0
        for cidx, cm in enumerate(char_matrices):
            if cm.taxon_set is not taxon_set:
                raise ValueError("Different `taxon_set` references in matrices to be merged")
            if len(cm) != len(taxon_set):
                raise ValueError("Number of sequences not equal to the number of taxa")
            if len(cm) != nseqs:
                raise ValueError("Different number of sequences across alignments: %d (expecting %d based on first matrix)" % (len(cm), nseqs))
            v1 = len(cm[0])
            for t, s in cm.items():
                if len(s) != v1:
                    raise ValueError("Unequal length sequences in character matrix %d".format(cidx+1))
            concatenated_chars.extend(cm,
                    extend_existing=True,
                    overwrite_existing=False)
            if cm.label is None:
                new_label = "locus%03d" % cidx
            else:
                new_label = cm.label
            cs_label = new_label
            i = 2
            while cs_label in concatenated_chars.character_subsets:
                label = "%s_%03d" % (new_label, i)
                i += 1
            character_indices = range(pos_start, pos_start + cm.vector_size)
            pos_start += cm.vector_size
            concatenated_chars.new_character_subset(character_indices=character_indices,
                    label=cs_label)
        return concatenated_chars
    concatenate = classmethod(concatenate)

    def concatenate_from_streams(cls, streams, schema, **kwargs):
        """
        Read a character matrix from each file object given in `streams`,
        assuming data format/schema `schema`, and passing any keyword arguments
        down to the underlying specialized reader. Merge the character matrices
        and return the combined character matrix. Component parts will be
        recorded as character subsets.
        """
        if 'taxon_set' not in kwargs:
            taxon_set = TaxonSet()
            kwargs["taxon_set"] = taxon_set
        else:
            taxon_set = kwargs["taxon_set"]
        char_matrices = []
        for stream in streams:
            char_matrices.append(cls.get_from_stream(stream,
                schema=schema, **kwargs))
        return cls.concatenate(char_matrices)
    concatenate_from_streams = classmethod(concatenate_from_streams)

    def concatenate_from_paths(cls, paths, schema, **kwargs):
        """
        Read a character matrix from each file path given in `paths`, assuming
        data format/schema `schema`, and passing any keyword arguments down to
        the underlying specialized reader. Merge the and return the combined
        character matrix. Component parts will be recorded as character
        subsets.
        """
        streams = [open(path, "rU") for path in paths]
        return cls.concatenate_from_streams(streams, schema, **kwargs)
    concatenate_from_paths = classmethod(concatenate_from_paths)

    def __init__(self, *args, **kwargs):
        """__init__ calls TaxonSetLinked.__init__ for handling of `oid`, `label` and `taxon_set` keyword arguments.

        Can be initialized with:

            - source keyword arguments (see Readable.process_source_kwargs), or
            - a single unnamed CharacterMatrix instance (which will be deep-copied).

        """
        TaxonSetLinked.__init__(self,
                                taxon_set=kwargs.get("taxon_set", None),
                                label=kwargs.get("label", None),
                                oid=kwargs.get("oid", None))
        self.taxon_seq_map = CharacterDataMap()
        self.character_types = []
        self.character_subsets = containers.OrderedCaselessDict()
        self.markup_as_sequences = True
        if len(args) > 1:
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        if len(args) == 1:
            if ("stream" in kwargs and kwargs["stream"] is not None) \
                    or ("schema" in kwargs and kwargs["schema"] is not None):
                raise error.MultipleInitializationSourceError(class_name=self.__class__.__name__, arg=args[0])
            if isinstance(args[0], self.__class__):
                self.clone_from(args[0])
            else:
                raise error.InvalidArgumentValueError(func_name=self.__class__.__name__, arg=args[0])
        else:
            self.process_source_kwargs(**kwargs)
        if "oid" in kwargs:
            self.oid = kwargs["oid"]
        if "label" in kwargs:
            self.label = kwargs["label"]


    def add_character_subset(self, char_subset):
        """
        Adds a CharacterSubset object. Raises an error if one already exists
        with the same label.
        """
        label = char_subset.label
        if label in self.character_subsets:
            raise ValueError("Character subset '%s' already defined" % label)
        self.character_subsets[label] = char_subset
        return self.character_subsets[label]

    def new_character_subset(self, label, character_indices):
        """
        Defines a set of character (columns) that make up a character set.
        Raises an error if one already exists with the same label. Column
        indices are 0-based.
        """
        cs = CharacterSubset(character_indices=character_indices, label=label)
        return self.add_character_subset(cs)

    def create_taxon_to_state_set_map(self, char_indices=None):
        """Returns a dictionary that maps taxon objects to lists of sets of state
        indices
        if `char_indices` is not None it should be a iterable collection of
        character indices to include.
        """
        taxon_to_state_indices = {}
        for t in self.taxon_seq_map.keys():
            cdv = self[t]
            if char_indices is None:
                ci = range(len(cdv))
            else:
                ci = char_indices
            v = []
            for char_index in ci:
                cell = cdv[char_index]
                cell_value = cell.value
                try:
                    state_alphabet = cell.character_type.state_alphabet
                except AttributeError:
                    state_alphabet = self.default_state_alphabet
                inds = [state_alphabet.index(i) for i in cell_value.fundamental_states]
                v.append(set(inds))
            taxon_to_state_indices[t] = v
        return taxon_to_state_indices

    def clone_from(self, *args):
        "\TODO: may need to check that we are not overwriting oid"
        if len(args) > 1:
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        elif len(args) == 1:
            if isinstance(args[0],  self.__class__):
                ca = copy.deepcopy(args[0])
                self.__dict__ = ca.__dict__
            else:
                raise error.InvalidArgumentValueError(func_name=self.__class__.__name__, arg=args[0])
        return self

    def read(self, stream, schema, **kwargs):
        """
        Populates objects of this type from `schema`-formatted
        data in the file-like object source `stream`, *replacing*
        all current data. If multiple character matrices are in the data
        source, a 0-based index of the character matrix to use can
        be specified using the `matrix_offset` keyword (defaults to 0, i.e., first
        character matrix).
        """
        from dendropy.dataobject.dataset import DataSet
        index = kwargs.get("matrix_offset", 0)
        kwargs["exclude_chars"] = False
        kwargs["exclude_trees"] = True
        if 'data_type' not in kwargs and 'char_matrix_type' not in kwargs:
            kwargs['char_matrix_type'] = self.__class__
        d = DataSet(stream=stream, schema=schema, **kwargs)
        if len(d.char_matrices) == 0:
            raise ValueError("No character data in data source")
        if index >= len(d.char_matrices):
            raise IndexError("Character matrix of offset %d specified, but data source only has %d matrices defined (max. index=%d)" \
                % (index, len(d.char_matrices), len(d.char_matrices)-1))
        if not isinstance(self, d.char_matrices[index].__class__):
            raise ValueError("Character data found was of type '%s' (object is of type '%s')" %
                    (d.char_matrices[index].__class__.__name__, self.__class__.__name__))
        self.clone_from(d.char_matrices[index])

    def write(self, stream, schema, **kwargs):
        """
        Writes out this object's data to a file-like object opened for writing
        `stream`.
        """
        from dendropy.dataobject.dataset import DataSet
        d = DataSet()
        d.add(self)
        d.write(stream=stream, schema=schema, **kwargs)

    def extend_characters(self, other_matrix):
        """
        Extends this matrix by adding characters from sequences of taxa
        in given matrix to sequences of taxa with correspond labels in
        this one. Taxa in the second matrix that do not exist in the
        current one are ignored.
        """
        self.taxon_seq_map.extend_characters(other_matrix.taxon_seq_map)

    def extend_map(self,
                      other_map,
                      overwrite_existing=False,
                      extend_existing=False):
        """
        Extends this matrix by adding taxa and characters from the given
        map to this one.  If `overwrite_existing` is True and a taxon
        in the other map is already present in the current one, then
        the sequence associated with the taxon in the second map
        replaces the sequence in the current one. If `extend_existing`
        is True and a taxon in the other matrix is already present in
        the current one, then the squence map with the taxon in
        the second map will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True,  and a taxon in the other map is already present in
        the current one, then the sequence is ignored.
        """
        self.taxon_seq_map.extend(other_map,
            overwrite_existing=overwrite_existing,
            extend_existing=extend_existing)
        self.update_taxon_set()

    def extend(self,
               other_matrix,
               overwrite_existing=False,
               extend_existing=False):
        """
        Extends this matrix by adding taxa and characters from the given
        matrix to this one.  If `overwrite_existing` is True and a taxon
        in the other matrix is already present in the current one, then
        the sequence associated with the taxon in the second matrix
        replaces the sequence in the current one. If `extend_existing`
        is True and a taxon in the other matrix is already present in
        the current one, then the sequence associated with the taxon in
        the second matrix will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True, and a taxon in the other matrix is already present in
        the current one, then the sequence is ignored.
        """
        self.taxon_seq_map.extend(other_matrix.taxon_seq_map,
            overwrite_existing=overwrite_existing,
            extend_existing=extend_existing)
        self.update_taxon_set()

    def export_character_subset(self, character_subset):
        """
        Returns a new CharacterMatrix (of the same type) consisting only
        of columns given by the CharacterSubset, `character_subset`.
        Note that this new matrix will still reference the same taxon set.
        """
        if isinstance(character_subset, str):
            if character_subset not in self.character_subsets:
                raise KeyError(character_subset)
            else:
                character_subset = self.character_subsets[character_subset]
        return self.export_character_indices(character_subset.character_indices)

    def export_character_indices(self, indices):
        """
        Returns a new CharacterMatrix (of the same type) consisting only
        of columns given by the 0-based indices in `indices`.
        Note that this new matrix will still reference the same taxon set.
        """
        clone = self.__class__()
        clone.clone_from(self)
        for vec in clone.taxon_seq_map.values():
            for cell_idx in range(len(vec)-1, -1, -1):
                if cell_idx not in indices:
                    del(vec[cell_idx])
        return clone

    def reindex_subcomponent_taxa(self):
        """
        Synchronizes `Taxon` objects of map to `taxon_set` of self.
        """
        ti_mutable = self.taxon_set._is_mutable
        self.taxon_set._is_mutable = True
        new_map = CharacterDataMap()
        for taxon, seq in self.taxon_seq_map.items():
            taxon = self.taxon_set.require_taxon(label=taxon.label)
            new_map[taxon] = seq
        self.taxon_set._is_mutable = ti_mutable
        self.taxon_seq_map = new_map

    def update_taxon_set(self):
        """
        Updates local taxa block by adding taxa not already managed.
        Mainly for use after map extension
        """
        assert self.taxon_set is not None
        for taxon in self.taxon_seq_map:
            if taxon not in self.taxon_set:
                self.taxon_set.add(taxon)

    def vectors(self):
        "Returns list of vectors.        "
        if self.taxon_set is not None and self.taxon_seq_map is not None:
            if len(self.taxon_seq_map) > 0:
                return [self.taxon_seq_map[t] for t in self.taxon_set]
            return []
        return None

    # following allows a CharacterMatrix object to simulate a dictionary
    # by `passing-through` calls to the underlying character map

    def __len__(self):
        "Dictionary interface implementation for direct access to character map."
        return len(self.taxon_seq_map)

    def _get_vector_size(self):
        return self.taxon_seq_map.vector_size
    vector_size = property(_get_vector_size, None, None, "Returns number of characters in *first* sequence")

    def __getitem__(self, key):
        "Dictionary interface implementation for direct access to character map."
        if isinstance(key, int):
            if key >= 0 and key < len(self.taxon_set):
                key = self.taxon_set[key]
            else:
                raise KeyError(key)
        elif isinstance(key, str):
            label = key
            key = None
            for t in self.taxon_set:
                if t.label == label:
                    key = t
                    break
            if key is None:
                raise KeyError(label)
        return self.taxon_seq_map[key]

    def __setitem__(self, key, value):
        "Dictionary interface implementation for direct access to character map."
        if isinstance(key, int):
            if key >= 0 and key < len(self.taxon_set):
                key = self.taxon_set[key]
            else:
                raise KeyError(key)
        if key not in self.taxon_set:
            self.taxon_set.add(key)
        self.taxon_seq_map[key] = value

    def iterkeys(self):
        "Dictionary interface implementation for direct access to character map."
        for t in self.taxon_set:
            if t in self.taxon_seq_map:
                yield t

    def itervalues(self):
        "Dictionary interface implementation for direct access to character map."
        for t in self.taxon_set:
            if t in self.taxon_seq_map:
                yield self.taxon_seq_map[t]

    def iteritems(self):
        "Returns an iterator over character map's values."
        for t in self.taxon_set:
            if t in self.taxon_seq_map:
                yield t, self.taxon_seq_map[t]

    def items(self):
        "Returns character map key, value pairs in key-order."
        return [(t, self.taxon_seq_map[t]) for t in self.taxon_set if t in self.taxon_seq_map]

    def values(self):
        "Returns list of values."
        return [self.taxon_seq_map[t] for t in self.taxon_set if t in self.taxon_seq_map]

    def __iter__(self):
        "Returns an iterator over character map's ordered keys."
        for t in self.taxon_set:
            if t in self.taxon_seq_map:
                yield t

    def __delitem__(self, key):
        "Remove item from character map with specified key."
        return self.taxon_seq_map.__delitem__(key)

    def __contains__(self, key):
        "Returns true if character map has key, regardless of case."
        return self.taxon_seq_map.__contains__(key)

    def pop(self, key, alt_val=None):
        "a.pop(k[, x]):  a[k] if k in a, else x (and remove k)"
        return self.taxon_seq_map.pop(key, alt_val)

    def popitem(self):
        "a.popitem()  remove and last (key, value) pair"
        return self.taxon_seq_map.popitem()

    def keys(self):
        "Returns a copy of the ordered list of character map keys."
        return list(self.taxon_seq_map.keys())

    def clear(self):
        "Deletes all items from the character map dictionary."
        self.taxon_seq_map.clear()

    def has_key(self, key):
        "Returns true if character map has key, regardless of case."
        return key in self.taxon_seq_map

    def get(self, key, def_val=None):
        "Gets an item from character map by its key, returning default if key not present."
        return self.taxon_seq_map.get(key, def_val)

    def setdefault(self, key, def_val=None):
        "Sets the default value to return if key not present."
        return self.taxon_seq_map.setdefault(key, def_val)

    def id_chartype_map(self):
        """
        Returns dictionary of element id to corresponding
        character definition.
        """
        map = {}
        for char in self.character_types:
            map[char.oid] = char
        return map

    def description(self, depth=1, indent=0, itemize="", output=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s)" % self.oid
        else:
            label = " (%s: '%s')" % (self.oid, self.label)
        output_strio.write('%s%s%s object at %s%s'
                % (indent*' ',
                   itemize,
                   self.__class__.__name__,
                   hex(id(self)),
                   label))
        if depth >= 1:
            output_strio.write(':  %d Sequences' % len(self))
            if depth >= 2:
                if self.taxon_set is not None:
                    tlead = "\n%s[Taxon Set]\n" % (" " * (indent+4))
                    output_strio.write(tlead)
                    self.taxon_set.description(depth=depth-1, indent=indent+8, itemize="", output=output_strio)
                tlead = "\n%s[Characters]\n" % (" " * (indent+4))
                output_strio.write(tlead)
                indent += 8
                maxlabel = max([len(str(t.label)) for t in self.taxon_set])
                for i, t in enumerate(self.taxon_set):
                    output_strio.write('%s%s%s : %s characters\n' \
                        % (" " * indent,
                           "[%d] " % i,
                           str(t.label),
                           len(self.taxon_seq_map[t])))

        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

###############################################################################
## Specialized Matrices

class ContinuousCharacterMatrix(CharacterMatrix):
    "Character data container/manager manager."

    def __init__(self, *args, **kwargs):
        "See CharacterMatrix.__init__ documentation"
        CharacterMatrix.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo):
        o = TaxonSetLinked.__deepcopy__(self, memo)
        for k, v in self.__dict__.iteritems():
            if k not in ["taxon_set",
                         "_oid",
                         "taxon_seq_map"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        for taxon, cdv in self.taxon_seq_map.items():
            otaxon = memo[id(taxon)]
            ocdv = CharacterDataVector(oid=cdv.oid, label=cdv.label, taxon=otaxon)
            for cell in cdv:
                if cell.character_type is not None:
                    character_type = memo[id(cell.character_type)]
                else:
                    character_type = None
                ocdv.append(CharacterDataCell(value=cell.value, character_type=character_type))
            o.taxon_seq_map[otaxon] = ocdv
            memo[id(self.taxon_seq_map[taxon])] = o.taxon_seq_map[otaxon]
        return o

class DiscreteCharacterMatrix(CharacterMatrix):
    """Character data container/manager manager.

    That adds the attributes self.state_alphabets (a list of alphabets)
    and self.default_state_alphabet
    """

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        CharacterMatrix.__init__(self, **kwargs)
        self.state_alphabets = []
        self.default_state_alphabet = None
        self._default_symbol_state_map = None
        if len(args) > 0:
            self.clone_from(*args)

    def _get_default_symbol_state_map(self):
        if self._default_symbol_state_map is None and self.default_state_alphabet is not None:
            self._default_symbol_state_map = self.default_state_alphabet.symbol_state_map()
        return self._default_symbol_state_map

    default_symbol_state_map = property(_get_default_symbol_state_map)

    def append_taxon_sequence(self, taxon, state_symbols):
        if taxon not in self:
            self[taxon] = CharacterDataVector(taxon=taxon)
        for value in state_symbols:
            if isinstance(value, str):
                symbol = value
            else:
                symbol = str(value)
            self[taxon].append(CharacterDataCell(value=self.default_symbol_state_map[symbol]))

class StandardCharacterMatrix(DiscreteCharacterMatrix):
    "`standard` data."

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        DiscreteCharacterMatrix.__init__(self, **kwargs)
        if len(args) > 0:
            self.clone_from(*args)

    def __deepcopy__(self, memo):
        o = TaxonSetLinked.__deepcopy__(self, memo)
        for k, v in self.__dict__.iteritems():
            if k not in ["taxon_set",
                         "_oid",
                         "taxon_seq_map"]:
                o.__dict__[k] = copy.deepcopy(v, memo)

        for taxon, cdv in self.taxon_seq_map.items():
            otaxon = memo[id(taxon)]
            ocdv = CharacterDataVector(oid=cdv.oid, label=cdv.label, taxon=otaxon)
            for cell in cdv:
                if cell.character_type is not None:
                    character_type = memo[id(cell.character_type)]
                else:
                    character_type = None
                ocdv.append(CharacterDataCell(value=memo[id(cell.value)], character_type=character_type))
            o.taxon_seq_map[otaxon] = ocdv
            memo[id(self.taxon_seq_map[taxon])] = o.taxon_seq_map[otaxon]

        return o

    def extend(self,
               other_matrix,
               overwrite_existing=False,
               extend_existing=False):
        """
        Extends this matrix by adding taxa and characters from the given
        matrix to this one.  If `overwrite_existing` is True and a taxon
        in the other matrix is already present in the current one, then
        the sequence associated with the taxon in the second matrix
        replaces the sequence in the current one. If `extend_existing`
        is True and a taxon in the other matrix is already present in
        the current one, then the sequence associated with the taxon in
        the second matrix will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True, and a taxon in the other matrix is already present in
        the current one, then the sequence is ignored.
        """
        CharacterMatrix.extend(self,
                other_matrix=other_matrix,
                overwrite_existing=overwrite_existing,
                extend_existing=extend_existing)
        for s in other_matrix.state_alphabets:
            if s not in self.state_alphabets:
                self.state_alphabets.append(s)

class FixedAlphabetCharacterMatrix(DiscreteCharacterMatrix):

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        DiscreteCharacterMatrix.__init__(self, **kwargs)
        if len(args) > 0:
            self.clone_from(*args)

    def __deepcopy__(self, memo):
        o = TaxonSetLinked.__deepcopy__(self, memo)
        o.state_alphabets = self.state_alphabets
        memo[id(self.state_alphabets)] = o.state_alphabets
        o.default_state_alphabet = self.default_state_alphabet
        memo[id(self.default_state_alphabet)] = o.default_state_alphabet
        o._default_symbol_state_map = self._default_symbol_state_map
        memo[id(self._default_symbol_state_map)] = o._default_symbol_state_map
        o.character_types = copy.deepcopy(self.character_types, memo)
        for taxon, cdv in self.taxon_seq_map.items():
            otaxon = memo[id(taxon)]
            ocdv = CharacterDataVector(oid=cdv.oid, label=cdv.label, taxon=otaxon)
            for cell in cdv:
                if cell.character_type is not None:
                    character_type = memo[id(cell.character_type)]
                else:
                    character_type = None
                ocdv.append(CharacterDataCell(value=cell.value, character_type=character_type))
            o.taxon_seq_map[otaxon] = ocdv
            memo[id(self.taxon_seq_map[taxon])] = o.taxon_seq_map[otaxon]
        for k, v in self.__dict__.iteritems():
            if k not in ["taxon_set",
                         "_oid",
                         "state_alphabets",
                         "default_state_alphabet",
                         "_default_symbol_state_map",
                         "taxon_seq_map",
                         "character_types"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        return o

class DnaCharacterMatrix(FixedAlphabetCharacterMatrix):
    "DNA nucleotide data."

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        FixedAlphabetCharacterMatrix.__init__(self, **kwargs)
        self.default_state_alphabet = DNA_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
        if len(args) > 0:
            self.clone_from(*args)

class RnaCharacterMatrix(FixedAlphabetCharacterMatrix):
    "RNA nucleotide data."

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        FixedAlphabetCharacterMatrix.__init__(self, **kwargs)
        self.default_state_alphabet = RNA_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
        if len(args) > 0:
            self.clone_from(*args)

class NucleotideCharacterMatrix(FixedAlphabetCharacterMatrix):
    "Generic nucleotide data."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        FixedAlphabetCharacterMatrix.__init__(self, **kwargs)
        self.default_state_alphabet = NUCLEOTIDE_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
        if len(args) > 0:
            self.clone_from(*args)

class ProteinCharacterMatrix(FixedAlphabetCharacterMatrix):
    "Protein / amino acid data."

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`.
        """
        FixedAlphabetCharacterMatrix.__init__(self, **kwargs)
        self.default_state_alphabet = PROTEIN_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
        if len(args) > 0:
            self.clone_from(*args)

class RestrictionSitesCharacterMatrix(FixedAlphabetCharacterMatrix):
    "Restriction sites data."

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        FixedAlphabetCharacterMatrix.__init__(self, **kwargs)
        self.default_state_alphabet = RESTRICTION_SITES_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
        if len(args) > 0:
            self.clone_from(*args)

class InfiniteSitesCharacterMatrix(FixedAlphabetCharacterMatrix):
    "Infinite sites data."

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        FixedAlphabetCharacterMatrix.__init__(self, **kwargs)
        self.default_state_alphabet = INFINITE_SITES_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
        if len(args) > 0:
            self.clone_from(*args)

###############################################################################
## Wrappers and Convenience Functions

character_data_type_label_map = {
    'continuous' : ContinuousCharacterMatrix,
    'dna' : DnaCharacterMatrix,
    'rna' : RnaCharacterMatrix,
    'protein' : ProteinCharacterMatrix,
    'standard' : StandardCharacterMatrix,
    'restriction' : RestrictionSitesCharacterMatrix,
    'infinite' : InfiniteSitesCharacterMatrix,
}

def get_state_alphabet_from_symbols(symbols, gap_symbol="-", missing_symbol="?"):
    sa = StateAlphabet()
    for symbol in symbols:
        sa.append(StateAlphabetElement(symbol=symbol))
    if gap_symbol:
        sa.append(StateAlphabetElement(symbol=gap_symbol,
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=sa.get_states(symbols=symbols)))
    if missing_symbol:
        sa.append(StateAlphabetElement(symbol=missing_symbol,
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=sa.get_states(symbols=symbols)))
    return sa

###############################################################################
## Specialized Matrices

class SitePatterns(object):
    """
    Tracks distinct site patterns in a character matrix.
    Useful for efficient computations.
    """

    def __init__(self, matrix=None):
        pass
