#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
This module handles the core definition of phylogenetic character data.
"""

import copy
from dendropy.dataobject.base import IdTagged, Annotated
from dendropy.dataobject.taxon import TaxonLinked, TaxonSetLinked

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

    def __repr__(self):
        return str([self.oid,
                    self.symbol,
                    '[' + (', '.join(self._get_fundamental_symbols())) + ']'])

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
        return set([state.symbol for state in self._get_fundamental_states() if state.symbol])

    fundamental_symbols = property(_get_fundamental_symbols)

    def _get_fundamental_tokens(self):
        "Returns set of tokens of all _get_fundamental states to which this state maps."
        return set([state.token for state in self._get_fundamental_states() if state.token])

    fundamental_tokens = property(_get_fundamental_tokens)

class StateAlphabet(IdTagged, list):
    "A set of states available for a particular character type/format."

    def __init__(self, *args, **kwargs):
        IdTagged.__init__(self, *args, **kwargs)
        list.__init__(self, *args)
        self.missing = None

    def get_state(self, attr_name, value):
        "Returns state in self in which attr_name equals value."
        for state in self:
            if getattr(state, attr_name) == value:
                return state
        raise Exception("State with %s of '%s' not defined" % (attr_name, str(value)))

    def state_index_for_symbol(self, symbol):
        """
        Returns index of the StateAlphabetElement object corresponding to
        the given symbol.
        """
        for idx, state in enumerate(self):
            if state.symbol == symbol:
                return idx
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

class FixedStateAlphabet(StateAlphabet):

    def __init__(self, *args, **kwargs):
        StateAlphabet.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo):
        o = self
        memo[id(self)] = o
        return o

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
        for sym in DnaStateAlphabet._states:
            self.append(StateAlphabetElement(symbol=sym))
        self.gap = self[-1]
        for a in DnaStateAlphabet._ambig:
            k, v = a[0], a[1]
            sae = StateAlphabetElement(symbol=k,
                                       multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                       member_states=self.get_states(symbols=v))
            self.append(sae)
            if k == '?':
                self.missing = sae
            elif k == 'N':
                self.any_residue = sae

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
        for sym in RnaStateAlphabet._states:
            self.append(StateAlphabetElement(symbol=sym))
        self.gap = self[-1]
        for a in RnaStateAlphabet._ambig:
            k, v = a[0], a[1]
            sae = StateAlphabetElement(symbol=k,
                                       multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                       member_states=self.get_states(symbols=v))
            self.append(sae)
            if k == '?':
                self.missing = sae
            elif k == 'N':
                self.any_residue = sae

class ProteinStateAlphabet(FixedStateAlphabet):
    _states = "ACDEFGHIKLMNPQRSTUVWY-"
    _ambig = (('B', ('D', 'N')),
               ('Z', ('E', 'Q')),
               ('X', ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y')),
               ("?", ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '-')),
              )
    unknown_state_symbol = 'X'
    def __init__(self, *args, **kwargs):
        FixedStateAlphabet.__init__(self, *args, **kwargs)
        for sym in ProteinStateAlphabet._states:
            self.append(StateAlphabetElement(symbol=sym))
        self.gap = self[-1]
        for a in ProteinStateAlphabet._ambig:
            k, v = a[0], a[1]
            sae = StateAlphabetElement(symbol=k,
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=v))
            self.append(sae)
            if k == '?':
                self.missing = sae
            elif k == 'X':
                self.any_residue = sae

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

### GLOBAL STATE ALPHABETS ###

DNA_STATE_ALPHABET = DnaStateAlphabet()
RNA_STATE_ALPHABET = RnaStateAlphabet()
PROTEIN_STATE_ALPHABET = ProteinStateAlphabet()
RESTRICTION_SITES_STATE_ALPHABET = RestrictionSitesStateAlphabet()
INFINITE_SITES_STATE_ALPHABET = InfiniteSitesStateAlphabet()

class ColumnType(IdTagged):
    """
    A character format or type of a particular column: i.e., maps
    a particular set of character state definitions to a column in a character array.
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
    A container for the state / state value for a particular cell in a array.
    """

    def __init__(self, value=None, column_type=None):
        Annotated.__init__(self)
        self.value = value
        self.column_type = column_type

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
    "A list of character data values."

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

    def values_as_string_list(self):
        return [str(cell.value) for cell in self]

    def values_as_string(self, sep=""):
        return sep.join(self.values_as_string_list())

    def __str__(self):
        return str(self.values_as_string_list())

class CharacterDataMap(dict, Annotated):
    """
    An annotable dictionary with Taxon objects as keys and
    CharacterDataVectors objects as values.
    """

    def __init__(self):
        dict.__init__(self)
        Annotated.__init__(self)

    def extend_characters(self, other_map):
        """
        Extends this array by adding characters from sequences of taxa
        in given array to sequences of taxa with correspond labels in
        this one. Taxa in the second array that do not exist in the
        current one are ignored.
        """
        label_taxon_map = dict([(taxon.label, taxon) for taxon in other_map])
        for taxon in self:
            if taxon.label in label_taxon_map:
                self[taxon].extend(other_map[label_taxon_map[taxon.label]])

    def extend(self,
        other_map,
        overwrite_existing=False,
        append_existing=False):
        """
        Extends this array by adding taxa and characters from the given
        array to this one.  If `overwrite_existing` is True and a taxon
        in the other array is already present in the current one, then
        the sequence associated with the taxon in the second array
        replaces the sequence in the current one. If `append_existing`
        is True and a taxon in the other array is already present in
        the current one, then the squence associated with the taxon in
        the second array will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True,  and a taxon in the other array is already present in
        the current one, then the sequence is ignored.
        Note that the containing CharacterArray taxa has to be normalized
        after this operation.
        """
        if overwrite_existing and append_existing:
            raise Exception("Can only specify to overwrite or append, not both")
        label_taxon_map = dict([(taxon.label, taxon) for taxon in self])
        for other_taxon in other_map:
            if other_taxon.label in label_taxon_map:
                this_taxon = label_taxon_map[other_taxon.label]
                if overwrite_existing:
                    self[this_taxon] = other_map[other_taxon]
                elif append_existing:
                    self[this_taxon].extend(other_map[other_taxon])
            else:
                self[other_taxon] = other_map[other_taxon]

class CharacterArray(TaxonSetLinked):
    "Character data container/manager manager."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        TaxonSetLinked.__init__(self,
                                taxon_set=kwargs.get("taxon_set", None),
                                label=kwargs.get("label", None),
                                oid=kwargs.get("oid", None))
        self.taxon_seq_map = CharacterDataMap()
        self.column_types = []
        self.markup_as_sequences = True
        if len(args) > 1:
            raise TypeError("CharacterArray() takes a maximum of 1 non-keyword argument, but %d given: %s"\
                % (len(args), str(args)))
        elif len(args) > 0 and isinstance(args[0], CharacterArray):
            ca = copy.deepcopy(args[0])
            self.__dict__ = ca.__dict__
        elif len(args) > 0:
            raise TypeError("Invalid non-keyword argument passed to CharacterArray(): %s" % arg[0])

    def extend_characters(self, other_array):
        """
        Extends this array by adding characters from sequences of taxa
        in given array to sequences of taxa with correspond labels in
        this one. Taxa in the second array that do not exist in the
        current one are ignored.
        """
        self.taxon_seq_map.extend_characters(other_array.taxon_seq_map)

    def extend_map(self,
                      other_map,
                      overwrite_existing=False,
                      append_existing=False):
        """
        Extends this array by adding taxa and characters from the given
        map to this one.  If `overwrite_existing` is True and a taxon
        in the other map is already present in the current one, then
        the sequence associated with the taxon in the second map
        replaces the sequence in the current one. If `append_existing`
        is True and a taxon in the other array is already present in
        the current one, then the squence map with the taxon in
        the second map will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True,  and a taxon in the other map is already present in
        the current one, then the sequence is ignored.
        """
        self.taxon_seq_map.extend(other_map,
            overwrite_existing=overwrite_existing,
            append_existing=append_existing)
        self.update_taxon_set()

    def extend(self,
               other_array,
               overwrite_existing=False,
               append_existing=False):
        """
        Extends this array by adding taxa and characters from the given
        array to this one.  If `overwrite_existing` is True and a taxon
        in the other array is already present in the current one, then
        the sequence associated with the taxon in the second array
        replaces the sequence in the current one. If `append_existing`
        is True and a taxon in the other array is already present in
        the current one, then the sequence associated with the taxon in
        the second array will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True, and a taxon in the other array is already present in
        the current one, then the sequence is ignored.
        """
        self.taxon_seq_map.extend(other_array.taxon_seq_map,
            overwrite_existing=overwrite_existing,
            append_existing=append_existing)
        self.update_taxon_set()

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
        for taxon in self:
            if taxon not in self.taxon_set:
                self.taxon_set.add(taxon)

    def vectors(self):
        "Returns list of vectors.        "
        if self.taxon_set is not None and self.taxon_seq_map is not None:
            if len(self.taxon_seq_map) > 0:
                return [self.taxon_seq_map[t] for t in self.taxon_set]
            return []
        return None

    # following allows a CharacterArray object to simulate a dictionary
    # by `passing-through` calls to the underlying character map


    def __len__(self):
        "Dictionary interface implementation for direct access to character map."
        return len(self.taxon_seq_map)

    def __getitem__(self, key):
        "Dictionary interface implementation for direct access to character map."
        if isinstance(key, int):
            if key >= 0 and key < len(self.taxon_set):
                key = self.taxon_set[key]
            else:
                raise KeyError(key)
        return self.taxon_seq_map[key]

    def __setitem__(self, key, value):
        "Dictionary interface implementation for direct access to character map."
        self.taxon_seq_map[key] = value

#     def __contains__(self, key):
#         """
#         Dictionary interface implementation for direct access to character map.
#         """
#         return key in self.taxon_seq_map
#
    def iterkeys(self):
        "Dictionary interface implementation for direct access to character map."
        for key in self.taxon_seq_map:
            yield(key)

    def itervalues(self):
        "Dictionary interface implementation for direct access to character map."
        for value in self.taxon_seq_map.values():
            yield(value)

    def iteritems(self):
        "Returns an iterator over character map's values."
        for key, value in self.taxon_seq_map.iteritems():
            yield (key, value)

    def items(self):
        "Returns character map key, value pairs in key-order."
        return [(key, self.taxon_seq_map[key]) for key in self.taxon_seq_map.iterkeys()]

    def values(self):
        "Returns list of character map key, value pairs."
        return [v for v in self.taxon_seq_map.itervalues()]

    def __iter__(self):
        "Returns an iterator over character map's ordered keys."
        return self.taxon_seq_map.iterkeys()

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

    def id_column_map(self):
        """
        Returns dictionary of element id to corresponding
        character definition.
        """
        map = {}
        for char in self.column_types:
            map[char.oid] = char
        return map

class ContinuousCharacterArray(CharacterArray):
    "Character data container/manager manager."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        CharacterArray.__init__(self, *args, **kwargs)

class DiscreteCharacterArray(CharacterArray):
    """Character data container/manager manager.

    That adds the attributes self.state_alphabets (a list of alphabets)
    and self.default_state_alphabet
    """

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        CharacterArray.__init__(self, *args, **kwargs)
        self.state_alphabets = []
        self.default_state_alphabet = None
        self._default_symbol_state_map = None

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

class StandardCharacterArray(DiscreteCharacterArray):
    "`standard` data."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        DiscreteCharacterArray.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo):
        o = self.__class__(taxon_set=self.taxon_set, oid=self.oid)
        memo[id(self)] = o
        memo[id(self.taxon_set)] = o.taxon_set
        for i, t in enumerate(self.taxon_set):
            memo[id(t)] = o.taxon_set[i]
        for k, v in self.__dict__.iteritems():
            if k not in ["taxon_set",
                         "_oid",
                         "state_alphabets",
                         "default_state_alphabet",
                         "_default_symbol_state_map",
                         "taxon_seq_map",
                         "column_types"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        return o

class FixedAlphabetCharacterArray(DiscreteCharacterArray):

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        DiscreteCharacterArray.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo):
        o = self.__class__(taxon_set=self.taxon_set, oid=self.oid)
        memo[id(self)] = o
        memo[id(self.taxon_set)] = o.taxon_set
        for i, t in enumerate(self.taxon_set):
            memo[id(t)] = o.taxon_set[i]
        o.state_alphabets = self.state_alphabets
        memo[id(self.state_alphabets)] = o.state_alphabets
        o.default_state_alphabet = self.default_state_alphabet
        memo[id(self.default_state_alphabet)] = o.default_state_alphabet
        o._default_symbol_state_map = self._default_symbol_state_map
        memo[id(self._default_symbol_state_map)] = o._default_symbol_state_map
        o.column_types = copy.deepcopy(self.column_types, memo)
        for taxon, cdv in self.taxon_seq_map.items():
            ocdv = CharacterDataVector(oid=cdv.oid, label=cdv.label, taxon=taxon)
            for cell in cdv:
                if cell.column_type is not None:
                    column_type = memo[id(cell.column_type)]
                else:
                    column_type = None
                ocdv.append(CharacterDataCell(value=cell.value, column_type=column_type))
            o.taxon_seq_map[taxon] = ocdv
        memo[id(self.taxon_seq_map[taxon])] = o.taxon_seq_map[taxon]
        for k, v in self.__dict__.iteritems():
            if k not in ["taxon_set",
                         "_oid",
                         "state_alphabets",
                         "default_state_alphabet",
                         "_default_symbol_state_map",
                         "taxon_seq_map",
                         "column_types"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        return o

class DnaCharacterArray(FixedAlphabetCharacterArray):
    "DNA nucleotide data."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        FixedAlphabetCharacterArray.__init__(self, *args, **kwargs)
        self.default_state_alphabet = DNA_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)

class RnaCharacterArray(FixedAlphabetCharacterArray):
    "RNA nucleotide data."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        FixedAlphabetCharacterArray.__init__(self, *args, **kwargs)
        self.default_state_alphabet = RNA_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)

class ProteinCharacterArray(FixedAlphabetCharacterArray):
    "Protein / amino acid data."

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`.
        """
        FixedAlphabetCharacterArray.__init__(self, *args, **kwargs)
        self.default_state_alphabet = PROTEIN_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)

class RestrictionSitesCharacterArray(FixedAlphabetCharacterArray):
    "Restriction sites data."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        FixedAlphabetCharacterArray.__init__(self, *args, **kwargs)
        self.default_state_alphabet = RESTRICTION_SITES_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)

class InfiniteSitesCharacterArray(FixedAlphabetCharacterArray):
    "Infinite sites data."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_set`."
        FixedAlphabetCharacterArray.__init__(self, *args, **kwargs)
        self.default_state_alphabet = INFINITE_SITES_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
