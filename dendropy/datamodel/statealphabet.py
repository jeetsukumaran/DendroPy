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
State alphabets for discrete character sequence data types.
"""

from dendropy.datamodel import base

###############################################################################
## State Alphabet Infrastructure

class StateAlphabetElement(
        base.DataObject,
        base.Annotable):
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
                 label=None,
                 symbol=None,
                 token=None,
                 multistate=SINGLE_STATE,
                 member_states=None):
        base.DataObject.__init__(self, label=label)
        self.symbol = symbol
        self.token = token
        self.multistate = multistate
        self._member_states = None
        self._fundamental_states = None
        self._fundamental_symbols = None
        self._fundamental_tokens = None
        self.member_states = member_states

    def __str__(self):
        # note that tests currently assume this particular string
        # representation (i.e., undecorated symbol)
        return str(self.symbol)

    def _is_single_state(self):
        return self.multistate == StateAlphabetElement.SINGLE_STATE
    is_single_state = property(_is_single_state)

    def _get_member_states(self):
        return self._member_states
    def _set_member_states(self, m):
        self._member_states = m
        self._fundamental_states = None
        self._fundamental_symbols = None
        self._fundamental_tokens = None
    member_states = property(_get_member_states, _set_member_states)

    def _get_fundamental_states(self):
        """
        Returns value of self in terms of a set of _get_fundamental states (i.e.,
        set of single states) that correspond to this state.
        """
        if self._fundamental_states is None:
            if self._member_states is None:
                return set([self])
            else:
                states = set()
                for state in self._member_states:
                    assert state is not self
                    states.update(state.fundamental_states)
                return states
        return self._fundamental_states
    fundamental_states = property(_get_fundamental_states)

    def _get_fundamental_symbols(self):
        "Returns set of symbols of all _get_fundamental states to which this state maps."
        if self._fundamental_symbols is None:
            self._fundamental_symbols = set([state.symbol for state in self.fundamental_states])
        return self._fundamental_symbols
    fundamental_symbols = property(_get_fundamental_symbols)

    def _get_fundamental_tokens(self):
        "Returns set of tokens of all _get_fundamental states to which this state maps."
        if self._fundamental_tokens is None:
            self._fundamental_tokens = set([state.token for state in self.fundamental_states])
    fundamental_tokens = property(_get_fundamental_tokens)

    def is_exact_correspondence(self, other):
        """
        Tries to determine if two StateAlphabetElement definitions
        are equivalent by matching symbols.
        """
        match = True
        if self.multistate != other.multistate:
            return False
        if self.multistate and other.multistate:
            xf1 = self.fundamental_states
            xf2 = other.fundamental_states
            if len(xf1) != len(xf2):
                match = False
            else:
                f1 = set(xf1)
                f2 = set(xf2)
                for m1 in f1:
                    member_match = False
                    for m2 in f2:
                        if m1.is_exact_correspondence(m2):
                            member_match = True
                            f2.remove(m2)
                            break
                    if not member_match:
                        match = False
                        break
                if match:
                    f1 = set(xf1)
                    f2 = set(xf2)
                    for m2 in f2:
                        member_match = False
                        for m1 in f1:
                            if m1.is_exact_correspondence(m2):
                                f1.remove(m1)
                                member_match = True
                                break
                        if not member_match:
                            match = False
                            break
            return match
        else:
            try:
                return self.symbol.upper() == other.symbol.upper()
            except AttributeError:
                return self.symbol == other.symbol

###############################################################################
## FixedStateAlphabetElement

class FixedStateAlphabetElement(StateAlphabetElement):
    """
    Specialized for fixed state alphabets (e.g. DNA).
    """

    def __init__(self,
                 label=None,
                 symbol=None,
                 token=None,
                 multistate=StateAlphabetElement.SINGLE_STATE,
                 member_states=None):
         StateAlphabetElement.__init__(self,
                 label=label,
                 symbol=symbol,
                 token=token,
                 multistate=multistate,
                 member_states=member_states)

    def __deepcopy__(self, memo):
        memo[id(self)] = self
        return self

###############################################################################
## StateAlphabet

class StateAlphabet(
        base.DataObject,
        base.Annotable,
        ):

    "A list of states available for a particular character type/format."

    ###########################################################################
    ### Life-Cycle and Identity

    def __init__(self, *args, **kwargs):
        base.DataObject.__init__(self,
                label=kwargs.pop("label", None))
        self._alphabet_elements = list(*args)
        self.missing = None
        self.symbol_synonyms = {}
        self.case_sensitive = kwargs.pop('case_sensitive', False)
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def __copy__(self, memo=None):
        raise TypeError("Cannot directly copy {}".format(self.__class__.__name__))

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot directly copy {}".format(self.__class__.__name__))

    def __deepcopy__(self, memo=None):
        return base.Annotable.__deepcopy__(self, memo=memo)

    ###########################################################################
    ### List Access and Interface

    def _get_alphabet_elements(self):
        return self._alphabet_elements
    def _set_alphabet_elements(self, elements):
        self._alphabet_elements = elements
    elements = property(_get_alphabet_elements, _set_alphabet_elements)

    def insert(self, index, element):
        return self._alphabet_elements.insert(index, element)

    def append(self, element):
        return self._alphabet_elements.append(element)

    def extend(self, other):
        self._alphabet_elements.extend(other)

    def __iadd__(self, other):
        raise NotImplementedError

    def __add__(self, other):
        raise NotImplementedError

    def __contains__(self, element):
        raise NotImplementedError

    def __delitem__(self, element):
        raise NotImplementedError

    def __iter__(self):
        return iter(self._alphabet_elements)

    def __reversed__(self):
        return reversed(self._alphabet_elements)

    def __len__(self):
        return len(self._alphabet_elements)

    def __getitem__(self, index):
        if isinstance(index, slice):
            raise NotImplementedError
        else:
            return self._alphabet_elements[index]

    def __setitem__(self, index, value):
        if isinstance(index, slice):
            raise NotImplementedError
        else:
            self._alphabet_elements[index] = value

    def clear(self):
        # list.clear() only with 3.4 or so ...
        self._alphabet_elements = []

    def index(self, element):
        return self._alphabet_elements.index(element)

    def pop(self, index=-1):
        return self._alphabet_elements.pop(index)

    def remove(self, element):
        self._alphabet_elements.remove(element)

    def reverse(self):
        self._alphabet_elements.reverse()

    def sort(self, key=None, reverse=False):
        self._alphabet_elements.sort(key=key, reverse=reverse)

    ###########################################################################
    ### Specialized Symbol Management

    def get_state(self, attr_name, value):
        "Returns state in self in which attr_name equals value."
        for state in self._alphabet_elements:
            if getattr(state, attr_name) == value:
                return state
        if attr_name == "symbol":
            if value in self.symbol_synonyms:
                return self.symbol_synonyms[value]
            if value.islower() and not self.case_sensitive:
                return self.get_state('symbol', value.upper())
        raise Exception("State with {} of '{}' not defined".format(attr_name, str(value)))

    def state_index_for_symbol(self, symbol):
        """
        Returns index of the StateAlphabetElement object corresponding to
        the given symbol.
        """
        for idx, state in enumerate(self._alphabet_elements):
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
        for state in self._alphabet_elements:
            map[state.symbol] = state
        map.update(self.symbol_synonyms)
        if not self.case_sensitive:
            for state in self._alphabet_elements:
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

    def get_states(self, symbols=None, tokens=None):
        """
        Returns list of states with ids/symbols/tokens equal to values
        given in a list of ids/symbols/tokens (exact matches, one-to-one
        correspondence between state and attribute value in list).
        """
        if symbols is not None:
            attr_name = 'symbol'
            values = symbols
        elif tokens is not None:
            attr_name = 'token'
            values = tokens
        else:
            raise Exception("Must specify ids, symbols or tokens")
        return [self.get_state(attr_name=attr_name, value=i) for i in values]

    def match_state(self, symbols=None, tokens=None):
        "Returns SINGLE state that has ids/symbols/tokens as member states."
        if symbols is not None:
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

        for state in self._alphabet_elements:
            if getattr(state, attr_name) == values:
                return state
        return None

    def fundamental_states(self):
        "Returns list of fundamental states of this alphabet"
        return [s for s in self._alphabet_elements if s.multistate == StateAlphabetElement.SINGLE_STATE]

    def multi_states(self):
        "Returns list of multistate states of this alphabet"
        return [s for s in self._alphabet_elements if s.multistate != StateAlphabetElement.SINGLE_STATE]

    def ambiguous_states(self):
        "Returns list of ambiguous states of this alphabet"
        return [s for s in self._alphabet_elements if s.multistate == StateAlphabetElement.AMBIGUOUS_STATE]

    def polymorphic_states(self):
        "Returns list of ambiguous states of this alphabet"
        return [s for s in self._alphabet_elements if s.multistate == StateAlphabetElement.POLYMORPHIC_STATE]

    def get_states_as_cells(self, symbols=None, tokens=None):
        """
        Returns (plain) list of CharacterDataCell objects with values set to
        states corresponding to symbols given by `symbols`.
        """
        return [CharacterDataCell(value=s) for \
            s in self.get_states(symbols=symbols, tokens=tokens)]

    def get_states_as_vector(self, symbols=None, tokens=None, **kwargs):
        """
        Returns CharacterDataVector object, with member CharacterDataCell objects
        with values set to states corresponding to symbols given by `symbols`.
        If `taxon` is given in keyword arguments, its value will be assigned
        to the `taxon` property of the CharacterDataVector.
        """
        return CharacterDataVector(self.get_states_as_cells(symbols=symbols, tokens=tokens), **kwargs)

    def is_gap_state(self, el):
        """
        Returns True if the Alphabet has an element designated as the gap "state"
            and `el` is this element.
        """
        try:
            return el is self.gap
        except:
            return False

    def is_exact_correspondence(self,
            other,
            accept_other_as_subset=False,
            symbols_to_ignore=None):
        n1 = len(self._alphabet_elements)
        n2 = len(other._alphabet_elements)
        if n1 > n2:
            if accept_other_as_subset:
                sa1 = self
                sa2 = other
            else:
                return False
        elif n1 < n2:
            return False
        else:
            sa1 = self
            sa2 = other
        match = True
        f1 = set(sa1._alphabet_elements)
        f2 = set(sa2._alphabet_elements)
        if symbols_to_ignore is None:
            symbols_to_ignore = []
        for m1 in f1:
            if m1.symbol in symbols_to_ignore:
                continue
            member_match = False
            for m2 in f2:
                if m1.is_exact_correspondence(m2):
                    member_match = True
                    f2.remove(m2)
                    break
            if not member_match:
                match = False
                break
        f1 = set(sa1._alphabet_elements)
        f2 = set(sa2._alphabet_elements)
        if match:
            for m2 in f2:
                if m2.symbol in symbols_to_ignore:
                    continue
                member_match = False
                for m1 in f1:
                    if m1.is_exact_correspondence(m2):
                        member_match = True
                        break
                if not member_match:
                    match = False
                    break

        return match

###############################################################################
## Pre-defined State Alphabets

class FixedStateAlphabet(StateAlphabet):

    def __init__(self, *args, **kwargs):
        StateAlphabet.__init__(self, *args, **kwargs)

    def __deepcopy__(self, memo):
        memo[id(self)] = self
        return self

def _add_iupac(alphabet, states, ambig):
    for sym in states:
        sae = FixedStateAlphabetElement(symbol=sym)
        alphabet.append(sae)
        if sym == '-':
            alphabet.gap = sae
        else:
            setattr(alphabet, sym, sae)

    for a in ambig:
        k, v = a[0], a[1]
        sae = FixedStateAlphabetElement(symbol=k,
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
              # Added 'X', 2012-07-17: this is how Rutger's BioPhylo treats 'X",
              # and this definition is to allow for state alphabet
              # equivalency with files produced by that library
              # have not checked how this interacts with 'symbol_synonyms' below
              ("X",('A', 'C', 'G', 'T')),
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
              # Added 'X', 2012-07-17: this is how Rutger's BioPhylo treats 'X",
              # and this definition is to allow for state alphabet
              # equivalency with files produced by that library
              # have not checked how this interacts with 'symbol_synonyms' below
              ("X",('A', 'C', 'G', 'U')),
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
        self.append(FixedStateAlphabetElement(symbol="0"))
        self.append(FixedStateAlphabetElement(symbol="1"))
        if kwargs.get("allow_gaps", False):
            self.append(FixedStateAlphabetElement(symbol="-"))
            self.gap = self[-1]
            if kwargs.get("allow_missing", False):
                self.missing = FixedStateAlphabetElement(symbol="?",
                                                   multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                                   member_states=self.get_states(symbols=['0', '1', '-']))
                self.append(self.missing)
        elif kwargs.get("allow_missing", False):
            self.missing = FixedStateAlphabetElement(symbol="?",
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

