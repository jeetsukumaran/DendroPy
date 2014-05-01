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
Character state definitions and alphabets. Certain state alphabets, such as
DNA, RNA, protein, etc.  are defined here. These are termed "fixed" state
alphabets, and for each distinct state alphabet concept (e.g., DNA), there is
one and only one instance of a representation of that concept (i.e., all
DNA-type data in DendroPy, regardless of the source, refer to the same instance
of the state alphabet and state alphabet elements).
"""

import collections
import itertools
from dendropy.datamodel import base

###############################################################################
## StateAlphabet
class StateAlphabet(
        base.DataObject,
        base.Annotable):

    FUNDAMENTAL_STATE = 0
    AMBIGUOUS_STATE = 1
    POLYMORPHIC_STATE = 2

    """
    A master registry mapping state symbols to their definitions.

    There are two classes or "denominations" of states:

        - fundamental states
            These are the basic, atomic, self-contained states of the alphabet,
            distinct and mutually-exclusive from every other fundamental state.
            E.g., for DNA: adenine, guanine, cytosine, and thymine.

        - multi-state states
            The states are second-level or "pseudo-states", in that they are
            not properly states in and of themselves, but rather each consist
            of a set of other states. That is, a multi-state state is a set of
            two or more fundamental states.
            Multi-state states are of one of two types: "ambiguous" and
            "polymorphic" states. "Ambiguous" states represent states in which
            the true fundamental state is unknown, but consists of one of the
            fundamental states to which the ambiguous states map. "Polymorphic"
            states represent states in which the entity actually has multiple
            fundamental states simultaneously. "Ambiguous" states are an
            expression of uncertainty or lack of knowledge about the identity
            of state. With "polymorphic" states, on the other hand, there is no
            uncertaintly or lack of knowledge about the state: the state is
            known definitively, and it consists of multiple fundamental states.
            An example of an ambiguous state code

    The fundamental states of a state alphabet are, in principle, immutable:
    they are defined at the initialization/construction of a state alphabet,
    and after this, both the set of state definition instances (i.e., the
    particular membership of state instances that make up the fundamental
    states of a given state alphabet) as well as the state definition instances
    themselves (i.e., the definition of each state instance and its attributes,
    such as symbol, index, etc.) are read-only.

    Multi-states of a state alphabet can be defined upon initialization, just
    like fundamental states. And just like fundamental states, the individual
    multi-state definitions cannot be changed or removed from a state alphabet
    (i.e., once defined and added to a state alphabet, a multi-state definition
    cannot be deleted, or can its symbol, etc. be changed). However, new
    multi-state states *can* be added to a state alphabet, as long as the
    symbol of the newly added multi-state do not clash with existing symbols.

    Note that multi-state states can be specified in terms of other multi-state
    states, but that upon instantiation, these member multi-states will be
    expanded to their fundamental states.

    Parameters
    ----------

    label : string, optional
        The name for this state alphabet.

    fundamental_states : iterable of strings
        An iterable of symbols defining the fundamental (i.e., non-ambiguous
        and non-polymorphic states of this alphabet), with a 1-to-1
        correspodence between symbols and states. Each state will also be
        automatically indexed base on its position in this list. For DNA, this
        would be something like: `'ACGT'` or `('A', 'C', 'G', T')`. For
        "standard" characters, this would be something like `'01'` or `('0',
        '1')`.

    ambiguous_states : iterable of tuples
        An iterable consisting of tuples expressing ambiguous state symbols and
        the set of symbols representing the fundamental states to which they
        map. The first element in the tuple is the symbol used to represent
        the ambiguous state; this can be blank (""). The second element is
        an iterable of fundamental state symbols to which this ambiguous
        state maps. The fundamental state symbols *must* have already been
        defined, i.e. given in the value passed to `fundamental_states`.

    polymorphic_states : iterable of tuples
        An iterable consisting of tuples expressing polymorphic state symbols and
        the set of symbols representing the fundamental states to which they
        map. The first element in the tuple is the symbol used to represent
        the polymorphic state; this can be blank (""). The second element is
        an iterable of fundamental state symbols to which this polymorphic
        state maps. The fundamental state symbols *must* have already been
        defined, i.e. given in the value passed to `fundamental_states`.

    symbol_synonyms : dictionary
        A mapping of symbols, with keys being the new symbols and values being
        (already-defined) symbols of states to which they map. This provides a
        mechanism by which states with multiple symbols can be managed. For
        example, an ambiguous state, "unknown", representing all fundamental
        states might be defined with '?' as its primary symbol, and a synonym
        symbol for this state might be 'X'.
    """

    ###########################################################################
    ### Life-Cycle and Identity

    def __init__(self,
            fundamental_states,
            ambiguous_states=None,
            polymorphic_states=None,
            symbol_synonyms=None,
            label=None,
            case_sensitive=True):

        base.DataObject.__init__(self, label=label)
        self._is_case_sensitive = case_sensitive

        # Core collection underlying alphabet
        self._fundamental_states = []
        self._ambiguous_states = []
        self._polymorphic_states = []

        # Cache invalidation flag
        self._is_dirty = True

        # Populate core collection
        for symbol in fundamental_states:
            self.new_fundamental_state(symbol)
        if ambiguous_states:
            for ss in ambiguous_states:
                self.new_ambiguous_state(symbol=ss[0], member_state_symbols=ss[1])
        if polymorphic_states:
            for ss in polymorphic_states:
                self.new_polymorphic_state(symbol=ss[0], member_state_symbols=ss[1])
        if symbol_synonyms:
            for k in symbol_synonyms:
                self.new_symbol_synonym(k, symbol_synonyms[v])

        # Build mappings
        self.compile_lookup_mappings()


    def __copy__(self, memo=None):
        raise TypeError("Cannot (shallow) copy {}".format(self.__class__.__name__))

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot (shallow) copy {}".format(self.__class__.__name__))

    def __deepcopy__(self, memo=None):
        return base.Annotable.__deepcopy__(self, memo=memo)

    ###########################################################################
    ### Symbol Management

    def _direct_get_state_for_symbol(self, symbol):
        """
        Returns the :class:`StateLetter` instance corresponding to `symbol`.
        """
        for state_symbol, state in self.symbol_state_pair_iter(include_synonyms=True):
            if state_symbol == symbol:
                return state
        raise KeyError(symbol)

    def _direct_get_state_set_for_symbols(self, symbols):
        """
        Returns the list of :class:`StateLetter` instances corresponding to
        the iterable of symbols given by `symbols`, with each element in
        `symbols` corresponding to a single symbol.
        """
        ss = []
        for symbol in symbols:
            state = self._direct_get_state_for_symbol(symbol)
            ss.extend(state.fundamental_states)
        return ss

    def _validate_new_symbol(self, symbol):
        if symbol is None or symbol == "":
            raise ValueError("Cannot validate empty symbol")
        symbol = str(symbol)
        for state_symbol, state in self.symbol_state_pair_iter(include_synonyms=True):
            if state_symbol == symbol:
                raise ValueError("State with symbol or symbol synonym of '{}' already defined is this alphabet".format(symbol))
            if not self._is_case_sensitive:
                if state_symbol.upper() == symbol.upper():
                    raise ValueError("State with symbol or symbol synonym of '{}' already defined is this alphabet".format(symbol))
        return symbol

    def new_fundamental_state(self, symbol):
        """
        Adds a new fundamental state to the collection
        of states.
        """
        if symbol is None or symbol == "":
            raise ValueError("Fundamental states cannot be defined without a valid symbol")
        symbol = self._validate_new_symbol(symbol)
        index = len(self._fundamental_states)
        new_state = StateLetter(
                symbol=symbol,
                index=index,
                state_denomination=StateAlphabet.FUNDAMENTAL_STATE,
                member_states=None)
        self._fundamental_states.append(new_state)
        self._is_dirty = True
        return new_state

    def new_ambiguous_state(self,
            symbol,
            member_state_symbols):
        if symbol is not None and symbol != "":
            symbol = self._validate_new_symbol(symbol)
        member_states = self._direct_get_state_set_for_symbols(member_state_symbols)
        member_states = frozenset(member_states)
        new_state = StateLetter(
                symbol=symbol,
                index=None,
                state_denomination=StateAlphabet.AMBIGUOUS_STATE,
                member_states=member_states)
        self._ambiguous_states.append(new_state)
        self._is_dirty = True
        return new_state

    def new_polymorphic_state(self,
            symbol,
            member_state_symbols):
        if symbol is not None and symbol != "":
            symbol = self._validate_new_symbol(symbol)
        member_states = self._direct_get_state_set_for_symbols(member_state_symbols)
        member_states = frozenset(member_states)
        new_state = StateLetter(
                symbol=symbol,
                index=None,
                state_denomination=StateAlphabet.POLYMORPHIC_STATE,
                member_states=member_states)
        self._polymorphic_states.append(new_state)
        self._is_dirty = True
        return new_state

    def new_symbol_synonym(self,
            symbol_synonym, referenced_symbol):
        if symbol is None or symbol == "":
            raise ValueError("Symbol synonym cannot be empty")
        if symbol is not None and symbol != "":
            symbol = self._validate_new_symbol(symbol)
        state = self._direct_get_state_for_symbol(referenced_symbol)
        state.symbol_synonyms.add(symbol_synonym)
        self._is_dirty = True
        return state

    def _set_symbol_mapping(self, d, symbol, state):
        if symbol is None or symbol == "":
            raise ValueError("Symbol synonym cannot be empty")
        assert symbol not in d
        d[symbol] = state
        if not self._is_case_sensitive:
            for s in (symbol.lower(), symbol.upper()):
                if s != symbol:
                    assert s not in d
                    d[s] = state

    def compile_lookup_mappings(self):
        self._canonical_symbol_state_map = collections.OrderedDict()
        self._full_symbol_state_map = collections.OrderedDict()
        self._index_state_map = collections.OrderedDict()
        self._fundamental_states_to_ambiguous_state_map = {}
        self._fundamental_states_to_polymorphic_state_map = {}
        for idx, state in enumerate(self.state_iter()):
            if state.symbol:
                self._set_symbol_mapping(
                        self._canonical_symbol_state_map,
                        state.symbol,
                        state)
                self._set_symbol_mapping(
                        self._full_symbol_state_map,
                        state.symbol,
                        state)
                if state.symbol_synonyms:
                    for ss in state.symbol_synonyms:
                        self._set_symbol_mapping(
                                self._full_symbol_state_map,
                                ss,
                                state)
            else:
                assert state.state_denomination != StateAlphabet.FUNDAMENTAL_STATE
            # if state in self._fundamental_states:
            #     assert idx == state._index
            # else:
            #     state._index = idx
            state._index = idx
            self._index_state_map[idx] = state
            if state.state_denomination == StateAlphabet.AMBIGUOUS_STATE:
                member_states = state.member_states
                if member_states in self._fundamental_states_to_ambiguous_state_map:
                    raise ValueError("Multiple definitions of ambiguous state with member states of '{}': {}, {}. Define a symbol synonym instead.".format(
                        state.member_states_str, self._fundamental_states_to_ambiguous_state_map[member_states], state))
                assert member_states not in self._fundamental_states_to_ambiguous_state_map
                self._fundamental_states_to_ambiguous_state_map[member_states] = state
            elif state.state_denomination == StateAlphabet.POLYMORPHIC_STATE:
                member_states = state.member_states
                if member_states in self._fundamental_states_to_polymorphic_state_map:
                    raise ValueError("Multiple definitions of polymorphic state with member states of '{}': {}, {}. Define a symbol synonym instead.".format(
                        state.member_states_str, self._fundamental_states_to_polymorphic_state_map[member_states], state))
                self._fundamental_states_to_polymorphic_state_map[member_states] = state
        self._is_dirty = False

    ###########################################################################
    ### Symbol Access

    def state_iter(self):
        """
        Returns a "raw" iterator over all state identities.
        """
        return itertools.chain(
                self._fundamental_states,
                self._ambiguous_states,
                self._polymorphic_states)

    def symbol_state_pair_iter(self, include_synonyms=True):
        """
        Returns a "raw" iterator over pairs of symbols and the states to which
        they correspond.
        """
        for state in self.state_iter():
            yield (state.symbol, state)
            if include_synonyms:
                for synonym in state.symbol_synonyms:
                    yield (synonym, state)

    def _get_canonical_symbol_state_map(self):
        if self._is_dirty:
            self.compile_lookup_mappings()
        return self._canonical_symbol_state_map
    canonical_symbol_state_map = property(_get_canonical_symbol_state_map)

    def _get_full_symbol_state_map(self):
        if self._is_dirty:
            self.compile_lookup_mappings()
        return self._full_symbol_state_map
    full_symbol_state_map = property(_get_full_symbol_state_map)



###############################################################################
## StateLetter

class StateLetter(
        base.DataObject,
        base.Annotable):
    """
    A character state definition, which can either be a fundamental state or
    a mapping to a set of other character states (for polymorphic or ambiguous
    characters).
    """

    def __init__(self,
            symbol=None,
            index=None,
            state_denomination=StateAlphabet.FUNDAMENTAL_STATE,
            member_states=None):
        """
        A state is immutable with respect to its definition and identity.
        Specifically, it 'symbol', 'index', 'multistate', and 'member_states'
        properties are set upon definition/creation, and after that are
        read-only.

        Parameters
        ----------
        symbol : string
            A text symbol or token representation of this character state.
            E.g., 'G' for the base guanine in a DNA state alphabet, or '1' for
            presence of a wing in a morphological data set.
        index : integer
            The (0-based) numeric index for this state in the state alphabet.
            E.g., for a DNA alphabet: 0 = 'A'/adenine, 1 = 'C'/cytosine, 2 =
            'G'/guanine, 3 = 'T'/thymine. Or for a "standard" alphabet: 0 =
            '0', 1 = '1'. Note that ambiguous and polymorphic state definitions
            typically are not indexed.
        state_denomination : 'enum'
            One of: `StateAlphabet.FUNDAMENTAL_STATE`,
            `StateAlphabet.AMBIGUOUS_STATE`, or
            `StateAlphabet.POLYMORPHIC_STATE`.
        member_states : iterable of :class:`StateLetter` instances.
            If a multi-state, then a collection of :class:`StateLetter`
            instances to which this state maps.
        """
        base.DataObject.__init__(self, label=None)
        self._symbol = symbol
        self._index = index
        self._state_denomination = state_denomination
        self._member_states = None
        self._fundamental_states = None
        self._fundamental_symbols = None
        self._fundamental_indexes = None
        self._partials_vector = None
        self._member_states = member_states
        self.symbol_synonyms = set()

    def __str__(self):
        if self._symbol:
            return str(self._symbol)
        elif self._state_denomination == StateAlphabet.FUNDAMENTAL_STATE:
            return ""
        else:
            return self.member_states_str

    def __repr__(self):
        s = str(self)
        return "<{} at {}: '{}'>".format(self.__class__.__name__,
                hex(id(self)), str(s))

    def _get_member_states_str(self):
        s = ",".join([m._symbol for m in self._member_states])
        if self._state_denomination == StateAlphabet.POLYMORPHIC_STATE:
            return "(" + s + ")"
        else:
            return "{" + s + "}"
    member_states_str = property(_get_member_states_str)

    def _get_symbol(self):
        return self._symbol
    symbol = property(_get_symbol)

    def _get_state_denomination(self):
        return self._state_denomination
    state_denomination = property(_get_state_denomination)

    def _is_single_state(self):
        return self._state_denomination == StateAlphabet.FUNDAMENTAL_STATE
    is_single_state = property(_is_single_state)

    def _is_fundamental_state(self):
        return self._state_denomination == StateAlphabet.FUNDAMENTAL_STATE
    is_fundamental_state = property(_is_fundamental_state)

    def _get_member_states(self):
        return self._member_states
    def _set_member_states(self, m):
        self._member_states = m
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

    def _get_fundamental_indexes(self):
        "Returns set of indexes of all _get_fundamental states to which this state maps."
        if self._fundamental_indexes is None:
            self._fundamental_indexes = set([state.index for state in self.fundamental_states])
        return self._fundamental_indexes
    fundamental_indexes = property(_get_fundamental_indexes)

    def is_exact_correspondence(self, other):
        """
        Tries to determine if two StateAlphabetElement definitions
        are equivalent by matching symbols.
        """
        match = True
        if self._state_denomination != other._state_denomination:
            return False
        if self._state_denomination != StateAlphabet.FUNDAMENTAL_STATE and other._state_denomination != StateAlphabet.FUNDAMENTAL_STATE:
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
            return self._symbol == other._symbol

if __name__ == "__main__":

    fundamental_states = "ACGT-"
    ambiguous_states = (
            ("?", ('A', 'C', 'G', 'T', '-')),
            ("N", ('A', 'C', 'G', 'T')),
            ("M", ('A', 'C')),
            ("R", ('A', 'G')),
            ("W", ('A', 'T')),
            ("S", ('C', 'G')),
            ("Y", ('C', 'T')),
            ("K", ('G', 'T')),
            ("V", ('A', 'C', 'G')),
            ("H", ('A', 'C', 'T')),
            ("D", ('A', 'G', 'T')),
            ("B", ('C', 'G', 'T')),
            )
    polymorphic_states = (
            ("ZZZ", ("S", "Y")),
            ("", ("A", "C", "G", "T")),
            )
    symbol_synonyms = {
            "U": "T",
            "X": "N",
            }
    unknown_state_symbol = 'N'

    dna = StateAlphabet(
        fundamental_states=fundamental_states,
        ambiguous_states=ambiguous_states,
        polymorphic_states=polymorphic_states,
        )
    for s in dna.state_iter():
        print(s, s.symbol, s.member_states)
    for s in dna.canonical_symbol_state_map:
        print(s, dna.canonical_symbol_state_map[s])


