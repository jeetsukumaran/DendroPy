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
from dendropy.datamodel import basemodel
from dendropy.utility import container

###############################################################################
## StateAlphabet

class StateAlphabet(
        basemodel.DataObject,
        basemodel.Annotable):
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
            An example of an ambiguous state would be 'N', representing any
            base in molecular sequence data. An example of a polymorphic state
            would be the range of a widespread species found in multiple
            geographic units.

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
        map. The first element in the tuple is the symbol used to represent the
        ambiguous state; this can be blank (""), but if not blank it needs to
        be unique across all symbols (including case-variants if the state
        alphabet is case-insensitive).  The second element is an
        iterable of fundamental state symbols to which this ambiguous state
        maps. The fundamental state symbols *must* have already been defined,
        i.e. given in the value passed to `fundamental_states`. Note: a
        dictionary may seem like a more tractable structure than iterable of
        tuples, but we may need to specify multiple anonymous or blank
        ambiguous states.

    polymorphic_states : iterable of tuples
        An iterable consisting of tuples expressing polymorphic state symbols and
        the set of symbols representing the fundamental states to which they
        map. The first element in the tuple is the symbol used to represent the
        polymorphic state; this can be blank (""), but if not blank it needs to
        be unique across all symbols (including case-variants if the state
        alphabet is case-insensitive).  The second element is an
        iterable of fundamental state symbols to which this polymorphic state
        maps. The fundamental state symbols *must* have already been defined,
        i.e. given in the value passed to `fundamental_states`. Note: a
        dictionary may seem like a more tractable structure than iterable of
        tuples, but we may need to specify multiple anonymous or blank
        polymorphic states.

    symbol_synonyms : dictionary
        A mapping of symbols, with keys being the new symbols and values being
        (already-defined) symbols of states to which they map. This provides a
        mechanism by which states with multiple symbols can be managed. For
        example, an ambiguous state, "unknown", representing all fundamental
        states might be defined with '?' as its primary symbol, and a synonym
        symbol for this state might be 'X'.
    """

    ###########################################################################
    ### CLass-level Constants

    FUNDAMENTAL_STATE = 0
    AMBIGUOUS_STATE = 1
    POLYMORPHIC_STATE = 2

    ###########################################################################
    ### Life-Cycle and Identity

    def __init__(self,
            fundamental_states=None,
            ambiguous_states=None,
            polymorphic_states=None,
            symbol_synonyms=None,
            label=None,
            case_sensitive=True):

        basemodel.DataObject.__init__(self, label=label)
        self._is_case_sensitive = case_sensitive

        # Core collection underlying alphabet
        self._fundamental_states = []
        self._ambiguous_states = []
        self._polymorphic_states = []

        # Look-up mappings
        self._canonical_symbol_state_map = None
        self._full_symbol_state_map = None
        self._index_state_map = None
        self._fundamental_states_to_ambiguous_state_map = None
        self._fundamental_states_to_polymorphic_state_map = None

        # Cache invalidation flag
        self._is_dirty = True

        # Populate core collection
        if fundamental_states:
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
                    self.new_symbol_synonym(k, symbol_synonyms[k])
            # Build mappings
            self.compile_lookup_mappings()

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return other is self

    def __copy__(self, memo=None):
        return self

    def taxon_namespace_scoped_copy(self, memo=None):
        return self

    def __deepcopy__(self, memo=None):
        return basemodel.Annotable.__deepcopy__(self, memo=memo)

    ###########################################################################
    ### Symbol Management

    def _direct_get_state_for_symbol(self, symbol):
        """
        Returns the :class:`StateIdentity` instance corresponding to `symbol`.
        """
        for state_symbol, state in self.symbol_state_pair_iter(include_synonyms=True):
            if state_symbol == symbol:
                return state
        raise KeyError(symbol)

    def _direct_get_fundamental_state_set_for_symbols_for_symbols(self, symbols):
        """
        Returns the list of :class:`StateIdentity` instances corresponding to
        the iterable of symbols given by `symbols`, with each element in
        `symbols` corresponding to a single symbol.
        """
        ss = []
        for symbol in symbols:
            state = self._direct_get_state_for_symbol(symbol)
            ss.extend(state.fundamental_states)
        return frozenset(ss)

    def _validate_new_symbol(self, symbol):
        if symbol is None or symbol == "":
            raise ValueError("Cannot validate empty symbol")
        symbol = str(symbol)
        for state_symbol, state in self.symbol_state_pair_iter(include_synonyms=True):
            if state_symbol == symbol:
                raise ValueError("State with symbol or symbol synonym of '{}' already defined is this alphabet".format(symbol))
        return symbol

    def new_fundamental_state(self, symbol):
        """
        Adds a new fundamental state to the collection
        of states in this alphabet.

        Parameters
        ----------
        symbol : string
            The symbol used to represent this state. Cannot have previously
            been used to refer to any other state, fundamental or otherwise, as
            a primary or synonymous symbol (including implicit synonyms given
            by case-variants if the state alphabet is not case-sensitive).
            Cannot be blank ("") or `None`.

        Returns
        -------
        s : :class:`StateIdentity`
            The new state created and added.
        """
        if symbol is None or symbol == "":
            raise ValueError("Fundamental states cannot be defined without a valid symbol")
        symbol = self._validate_new_symbol(symbol)
        index = len(self._fundamental_states)
        new_state = StateIdentity(
                symbol=symbol,
                index=index,
                state_denomination=StateAlphabet.FUNDAMENTAL_STATE,
                member_states=None)
        self._fundamental_states.append(new_state)
        if not self._is_case_sensitive:
            for s in (symbol.upper(), symbol.lower()):
                if s != symbol:
                    self.new_symbol_synonym(s, symbol)
        self._is_dirty = True
        return new_state

    def new_ambiguous_state(self,
            symbol,
            member_state_symbols):
        """
        Adds a new ambiguous state to the collection
        of states in this alphabet.

        Parameters
        ----------
        symbol : string or None
            The symbol used to represent this state. Cannot have previously
            been used to refer to any other state, fundamental or otherwise, as
            a primary or synonymous symbol (including implicit synonyms given
            by case-variants if the state alphabet is not case-sensitive). Can
            be blank ("") or `None` if there.

        member_states : iterable of strings
            List of symbols representing states to which this state maps. Symbols
            representing multistates will taken to refer to the set of
            fundamental states to which they, in turn, map.

        Returns
        -------
        s : :class:`StateIdentity`
            The new state created and added.
        """
        if symbol is not None and symbol != "":
            symbol = self._validate_new_symbol(symbol)
        member_states = self._direct_get_fundamental_state_set_for_symbols_for_symbols(member_state_symbols)
        new_state = StateIdentity(
                symbol=symbol,
                index=None,
                state_denomination=StateAlphabet.AMBIGUOUS_STATE,
                member_states=member_states)
        self._ambiguous_states.append(new_state)
        if symbol and not self._is_case_sensitive:
            for s in (symbol.upper(), symbol.lower()):
                if s != symbol:
                    self.new_symbol_synonym(s, symbol)
        self._is_dirty = True
        return new_state

    def new_polymorphic_state(self,
            symbol,
            member_state_symbols):
        """
        Adds a new polymorphic state to the collection
        of states in this alphabet.

        Parameters
        ----------
        symbol : string or None
            The symbol used to represent this state. Cannot have previously
            been used to refer to any other state, fundamental or otherwise, as
            a primary or synonymous symbol (including implicit synonyms given
            by case-variants if the state alphabet is not case-sensitive). Can
            be blank ("") or `None` if there.

        member_states : iterable of strings
            List of symbols representing states to which this state maps. Symbols
            representing multistates will taken to refer to the set of
            fundamental states to which they, in turn, map.

        Returns
        -------
        s : :class:`StateIdentity`
            The new state created and added.
        """
        if symbol is not None and symbol != "":
            symbol = self._validate_new_symbol(symbol)
        member_states = self._direct_get_fundamental_state_set_for_symbols_for_symbols(member_state_symbols)
        new_state = StateIdentity(
                symbol=symbol,
                index=None,
                state_denomination=StateAlphabet.POLYMORPHIC_STATE,
                member_states=member_states)
        self._polymorphic_states.append(new_state)
        if symbol and not self._is_case_sensitive:
            for s in (symbol.upper(), symbol.lower()):
                if s != symbol:
                    self.new_symbol_synonym(s, symbol)
        self._is_dirty = True
        return new_state

    def new_symbol_synonym(self,
            symbol_synonym, referenced_symbol):
        """
        Defines an alternative symbol mapping for an existing state.

        Parameters
        ----------
        symbol_synonym : string
            The (new) alternative symbol.

        referenced_symbol : string
            The symbol for the state to which the alternative symbol will also
            map.

        Returns
        -------
        s : :class:`StateIdentity`
            The state to which this synonym maps.
        ------
        """
        if symbol_synonym is None or symbol_synonym == "":
            raise ValueError("Symbol synonym cannot be empty")
        symbol_synonym = self._validate_new_symbol(symbol_synonym)
        state = self._direct_get_state_for_symbol(referenced_symbol)
        state.symbol_synonyms.add(symbol_synonym)
        self._is_dirty = True
        return state

    def _set_symbol_mapping(self, d, symbol, state):
        if symbol is None or symbol == "":
            raise ValueError("Symbol synonym cannot be empty")
        assert symbol not in d
        d[symbol] = state

    def compile_lookup_mappings(self):
        """
        Builds lookup tables/mappings for quick referencing and dereferencing
        of symbols/states.
        """
        temp_canonical_symbol_state_map = collections.OrderedDict()
        temp_full_symbol_state_map = collections.OrderedDict()
        temp_index_state_map = collections.OrderedDict()
        temp_fundamental_states_to_ambiguous_state_map = {}
        temp_fundamental_states_to_polymorphic_state_map = {}
        for idx, state in enumerate(self.state_iter()):
            if state.symbol:
                assert state.symbol not in temp_canonical_symbol_state_map
                temp_canonical_symbol_state_map[state.symbol] = state
                self._set_symbol_mapping(
                        temp_full_symbol_state_map,
                        state.symbol,
                        state)
                if state.symbol_synonyms:
                    for ss in state.symbol_synonyms:
                        self._set_symbol_mapping(
                                temp_full_symbol_state_map,
                                ss,
                                state)
            else:
                assert state.state_denomination != StateAlphabet.FUNDAMENTAL_STATE
            # if state in temp_fundamental_states:
            #     assert idx == state._index
            # else:
            #     state._index = idx
            state._index = idx
            temp_index_state_map[idx] = state
            if state.state_denomination == StateAlphabet.AMBIGUOUS_STATE:
                member_states = state.member_states
                if member_states in temp_fundamental_states_to_ambiguous_state_map:
                    raise ValueError("Multiple definitions of ambiguous state with member states of '{}': {}, {}. Define a symbol synonym instead.".format(
                        state.member_states_str, temp_fundamental_states_to_ambiguous_state_map[member_states], state))
                assert member_states not in temp_fundamental_states_to_ambiguous_state_map
                temp_fundamental_states_to_ambiguous_state_map[member_states] = state
            elif state.state_denomination == StateAlphabet.POLYMORPHIC_STATE:
                member_states = state.member_states
                if member_states in temp_fundamental_states_to_polymorphic_state_map:
                    raise ValueError("Multiple definitions of polymorphic state with member states of '{}': {}, {}. Define a symbol synonym instead.".format(
                        state.member_states_str, temp_fundamental_states_to_polymorphic_state_map[member_states], state))
                temp_fundamental_states_to_polymorphic_state_map[member_states] = state
        self._canonical_symbol_state_map = container.FrozenOrderedDict(temp_canonical_symbol_state_map)
        self._full_symbol_state_map = container.FrozenOrderedDict(temp_full_symbol_state_map)
        self._index_state_map = container.FrozenOrderedDict(temp_index_state_map)
        self._fundamental_states_to_ambiguous_state_map = container.FrozenOrderedDict(temp_fundamental_states_to_ambiguous_state_map)
        self._fundamental_states_to_polymorphic_state_map = container.FrozenOrderedDict(temp_fundamental_states_to_polymorphic_state_map)
        temp_is_dirty = False

    def set_state_as_attribute(self, state, attr_name=None):
        """
        Sets the given state as an attribute of this alphabet.
        The name of the attribute will be `attr_name` if specified,
        or the state symbol otherwise.

        Parameters
        ----------
        state : :class:`StateIdentity`
            The state to be made an attribute of this alphabet.
        attr_name : string
            The name of the attribute. If not specified, the state
            symbol will be used.
        """
        if (state not in self._fundamental_states
                and state not in self._ambiguous_states
                and state not in self._polymorphic_states):
            raise ValueError("State {} not defined in current alphabet".format(state))
        if attr_name is None:
            attr_name = state.symbol
        if attr_name is None:
            raise TypeError("Cannot set attribute: non-None symbol needed for state or non-None attribute name needs to be provided")
        setattr(self, attr_name, state)

    ###########################################################################
    ### Symbol Access

    def __len__(self):
        """
        Number of states.
        """
        return ( len(self._fundamental_states)
                + len(self._ambiguous_states)
                + len(self._polymorphic_states) )

    def state_iter(self):
        """
        Returns an iterator over all state identities.
        """
        return itertools.chain(
                self._fundamental_states,
                self._ambiguous_states,
                self._polymorphic_states)

    def fundamental_state_iter(self):
        """
        Returns an iterator over all fundamental state identities.
        """
        return itertools.chain(self._fundamental_states)

    def ambiguous_state_iter(self):
        """
        Returns an iterator over all ambiguous state identities.
        """
        return itertools.chain(self._ambiguous_states)

    def polymorphic_state_iter(self):
        """
        Returns an iterator over all polymorphic state identities.
        """
        return itertools.chain(self._polymorphic_states)

    def multistate_state_iter(self):
        """
        Returns an iterator over all ambiguous and polymorphic state
        identities.
        """
        return itertools.chain(self._ambiguous_states, self._polymorphic_states)

    def fundamental_symbol_iter(self, include_synonyms=True):
        """
        Returns an iterator over all symbols (including synonyms, unless
        `include_synonyms` is `False`) that map to fundamental states.
        """
        for state in self.fundamental_state_iter():
            yield state.symbol
            if state.symbol_synonyms and include_synonyms:
                for symbol in state.symbol_synonyms:
                    yield symbol

    def ambiguous_symbol_iter(self, include_synonyms=True):
        """
        Returns an iterator over all symbols (including synonyms, unless
        `include_synonyms` is `False`) that map to ambiguous states.
        """
        for state in self.ambiguous_state_iter():
            yield state.symbol
            if state.symbol_synonyms and include_synonyms:
                for symbol in state.symbol_synonyms:
                    yield symbol

    def polymorphic_symbol_iter(self, include_synonyms=True):
        """
        Returns an iterator over all symbols (including synonyms, unless
        `include_synonyms` is `False`) that map to polymorphic states.
        """
        for state in self.polymorphic_state_iter():
            yield state.symbol
            if state.symbol_synonyms and include_synonyms:
                for symbol in state.symbol_synonyms:
                    yield symbol

    def multistate_symbol_iter(self, include_synonyms=True):
        """
        Returns an iterator over all symbols (including synonyms, unless
        `include_synonyms` is `False`) that map to multistate states.
        """
        for state in self.multistate_state_iter():
            yield state.symbol
            if state.symbol_synonyms and include_synonyms:
                for symbol in state.symbol_synonyms:
                    yield symbol

    def symbol_state_pair_iter(self, include_synonyms=True):
        """
        Returns an iterator over all symbols paired with the state to which the
        they symbols map.
        """
        for state in self.state_iter():
            yield (state.symbol, state)
            if include_synonyms:
                for synonym in state.symbol_synonyms:
                    yield (synonym, state)

    def _get_canonical_symbol_state_map(self):
        """
        Dictionary with state symbols as keys and states as values. Does not
        include symbol synonyms or case variations.
        """
        if self._is_dirty:
            self.compile_lookup_mappings()
        return self._canonical_symbol_state_map
    canonical_symbol_state_map = property(_get_canonical_symbol_state_map, __doc__)

    def _get_full_symbol_state_map(self):
        """
        Dictionary with state symbols as keys and states as values.
        Includes symbol synonyms or case variations.
        """
        if self._is_dirty:
            self.compile_lookup_mappings()
        return self._full_symbol_state_map
    full_symbol_state_map = property(_get_full_symbol_state_map, __doc__)

    def __getitem__(self, key):
        """
        Returns state identity corresponding to `key`.

        Parameters
        ----------
        key : integer or string
            If and integer value, looks up and returns state identity by index.
            If a string value, looks up and returns state identity by symbol.

        Returns
        -------
        s : :class:`StateIdentity` instance
            Returns a :class:`StateIdentity` corresponding to `key`.

        Raises
        ------
        KeyError if `key` is not valid.

        """
        if self._is_dirty:
            self.compile_lookup_mappings()
        if isinstance(key, int):
            return self._index_state_map[key]
        else:
            return self._full_symbol_state_map[key]

    def get_states_for_symbols(self, symbols):
        """
        Returns list of states corresponding to symbols.

        Parameters
        ----------
        symbols : iterable of symbols

        Returns
        -------
        s : list of :class:`StateIdentity`
            A list of :class:`StateIdentity` instances corresponding to symbols
            given in `symbols`.
        """
        states = [self.full_symbol_state_map[s] for s in symbols]
        return states

    def get_fundamental_states_for_symbols(self, symbols):
        """
        Returns list of *fundamental* states corresponding to symbols.

        Parameters
        ----------
        symbols : iterable of symbols

        Returns
        -------
        s : list of :class:`StateIdentity`
            A list of fundamental :class:`StateIdentity` instances corresponding
            to symbols given in `symbols`, with multi-state states expanded
            into their fundamental symbols.
        """
        states = []
        for symbol in symbols:
            state = self._full_symbol_state_map[symbol]
            states.extend(state.fundamental_states)
        return states

    def get_canonical_symbol_for_symbol(self, symbol):
        """
        Returns the canonical state symbol for the state to which `symbol`
        maps. E.g., in a DNA alphabet, return 'A' for 'a'.

        Parameters
        ----------
        symbol : string

        Returns
        -------
        s : string
            Canonical symbol for state with symbol or synonym symbol of
            `symbol`.
        """
        return self[symbol].symbol

    def match_ambiguous_state(self, symbols):
        """
        Returns ambiguous state with fundamental member states
        represented by symbols given in `symbols`.

        Parameters
        ----------
        symbols : iterable of symbols

        Returns
        -------
        s : :class:`StateIdentity` instance
        """
        states = frozenset(self.get_fundamental_states_for_symbols(symbols))
        return self._fundamental_states_to_ambiguous_state_map[states]

    def match_polymorphic_state(self, symbols):
        """
        Returns polymorphic state with fundamental member states
        represented by symbols given in `symbols`.

        Parameters
        ----------
        symbols : iterable of symbols

        Returns
        -------
        s : :class:`StateIdentity` instance
        """
        states = frozenset(self.get_fundamental_states_for_symbols(symbols))
        return self._fundamental_states_to_polymorphic_state_map[states]

    def match_state(self, symbols):
        """
        No longer supported: use :meth:`StateAlphabet.match_ambiguous_state()` or
        :meth:`StateAlphabet.match_polymorphic_state()` instead.
        """
        raise NotImplementedError(__doc__)

###############################################################################
## StateIdentity

class StateIdentity(
        basemodel.DataObject,
        basemodel.Annotable):
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
        member_states : iterable of :class:`StateIdentity` instances.
            If a multi-state, then a collection of :class:`StateIdentity`
            instances to which this state maps.
        """
        basemodel.DataObject.__init__(self, label=symbol)
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

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return other is self

    def __copy__(self, memo=None):
        return self

    def taxon_namespace_scoped_copy(self, memo=None):
        return self

    def __deepcopy__(self, memo=None):
        return basemodel.Annotable.__deepcopy__(self, memo=memo)

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
        Tries to determine if two StateIdentity definitions
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

###############################################################################
## DnaStateAlphabet

class DnaStateAlphabet(StateAlphabet):

    def __init__(self):
        fundamental_states = "ACGT-"
        polymorphic_states = None
        ambiguous_states = (
                ("?", "ACGT-"),
                ("N", "ACGT"),
                ("R", "AG"  ),
                ("Y", "CT"  ),
                ("M", "AC"  ),
                ("W", "AT"  ),
                ("S", "CG"  ),
                ("K", "GT"  ),
                ("V", "ACG" ),
                ("H", "ACT" ),
                ("D", "AGT" ),
                ("B", "CGT" ),
                )
        symbol_synonyms = {"X": "N"}
        StateAlphabet.__init__(self,
                fundamental_states=fundamental_states,
                polymorphic_states=polymorphic_states,
                ambiguous_states=ambiguous_states,
                symbol_synonyms=symbol_synonyms,
                label="DNA",
                case_sensitive=False)
        for state in self.state_iter():
            if state.symbol == "-":
                attr_name = "gap"
            elif state.symbol == "?":
                attr_name = "missing"
            else:
                attr_name = state.symbol
            self.set_state_as_attribute(state, attr_name)
        self.any_residue = self.N
        self.unknown_state_symbol = 'N'

###############################################################################
## RnaStateAlphabet

class RnaStateAlphabet(StateAlphabet):

    def __init__(self):
        fundamental_states = "ACGU-"
        polymorphic_states = None
        ambiguous_states = (
                ("?", "ACGU-"),
                ("N", "ACGU"),
                ("R", "AG"  ),
                ("Y", "CU"  ),
                ("M", "AC"  ),
                ("W", "AU"  ),
                ("S", "CG"  ),
                ("K", "GU"  ),
                ("V", "ACG" ),
                ("H", "ACU" ),
                ("D", "AGU" ),
                ("B", "CGU" ),
                )
        symbol_synonyms = {"X": "N"}
        StateAlphabet.__init__(self,
                fundamental_states=fundamental_states,
                polymorphic_states=polymorphic_states,
                ambiguous_states=ambiguous_states,
                symbol_synonyms=symbol_synonyms,
                label="RNA",
                case_sensitive=False)
        for state in self.state_iter():
            if state.symbol == "-":
                attr_name = "gap"
            elif state.symbol == "?":
                attr_name = "missing"
            else:
                attr_name = state.symbol
            self.set_state_as_attribute(state, attr_name)
        self.any_residue = self.N
        self.unknown_state_symbol = 'N'

###############################################################################
## NucleotideStateAlphabet

class NucleotideStateAlphabet(StateAlphabet):

    def __init__(self):
        fundamental_states = "ACGTU-"
        polymorphic_states = None
        ambiguous_states = (
                ("?", "ACGTU-"),
                ("N", "ACGTU"),
                ("R", "AG"  ),
                ("Y", "CTU"  ),
                ("M", "AC"  ),
                ("W", "ATU"  ),
                ("S", "CG"  ),
                ("K", "GTU"  ),
                ("V", "ACG" ),
                ("H", "ACTU" ),
                ("D", "AGTU" ),
                ("B", "CGTU" ),
                )
        symbol_synonyms = {"X": "N"}
        StateAlphabet.__init__(self,
                fundamental_states=fundamental_states,
                polymorphic_states=polymorphic_states,
                ambiguous_states=ambiguous_states,
                symbol_synonyms=symbol_synonyms,
                label="Nucleotide",
                case_sensitive=False)
        for state in self.state_iter():
            if state.symbol == "-":
                attr_name = "gap"
            elif state.symbol == "?":
                attr_name = "missing"
            else:
                attr_name = state.symbol
            self.set_state_as_attribute(state, attr_name)
        self.any_residue = self.N
        self.unknown_state_symbol = 'N'

###############################################################################
## ProteinStateAlphabet

class ProteinStateAlphabet(StateAlphabet):

    def __init__(self):
        fundamental_states = "ACDEFGHIKLMNPQRSTUVWY*-"
        polymorphic_states = None
        ambiguous_states = (
                ("B", "DN"),
                ("Z", "EQ"),
                ("X", "ACDEFGHIKLMNPQRSTUVWY*"),
                ("?", "ACDEFGHIKLMNPQRSTUVWY*-"),
                )
        symbol_synonyms = {}
        StateAlphabet.__init__(self,
                fundamental_states=fundamental_states,
                polymorphic_states=polymorphic_states,
                ambiguous_states=ambiguous_states,
                symbol_synonyms=symbol_synonyms,
                label="Protein",
                case_sensitive=False)
        for state in self.state_iter():
            if state.symbol == "-":
                attr_name = "gap"
            elif state.symbol == "?":
                attr_name = "missing"
            elif state.symbol == "*":
                attr_name = "stop"
            else:
                attr_name = state.symbol
            self.set_state_as_attribute(state, attr_name)
        self.any_residue = self.X
        self.unknown_state_symbol = 'X'

###############################################################################
## BinaryStateAlphabet

class BinaryStateAlphabet(StateAlphabet):

    def __init__(self, allow_gaps=False, allow_missing=False):
        fundamental_states = "10"
        if allow_gaps:
            fundamental_states += "-"
        polymorphic_states = None
        ambiguous_states = []
        if allow_missing:
            ambiguous_states.append( ("?", fundamental_states) )
        symbol_synonyms = {}
        StateAlphabet.__init__(self,
                fundamental_states=fundamental_states,
                polymorphic_states=polymorphic_states,
                ambiguous_states=ambiguous_states,
                symbol_synonyms=symbol_synonyms,
                label="Binary",
                case_sensitive=False)
        for state in self.state_iter():
            if state.symbol == "-":
                attr_name = "gap"
            elif state.symbol == "?":
                attr_name = "missing"
            elif state.symbol == "*":
                attr_name = "stop"
            else:
                attr_name = state.symbol
            self.set_state_as_attribute(state, attr_name)

###############################################################################
## RestrictionSitesStateAlphabet

class RestrictionSitesStateAlphabet(BinaryStateAlphabet):

    def __init__(self, allow_gaps=False, allow_missing=False):
        BinaryStateAlphabet.__init__(self, allow_gaps=allow_gaps, allow_missing=allow_missing)

###############################################################################
## InfiniteSitesStateAlphabet

class InfiniteSitesStateAlphabet(BinaryStateAlphabet):

    def __init__(self, allow_gaps=False, allow_missing=False):
        BinaryStateAlphabet.__init__(self, allow_gaps=allow_gaps, allow_missing=allow_missing)


###############################################################################
## GLOBAL STATE ALPHABETS

DNA_STATE_ALPHABET                =  DnaStateAlphabet()
RNA_STATE_ALPHABET                =  RnaStateAlphabet()
NUCLEOTIDE_STATE_ALPHABET         =  NucleotideStateAlphabet()
BINARY_STATE_ALPHABET             =  BinaryStateAlphabet()
PROTEIN_STATE_ALPHABET            =  ProteinStateAlphabet()
RESTRICTION_SITES_STATE_ALPHABET  =  RestrictionSitesStateAlphabet()
INFINITE_SITES_STATE_ALPHABET     =  InfiniteSitesStateAlphabet()

