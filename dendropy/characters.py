#! /usr/bin/env python

############################################################################
##  characters.py
##
##  Part of the DendroPy phylogenetic computation library.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
############################################################################

"""
This module handles the core definitions of character data types, as well as 
specializations to handle nucleotide, etc. character types.
"""     

from dendropy import base
from dendropy import taxa
from dendropy import utils


class StateAlphabetElement(base.IdTagged):
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
                 elem_id=None, 
                 label=None,
                 symbol=None, 
                 token=None, 
                 multistate=SINGLE_STATE, 
                 member_states=None):
        base.IdTagged.__init__(self, elem_id=elem_id, label=label)
        self.symbol = symbol
        self.token = token
        self.multistate = multistate
        self.member_states = member_states
        
    def __str__(self):
        return str(self.symbol)
        
    def __repr__(self):
        return str([self.elem_id, 
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
        """
        Returns set of id's of all _get_fundamental states to which this state maps.
        """
        return set([state.elem_id for state in self._get_fundamental_states()])
        
    fundamental_ids = property(_get_fundamental_ids)          
        
    def _get_fundamental_symbols(self):
        """
        Returns set of symbols of all _get_fundamental states to which this state maps.
        """
        return set([state.symbol for state in self._get_fundamental_states() if state.symbol]) 
        
    fundamental_symbols = property(_get_fundamental_symbols)
        
    def _get_fundamental_tokens(self):
        """
        Returns set of tokens of all _get_fundamental states to which this state maps.
        """
        return set([state.token for state in self._get_fundamental_states() if state.token])        
        
    fundamental_tokens = property(_get_fundamental_tokens)   
             
class StateAlphabet(base.IdTagged, list):
    """
    A set of states available for a particular character type/format.
    """
    
    def __init__(self, elem_id=None, label=None):
        base.IdTagged.__init__(self, elem_id=elem_id, label=label)  
        list.__init__(self)
        
    def get_state(self, attr_name, value):
        """
        Returns state in self in which attr_name equals value.
        """
        for state in self:
            if getattr(state, attr_name) == value:
                return state
        raise Exception("State with %s value of '%s' not defined" % (attr_name, str(value)))
    
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
        """
        Returns a StateAlphabetElement object corresponding to given symbol.
        """
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
    
    def get_states(self, elem_ids=None, symbols=None, tokens=None):
        """
        Returns list of states with ids/symbols/tokens equal to values
        given in a list of ids/symbols/tokens (exact matches, one-to-one
        correspondence between state and attribute value in list).
        """
        if elem_ids is not None:
            attr_name = 'elem_id'
            values = elem_ids
        elif symbols is not None:
            attr_name = 'symbol'
            values = symbols
        elif tokens is not None:
            attr_name = 'token'
            values = tokens
        else:
            raise Exception("Must specify ids, symbols or tokens")    
        states = []
        for value in values:
            states.append(self.get_state(attr_name=attr_name, value=value))
        return states
    
    def match_state(self, elem_ids=None, symbols=None, tokens=None):
        """
        Returns SINGLE state that has ids/symbols/tokens as member states.
        """
        if elem_ids is not None:
            attr_name = 'fundamental_ids'
            values = elem_ids
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
        """
        Returns dictionary of element id's to state objects.
        """
        map = {}
        for state in self:
            map[state.elem_id] = state
        return map
           
class DnaStateAlphabet(StateAlphabet):

    def __init__(self, elem_id=None, label=None):
        StateAlphabet.__init__(self, elem_id=elem_id, label=label)
        self.append(StateAlphabetElement(symbol="A"))
        self.append(StateAlphabetElement(symbol="C"))     
        self.append(StateAlphabetElement(symbol="G"))
        self.append(StateAlphabetElement(symbol="T"))
        self.append(StateAlphabetElement(symbol="-")) 
        self.append(StateAlphabetElement(symbol="?",
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'G', 'T', '-'])))
        self.append(StateAlphabetElement(symbol="N",
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'G', 'T'])))
        self.append(StateAlphabetElement(symbol="M", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C'])))                                            
        self.append(StateAlphabetElement(symbol="R", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'G'])))
        self.append(StateAlphabetElement(symbol="W",
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'T'])))
        self.append(StateAlphabetElement(symbol="S", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['C', 'G'])))                                            
        self.append(StateAlphabetElement(symbol="Y", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['C', 'T'])))   
        self.append(StateAlphabetElement(symbol="K", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['G', 'T'])))                                            
        self.append(StateAlphabetElement(symbol="V", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'G'])))
        self.append(StateAlphabetElement(symbol="H",
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'T'])))
        self.append(StateAlphabetElement(symbol="D", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'G', 'T'])))                                            
        self.append(StateAlphabetElement(symbol="B", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['C', 'G', 'T'])))
        
class RnaStateAlphabet(StateAlphabet):

    def __init__(self, elem_id=None, label=None):
        StateAlphabet.__init__(self, elem_id=elem_id, label=label)
        self.append(StateAlphabetElement(symbol="A"))
        self.append(StateAlphabetElement(symbol="C"))     
        self.append(StateAlphabetElement(symbol="G"))
        self.append(StateAlphabetElement(symbol="U"))
        self.append(StateAlphabetElement(symbol="-")) 
        self.append(StateAlphabetElement(symbol="?",
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'G', 'U', '-'])))
        self.append(StateAlphabetElement(symbol="N",
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'G', 'U'])))
        self.append(StateAlphabetElement(symbol="M", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C'])))                                            
        self.append(StateAlphabetElement(symbol="R", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'G'])))
        self.append(StateAlphabetElement(symbol="W",
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'U'])))
        self.append(StateAlphabetElement(symbol="S", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['C', 'G'])))                                            
        self.append(StateAlphabetElement(symbol="Y", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['C', 'U'])))   
        self.append(StateAlphabetElement(symbol="K", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['G', 'U'])))                                            
        self.append(StateAlphabetElement(symbol="V", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'G'])))
        self.append(StateAlphabetElement(symbol="H",
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'U'])))
        self.append(StateAlphabetElement(symbol="D", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'G', 'U'])))                                            
        self.append(StateAlphabetElement(symbol="B", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['C', 'G', 'U'])))

class ProteinStateAlphabet(StateAlphabet):

    def __init__(self, elem_id=None, label=None):
        StateAlphabet.__init__(self, elem_id=elem_id, label=label)
        self.append(StateAlphabetElement(symbol="A"))
        self.append(StateAlphabetElement(symbol="C"))
        self.append(StateAlphabetElement(symbol="D"))
        self.append(StateAlphabetElement(symbol="E"))
        self.append(StateAlphabetElement(symbol="F"))
        self.append(StateAlphabetElement(symbol="G"))
        self.append(StateAlphabetElement(symbol="H"))
        self.append(StateAlphabetElement(symbol="I"))
        self.append(StateAlphabetElement(symbol="K"))
        self.append(StateAlphabetElement(symbol="L"))
        self.append(StateAlphabetElement(symbol="M"))
        self.append(StateAlphabetElement(symbol="N"))
        self.append(StateAlphabetElement(symbol="P"))
        self.append(StateAlphabetElement(symbol="Q"))
        self.append(StateAlphabetElement(symbol="R"))
        self.append(StateAlphabetElement(symbol="S"))
        self.append(StateAlphabetElement(symbol="T"))
        self.append(StateAlphabetElement(symbol="U"))
        self.append(StateAlphabetElement(symbol="V"))
        self.append(StateAlphabetElement(symbol="W"))
        self.append(StateAlphabetElement(symbol="Y"))
        self.append(StateAlphabetElement(symbol="-"))        
        self.append(StateAlphabetElement(symbol="B", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['D', 'N'])))     
        self.append(StateAlphabetElement(symbol="Z", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['E', 'Q'])))
        self.append(StateAlphabetElement(symbol="X", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'D', 'E',
                                           'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
                                           'T', 'U', 'V', 'W', 'Y'])))   
        self.append(StateAlphabetElement(symbol="?", 
                                           multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                           member_states=self.get_states(symbols=['A', 'C', 'D', 'E',
                                           'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
                                           'T', 'U', 'V', 'W', 'Y', '-'])))            

class BinaryStateAlphabet(StateAlphabet):

    def __init__(self, elem_id=None, label=None, allow_gaps=True, allow_missing=True):
        StateAlphabet.__init__(self, elem_id=elem_id, label=label)
        self.append(StateAlphabetElement(symbol="0"))
        self.append(StateAlphabetElement(symbol="1"))
        if allow_gaps:
            self.append(StateAlphabetElement(symbol="-")) 
            if allow_missing:
                self.append(StateAlphabetElement(symbol="?", 
                                                   multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                                   member_states=self.get_states(symbols=['0', '1', '-'])))
        elif allow_missing:
            self.append(StateAlphabetElement(symbol="?", 
                                               multistate=StateAlphabetElement.AMBIGUOUS_STATE,
                                               member_states=self.get_states(symbols=['0', '1'])))
                        
class RestrictionSitesStateAlphabet(BinaryStateAlphabet):

    def __init__(self, elem_id=None, label=None):
        BinaryStateAlphabet.__init__(self, elem_id=elem_id, label=label)    
        
class InfiniteSitesStateAlphabet(BinaryStateAlphabet):

    def __init__(self, elem_id=None, label=None):
        BinaryStateAlphabet.__init__(self, elem_id=elem_id, label=label)

### GLOBAL STATE ALPHABETS ###                                       

DNA_STATE_ALPHABET = DnaStateAlphabet()                                           
RNA_STATE_ALPHABET = RnaStateAlphabet() 
PROTEIN_STATE_ALPHABET = ProteinStateAlphabet()                                            
RESTRICTION_SITES_STATE_ALPHABET = RestrictionSitesStateAlphabet()         
INFINITE_SITES_STATE_ALPHABET = InfiniteSitesStateAlphabet()
                
class ColumnType(base.IdTagged):
    """                                                                                                                                                                                                                                                                                                                                                                           
    A character format or type of a particular column: i.e., maps
    a particular set of character state definitions to a column in a character char_block.
    """
  
    def __init__(self, elem_id=None,label=None, state_alphabet=None):
        base.IdTagged.__init__(self, elem_id=elem_id, label=label)
        self.__state_alphabet = None
        self.id_state_map = None
        self.state_alphabet = state_alphabet
        
    def _set_state_alphabet(self, value):
        self.__state_alphabet = value
        if self.__state_alphabet is not None:
            self.id_state_map = self.__state_alphabet.id_state_map()
        else:
            self.id_state_map = None
            
    def _get_state_alphabet(self):
        return self.__state_alphabet
        
    state_alphabet = property(_get_state_alphabet, _set_state_alphabet)
      
class CharacterDataCell(base.Annotated):
    """                                                                                                                                                                                                                                                                                                                                                                           
    A container for the state / state value for a particular cell in a char_block.
    """
  
    def __init__(self, value=None, column_type=None):
        base.Annotated.__init__(self)
        self.value = value
        self.column_type = column_type 
        
    def __str__(self):
        return str(self.value)
        
    def __eq__(self, other):
        if isinstance(other, CharacterDataCell):
            return self.value == other.value
#         elif isinstance(other, self.value):
#             
        else:
            return NotImplemented
            
    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result
                
class CharacterDataVector(list, taxa.TaxonLinked):
    """
    A list of character data values.
    """

    def __init__(self, elem_id=None, label=None, taxon=None):
        list.__init__(self)
        taxa.TaxonLinked.__init__(self, elem_id=elem_id, label=label, taxon=taxon)
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
        if self:
            return [str(cell.value) for cell in self]
        else:
            return []
            
    def values_as_string(self, sep=""):
        if self:
            return sep.join([str(cell.value) for cell in self])
        else:
            return ""        
            
    def __str__(self):
        return str(self.values_as_string_list())

class CharacterDataMatrix(dict, base.Annotated):
    """
    An annotable dictionary with Taxon objects as keys and 
    CharacterDataVectors objects as values.
    """

    def __init__(self):
        #utils.OrderedCaselessDict.__init__(self) 
        dict.__init__(self)
        base.Annotated.__init__(self)
        
    def extend_characters(self, other_matrix):
        """
        Extends this char_block by adding characters from sequences of taxa
        in given char_block to sequences of taxa with correspond labels in
        this one. Taxa in the second char_block that do not exist in the
        current one are ignored.
        """
        label_taxon_map = dict([(taxon.label, taxon) for taxon in other_matrix])
        for taxon in self:
            if taxon.label in label_taxon_map:
                self[taxon].extend(other_matrix[label_taxon_map[taxon.label]])
                
    def extend(self, 
        other_matrix, 
        overwrite_existing=False, 
        append_existing=False):
        """
        Extends this char_block by adding taxa and characters from the given
        char_block to this one.  If `overwrite_existing` is True and a taxon
        in the other char_block is already present in the current one, then
        the sequence associated with the taxon in the second char_block
        replaces the sequence in the current one. If `append_existing`
        is True and a taxon in the other char_block is already present in
        the current one, then the squence associated with the taxon in
        the second char_block will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True,  and a taxon in the other char_block is already present in
        the current one, then the sequence is ignored.
        Note that the containing CharactersBlock taxa_block has to be normalized 
        after this operation.
        """
        if overwrite_existing and append_existing:
            raise Exception("Can only specify to overwrite or append, not both")
        label_taxon_map = dict([(taxon.label, taxon) for taxon in self])
        for other_taxon in other_matrix:
            if other_taxon.label in label_taxon_map:
                this_taxon = label_taxon_map[other_taxon.label]
                if overwrite_existing:
                    self[this_taxon] = other_matrix[other_taxon]
                elif append_existing:
                    self[this_taxon].extend(other_matrix[other_taxon])
            else:
                self[other_taxon] = other_matrix[other_taxon]
                                                                                                                   
class CharactersBlock(taxa.TaxaLinked):
    """
    Character data container/manager manager.
    """

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        taxa.TaxaLinked.__init__(self, *args, **kwargs)
        self.matrix = CharacterDataMatrix()
        self.column_types = []
        self.markup_as_sequences = True
        
    def extend_characters(self, other_char_block):
        """
        Extends this char_block by adding characters from sequences of taxa
        in given char_block to sequences of taxa with correspond labels in
        this one. Taxa in the second char_block that do not exist in the
        current one are ignored.
        """
        self.matrix.extend_characters(other_char_block.matrix)
                
    def extend(self, 
        other_char_block, 
        overwrite_existing=False, 
        append_existing=False):
        """
        Extends this char_block by adding taxa and characters from the given
        char_block to this one.  If `overwrite_existing` is True and a taxon
        in the other char_block is already present in the current one, then
        the sequence associated with the taxon in the second char_block
        replaces the sequence in the current one. If `append_existing`
        is True and a taxon in the other char_block is already present in
        the current one, then the squence associated with the taxon in
        the second char_block will be added to the sequence in the current
        one. If both are True, then an exception is raised. If neither
        are True,  and a taxon in the other char_block is already present in
        the current one, then the sequence is ignored.
        """
        self.matrix.extend(other_char_block.matrix, 
            overwrite_existing=overwrite_existing, 
            append_existing=append_existing)
        for taxon in self:
            if taxon not in self.taxa_block:
                self.taxa_block.append(taxon)
        
    def vectors(self):
        """
        Returns list of vectors.        
        """
        if self.taxa_block is not None and self.matrix is not None:
            if len(self.matrix) > 0:
                return [self.matrix[t] for t in self.taxa_block]
            else:
                return []
        else:
            return None
            
    # following allows a CharactersBlock object to simulate a dictionary            
    # by `passing-through` calls to the underlying matrix
        
        
    def __len__(self):
        """
        Dictionary interface implementation for direct access to matrix.
        """    
        return len(self.matrix)
        
    def __getitem__(self, key):
        """
        Dictionary interface implementation for direct access to matrix.
        """
        return self.matrix[key]
        
    def __setitem__(self, key, value):
        """
        Dictionary interface implementation for direct access to matrix.
        """
        self.matrix[key] = value
                
    def __contains__(self, key):
        """
        Dictionary interface implementation for direct access to matrix.
        """
        return key in self.matrix        
                
    def iterkeys(self):
        """
        Dictionary interface implementation for direct access to matrix.
        """
        for key in self.matrix:
            yield(key)
    
    def itervalues(self):
        """
        Dictionary interface implementation for direct access to matrix.
        """
        for value in self.matrix.values():
            yield(value)
    
    def iteritems(self):
        """
        Returns an iterator over matrix's values.
        """
        for key, value in self.matrix.iteritems():
            yield (key, value)

    def items(self):
        """
        Returns matrix key, value pairs in key-order.
        """
        return [(key, self.matrix[key]) for key in self.matrix.iterkeys()]

    def values(self):
        """
        Returns list of matrix key, value pairs.
        """
        return [v for v in self.matrix.itervalues()]
    
    def __iter__(self):
        """
        Returns an iterator over matrix's ordered keys.
        """
        return self.matrix.iterkeys()
    
    def __delitem__(self, key):
        """
        Remove item from matrix with specified key.
        """
        super(dict, self.matrix).__delitem__(key)                

    def __contains__(self, key):
        """
        Returns true if matrix has key, regardless of case.
        """
        return dict(self.matrix).__contains__(key)

    def pop(self, key, alt_val=None):
        """
        a.pop(k[, x]):  a[k] if k in a, else x (and remove k)
        """
        if key in self.matrix:
            val = self.matrix[key]
            self.matrix.__delitem__(key)
            return val
        else:
            return alt_val
        
    def popitem(self):
        """
        a.popitem()  remove and last (key, value) pair
        """
        key = self.matrix[-1]
        item = (key, self.matrix[key])
        self.matrix.__delitem__(key)
        return item

    def index(self, key):
        """
        Return the index of (caseless) key.
        Raise KeyError if not found.
        """
        count = 0
        for k in self.matrix.keys():
            if k == key:
                return count
            count = count + 1
        raise KeyError(key)

    def keys(self):
        """
        Returns a copy of the ordered list of matrix keys.
        """
        return list(self.matrix.keys())

    def clear(self):
        """
        Deletes all items from the matrix dictionary.
        """
        self.matrix.clear()

    def has_key(self, key):
        """
        Returns true if matrix has key, regardless of case.
        """
        return key in self.matrix

    def get(self, key, def_val=None):
        """
        Gets an item from matrix by its key, returning default if key not present.
        """
        return super(dict, self.matrix).get(key, def_val)

    def setdefault(self, key, def_val=None):
        """
        Sets the default value to return if key not present.
        """
        return super(dict, self.matrix).setdefault(key, def_val)
      
    def id_column_map(self):
        """
        Returns dictionary of element id to corresponding
        character definition.
        """
        map = {}
        for char in self.column_types:
            map[char.elem_id] = char
        return map
        
class ContinuousCharactersBlock(CharactersBlock):
    """
    Character data container/manager manager.
    """

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        CharactersBlock.__init__(self, *args, **kwargs)

class DiscreteCharactersBlock(CharactersBlock):
    """
    Character data container/manager manager.
    """

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        CharactersBlock.__init__(self, *args, **kwargs)
        self.state_alphabets = []
        self.default_state_alphabet = None
        self.__default_symbol_state_map = None
        
    def _get_default_symbol_state_map(self):
        if self.__default_symbol_state_map is None and self.default_state_alphabet is not None:
            self.__default_symbol_state_map = self.default_state_alphabet.symbol_state_map()
        return self.__default_symbol_state_map
    
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
                    
class StandardCharactersBlock(DiscreteCharactersBlock):
    """
    `standard` data.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        DiscreteCharactersBlock.__init__(self, *args, **kwargs)
                         
class DnaCharactersBlock(DiscreteCharactersBlock):
    """
    DNA nucleotide data.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        DiscreteCharactersBlock.__init__(self, *args, **kwargs)
        self.default_state_alphabet = DNA_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)         

class RnaCharactersBlock(DiscreteCharactersBlock):
    """
    RNA nucleotide data.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        DiscreteCharactersBlock.__init__(self, *args, **kwargs)
        self.default_state_alphabet = RNA_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)      
        
class ProteinCharactersBlock(DiscreteCharactersBlock):
    """
    Protein / amino acid data.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """        
        DiscreteCharactersBlock.__init__(self, *args, **kwargs)
        self.default_state_alphabet = PROTEIN_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)               
        
class RestrictionSitesCharactersBlock(DiscreteCharactersBlock):
    """
    Restriction sites data.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        DiscreteCharactersBlock.__init__(self, *args, **kwargs)
        self.default_state_alphabet = RESTRICTION_SITES_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)         
        
class InfiniteSitesCharactersBlock(DiscreteCharactersBlock):
    """
    Infinite sites data.
    """
    
    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        DiscreteCharactersBlock.__init__(self, *args, **kwargs)
        self.default_state_alphabet = INFINITE_SITES_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)         

if __name__ == "__main__":
    dna = DnaStateAlphabets()
    for s in dna:
        print repr(s)
    print
    print
    input = "\n"
    while input:
        input = raw_input("Enter symbol(s): ")
        if input:
            input = input.upper()
            print repr(dna.match_state(symbols=input))