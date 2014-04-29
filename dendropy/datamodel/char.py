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
Character and character-sequence data structures.
"""

import warnings
import copy
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.utility import error
from dendropy.utility import container
from dendropy.datamodel.statealphabet import DNA_STATE_ALPHABET
from dendropy.datamodel.statealphabet import RNA_STATE_ALPHABET
from dendropy.datamodel.statealphabet import NUCLEOTIDE_STATE_ALPHABET
from dendropy.datamodel.statealphabet import PROTEIN_STATE_ALPHABET
from dendropy.datamodel.statealphabet import RESTRICTION_SITES_STATE_ALPHABET
from dendropy.datamodel.statealphabet import INFINITE_SITES_STATE_ALPHABET
from dendropy.datamodel import base
from dendropy.datamodel import taxon
from dendropy import dataio

###############################################################################
## Continuous Characters

class ContinuousCharElement(
        base.DataObject,
        base.Annotable):
    def __init__(self, value, column_def,  **kwargs):
        base.DataObject.__init__(self,
                label=kwargs.pop("label", None))
        self.value = value
        self.column_def = column_def
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

###############################################################################
## Data Containers

class CharacterType(
        base.DataObject,
        base.Annotable):
    """
    A character format or type of a particular column: i.e., maps
    a particular set of character state definitions to a column in a character matrix.
    """

    def __init__(self, **kwargs):
        base.DataObject.__init__(self,
                label=kwargs.pop("label", None))
        self._state_alphabet = None
        self.state_alphabet = kwargs.pop("state_alphabet", None)
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def _set_state_alphabet(self, value):
        self._state_alphabet = value
    def _get_state_alphabet(self):
        return self._state_alphabet
    state_alphabet = property(_get_state_alphabet, _set_state_alphabet)

    def __copy__(self, memo=None):
        raise TypeError("Cannot directly copy {}".format(self.__class__.__name__))

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot directly copy {}".format(self.__class__.__name__))

    def __deepcopy__(self, memo=None):
        return base.Annotable.__deepcopy__(self, memo=memo)

class CharacterDataVector(object):
    """A list of character data values for a taxon -- a row of a Character Matrix.

    The CharacterDataVector typically contains elements that are instances of
    CharacterDataCell
    """

    def __init__(self, *args, **kwargs):
        self._values = list(*args)
        self.taxon = kwargs.pop("taxon", None)
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def __deepcopy__(self, memo):
        other = CharacterDataVector()
        other.taxon = copy.deepcopy(self.taxon, memo)
        memo[id(self.taxon)] = other.taxon
        for v in self._values:
            v2 = copy.deepcopy(v, memo)
            o._values.append(v2)
            memo[id(v)] = v2
        return o

    def values(self):
        return list(self._values)

    def symbols_as_list(self):
        return [str(v) for v in self._values]

    def symbols_as_string(self, sep=""):
        return sep.join(self.symbols_as_list())

    def __str__(self):
        return str(self.symbols_as_string())

    def insert(self, index, element):
        return self._values.insert(index, element)

    def append(self, element):
        return self._values.append(element)

    def extend(self, other):
        self._values.extend(other)

    def __iadd__(self, other):
        raise NotImplementedError

    def __add__(self, other):
        raise NotImplementedError

    def __contains__(self, element):
        raise NotImplementedError

    def __delitem__(self, element):
        raise NotImplementedError

    def __iter__(self):
        return iter(self._values)

    def __reversed__(self):
        return reversed(self._values)

    def __len__(self):
        return len(self._values)

    def __getitem__(self, index):
        if isinstance(index, slice):
            raise NotImplementedError
        else:
            return self._values[index]

    def __setitem__(self, index, value):
        if isinstance(index, slice):
            raise NotImplementedError
        else:
            self._values[index] = value

    def clear(self):
        # list.clear() only with 3.4 or so ...
        self._values = []

    def index(self, element):
        return self._values.index(element)

    def pop(self, index=-1):
        return self._values.pop(index)

    def remove(self, element):
        self._values.remove(element)

    def reverse(self):
        self._values.reverse()

    def sort(self, key=None, reverse=False):
        self._values.sort(key=key, reverse=reverse)

class CharacterDataMap(dict,
        base.DataObject,
        base.Annotable):
    """
    An annotable dictionary with Taxon objects as keys and
    CharacterDataVectors objects as values.
    """

    def __init__(self, label=None):
        dict.__init__(self)
        base.DataObject.__init__(self, label=label)

    def __copy__(self, memo=None):
        raise TypeError("Cannot directly copy {}".format(self.__class__.__name__))

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot directly copy {}".format(self.__class__.__name__))

    def __deepcopy__(self, memo=None):
        return base.Annotable.__deepcopy__(self, memo=memo)

    def _get_vector_size(self):
        """
        Returns number of characters in *first* sequence.
        """
        if len(self):
            # yuck, but len(self.values())
            # means we have to create and populatre a list ...
            return len(self[next(iter(self))])
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

class CharacterSubset(
        base.DataObject,
        base.Annotable,
        ):
    """
    Tracks definition of a subset of characters.
    """

    def __init__(self, label=None, character_indices=None):
        """
        Keyword arguments:

            - `label`: name of this subset
            - `character_indices`: list of 0-based (integer) indices
               of column positions that constitute this subset.

        """
        base.DataObject.__init__(self, label=label)
        if character_indices is None:
            self.character_indices = set()
        else:
            self.character_indices = set(character_indices)

    def __len__(self):
        return len(self.character_indices)

    def __iter__(self):
        return iter(self.character_indices)

    def __deepcopy__(self, memo):
        return base.Annotable.__deepcopy__(self, memo=memo)

###############################################################################
## Base Character Matrix

class CharacterMatrix(
        taxon.TaxonNamespaceAssociated,
        base.Annotable,
        base.Readable,
        base.Writeable,
        base.DataObject):
    "Character data container/manager manager."

    data_type = None

    def _parse_from_stream(cls,
            stream,
            schema,
            matrix_offset=0,
            **kwargs):
        taxon_namespace = taxon.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        if taxon_namespace is None:
            taxon_namespace = taxon.TaxonNamespace()
        tns_factory = lambda label: taxon_namespace
        label = kwargs.pop("label", None)
        matrix_type = get_char_matrix_type(data_type=cls.data_type)
        reader = dataio.get_reader(schema, **kwargs)
        char_matrices = reader.read_char_matrices(
                stream=stream,
                taxon_namespace_factory=tns_factory,
                char_matrix_factory=matrix_type,
                global_annotations_target=None)
        if len(char_matrices) == 0:
            raise ValueError("No character data in data source")
        # if matrix_offset >= len(char_matrices):
        #     raise IndexError(
        #             "Character matrix of offset {} specified, but data source only has"
        #             " {} matrices defined (maximum offset = {})".format(
        #             offset, len(char_matrices), len(char_matrices)-1))
        char_matrix = char_matrices[matrix_offset]
        # if not isinstance(d.char_matrices[index], cls):
        #     raise ValueError("Character data found was of type '%s' (object is of type '%s')" %
        #             (d.char_matrices[index].__class__.__name__, cls.__name__))
        # return d.char_matrices[index]
        return char_matrix
    _parse_from_stream = classmethod(_parse_from_stream)

    def concatenate(cls, char_matrices):
        """
        Creates and returns a single character matrix from multiple
        CharacterMatrix objects specified as a list, 'char_matrices'.
        All the CharacterMatrix objects in the list must be of the
        same type, and share the same TaxonNamespace reference. All taxa
        must be present in all alignments, all all alignments must
        be of the same length. Component parts will be recorded as
        character subsets.
        """
        taxon_namespace = char_matrices[0].taxon_namespace
        nseqs = len(char_matrices[0])
        concatenated_chars = cls(taxon_namespace=taxon_namespace)
        pos_start = 0
        for cidx, cm in enumerate(char_matrices):
            if cm.taxon_namespace is not taxon_namespace:
                raise ValueError("Different `taxon_namespace` references in matrices to be merged")
            if len(cm) != len(taxon_namespace):
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
        taxon_namespace = taxon.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        if taxon_namespace is None:
            taxon_namespace = TaxonNamespace()
        kwargs["taxon_namespace"] = taxon_namespace
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
        base.DataObject.__init__(self, label=kwargs.pop("label", None))
        taxon.TaxonNamespaceAssociated.__init__(self,
                taxon_namespace=taxon.process_kwargs_dict_for_taxon_namespace(kwargs, None))
        self.taxon_seq_map = CharacterDataMap()
        self.character_types = []
        self.character_subsets = container.OrderedCaselessDict()
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
        "TODO: may need to check that we are not overwriting oid"
        if len(args) > 1:
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        elif len(args) == 1:
            if isinstance(args[0],  self.__class__):
                ca = copy.deepcopy(args[0])
                for k, v in ca.__dict__.iteritems():
                    if k not in ["_annotations"]:
                        self.__dict__[k] = v
                self.annotations = ca.annotations
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
        warnings.warn("Repopulating a CharacterMatrix is now deprecated. Instantiate a new instance from the source instead.",
                DeprecationWarning)
        m = self.__class__._parse_from_stream(stream=stream,
                schema=schema,
                **kwargs)
        return self.clone_from(m)

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
        self.update_taxon_namespace()

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
        self.update_taxon_namespace()

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
        clone = self.__class__(self)
        # clone.clone_from(self)
        for vec in clone.taxon_seq_map.values():
            for cell_idx in range(len(vec)-1, -1, -1):
                if cell_idx not in indices:
                    del(vec[cell_idx])
        return clone

    def reindex_subcomponent_taxa(self):
        """
        Synchronizes `Taxon` objects of map to `taxon_namespace` of self.
        """
        ti_mutable = self.taxon_namespace._is_mutable
        self.taxon_namespace._is_mutable = True
        new_map = CharacterDataMap()
        for taxon, seq in self.taxon_seq_map.items():
            taxon = self.taxon_namespace.require_taxon(label=taxon.label)
            new_map[taxon] = seq
        self.taxon_namespace._is_mutable = ti_mutable
        self.taxon_seq_map = new_map

    def update_taxon_namespace(self):
        """
        Updates local taxa block by adding taxa not already managed.
        Mainly for use after map extension
        """
        assert self.taxon_namespace is not None
        for taxon in self.taxon_seq_map:
            if taxon not in self.taxon_namespace:
                self.taxon_namespace.add(taxon)

    def vectors(self):
        "Returns list of vectors.        "
        if self.taxon_namespace is not None and self.taxon_seq_map is not None:
            if len(self.taxon_seq_map) > 0:
                return [self.taxon_seq_map[t] for t in self.taxon_namespace]
            return []
        return None

    def prune_taxa(self, taxa, update_taxon_namespace=False):
        """
        Removes given taxa from matrix. If `preserve_taxon_namespace` is
        `True`, then the taxa are removed from the associated TaxonNamespace
        object as well. Otherwise this is not modified (default).
        """
        for taxon in taxa:
            if taxon in self.taxon_seq_map:
                del self.taxon_seq_map[taxon]
                if update_taxon_namespace and taxon in self.taxon_namespace:
                    self.taxon_namespace.remove(taxon)

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
            if key >= 0 and key < len(self.taxon_namespace):
                key = self.taxon_namespace[key]
            else:
                raise KeyError(key)
        elif isinstance(key, str):
            label = key
            key = None
            for t in self.taxon_namespace:
                if t.label == label:
                    key = t
                    break
            if key is None:
                raise KeyError(label)
        return self.taxon_seq_map[key]

    def __setitem__(self, key, value):
        "Dictionary interface implementation for direct access to character map."
        if isinstance(key, int):
            if key >= 0 and key < len(self.taxon_namespace):
                key = self.taxon_namespace[key]
            else:
                raise KeyError(key)
        if key not in self.taxon_namespace:
            self.taxon_namespace.add(key)
        self.taxon_seq_map[key] = value

    def iterkeys(self):
        "Dictionary interface implementation for direct access to character map."
        for t in self.taxon_namespace:
            if t in self.taxon_seq_map:
                yield t

    def itervalues(self):
        "Dictionary interface implementation for direct access to character map."
        for t in self.taxon_namespace:
            if t in self.taxon_seq_map:
                yield self.taxon_seq_map[t]

    def iteritems(self):
        "Returns an iterator over character map's values."
        for t in self.taxon_namespace:
            if t in self.taxon_seq_map:
                yield t, self.taxon_seq_map[t]

    def items(self):
        "Returns character map key, value pairs in key-order."
        return [(t, self.taxon_seq_map[t]) for t in self.taxon_namespace if t in self.taxon_seq_map]

    def values(self):
        "Returns list of values."
        return [self.taxon_seq_map[t] for t in self.taxon_namespace if t in self.taxon_seq_map]

    def __iter__(self):
        "Returns an iterator over character map's ordered keys."
        for t in self.taxon_namespace:
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
                if self.taxon_namespace is not None:
                    tlead = "\n%s[Taxon Set]\n" % (" " * (indent+4))
                    output_strio.write(tlead)
                    self.taxon_namespace.description(depth=depth-1, indent=indent+8, itemize="", output=output_strio)
                tlead = "\n%s[Characters]\n" % (" " * (indent+4))
                output_strio.write(tlead)
                indent += 8
                maxlabel = max([len(str(t.label)) for t in self.taxon_namespace])
                for i, t in enumerate(self.taxon_namespace):
                    output_strio.write('%s%s%s : %s characters\n' \
                        % (" " * indent,
                           "[%d] " % i,
                           str(t.label),
                           len(self.taxon_seq_map[t])))

        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    def new_character_data_vector(self, **kwargs):
        return CharacterDataVector(**kwargs)

###############################################################################
## Specialized Matrices

class ContinuousCharacterMatrix(CharacterMatrix):
    "Character data container/manager manager."

    data_type = "continuous"

    def __init__(self, *args, **kwargs):
        "See CharacterMatrix.__init__ documentation"
        CharacterMatrix.__init__(self, *args, **kwargs)

class DiscreteCharacterMatrix(CharacterMatrix):
    """Character data container/manager manager.

    That adds the attributes self.state_alphabets (a list of alphabets)
    and self.default_state_alphabet
    """

    data_type = "fixed"

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

    def remap_to_state_alphabet_by_symbol(self,
            state_alphabet,
            purge_other_state_alphabets=True):
        """
        All entities with any reference to a state alphabet will be have the
        reference reassigned to state alphabet ``sa``, and all entities with
        any reference to a state alphabet element will be have the reference
        reassigned to any state alphabet element in ``sa`` that has the same
        symbol. Raises KeyError if no matching symbol can be found.
        """
        symbol_state_map = state_alphabet.symbol_state_map()
        for vi, vec in enumerate(self.taxon_seq_map.values()):
            for ci, cell in enumerate(vec):
                cell.value = symbol_state_map[cell.value.symbol]
        for ct in self.character_types:
            ct.state_alphabet = state_alphabet
        if purge_other_state_alphabets:
            self.state_alphabets = [state_alphabet]
            self.default_state_alphabet = state_alphabet

    def remap_to_default_state_alphabet_by_symbol(self,
            purge_other_state_alphabets=True):
        """
        All entities with any reference to a state alphabet will be have the
        reference reassigned to the default state alphabet, and all entities
        with any reference to a state alphabet element will be have the
        reference reassigned to any state alphabet element in the default
        state alphabet that has the same symbol. Raises ValueError if no
        matching symbol can be found.
        """
        self.remap_to_state_alphabet_by_symbol(
                state_alphabet=self.default_state_alphabet,
                purge_other_state_alphabets=purge_other_state_alphabets)

class StandardCharacterMatrix(DiscreteCharacterMatrix):
    "`standard` data."

    data_type = "standard"

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        DiscreteCharacterMatrix.__init__(self, **kwargs)
        if len(args) > 0:
            self.clone_from(*args)

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

    data_type = "fixed"

    def __init__(self, *args, **kwargs):
        """See CharacterMatrix.__init__ documentation for kwargs.

        Unnamed args are passed to clone_from.
        """
        DiscreteCharacterMatrix.__init__(self, **kwargs)
        if len(args) > 0:
            self.clone_from(*args)

    def __deepcopy__(self, memo):
        memo[id(self.default_state_alphabet)] = self.default_state_alphabet
        memo[id(self.state_alphabets)] = list(self.state_alphabets)
        # memo[id(self._default_symbol_state_map)] = self._default_symbol_state_map
        return DiscreteCharacterMatrix.__deepcopy__(self, memo)

class DnaCharacterMatrix(FixedAlphabetCharacterMatrix):
    "DNA nucleotide data."

    data_type = "dna"

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

    data_type = "rna"

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

    data_type = "nucleotide"

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxon_namespace`."
        FixedAlphabetCharacterMatrix.__init__(self, **kwargs)
        self.default_state_alphabet = NUCLEOTIDE_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
        if len(args) > 0:
            self.clone_from(*args)

class ProteinCharacterMatrix(FixedAlphabetCharacterMatrix):
    "Protein / amino acid data."

    data_type = "protein"

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `oid`, `label` and `taxon_namespace`.
        """
        FixedAlphabetCharacterMatrix.__init__(self, **kwargs)
        self.default_state_alphabet = PROTEIN_STATE_ALPHABET
        self.state_alphabets.append(self.default_state_alphabet)
        if len(args) > 0:
            self.clone_from(*args)

class RestrictionSitesCharacterMatrix(FixedAlphabetCharacterMatrix):
    "Restriction sites data."

    data_type = "restriction"

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

    data_type = "infinite"

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
## Main Character Matrix Factory Function

data_type_matrix_map = {
    'continuous' : ContinuousCharacterMatrix,
    'dna' : DnaCharacterMatrix,
    'rna' : RnaCharacterMatrix,
    'protein' : ProteinCharacterMatrix,
    'standard' : StandardCharacterMatrix,
    'restriction' : RestrictionSitesCharacterMatrix,
    'infinite' : InfiniteSitesCharacterMatrix,
}

def get_char_matrix_type(data_type):
    if data_type is None:
        raise TypeError("'data_type' must be specified")
    matrix_type = data_type_matrix_map.get(data_type, None)
    if matrix_type is None:
        raise KeyError("Unrecognized data type specification: '{}'".format(data_type,
            sorted(data_type_matrix_map.keys())))
    return matrix_type

def new_char_matrix(data_type, **kwargs):
    matrix_type = get_char_matrix(data_type=data_type)
    m = matrix_type(**kwargs)
    return m

###############################################################################
## Wrappers and Convenience Functions

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
