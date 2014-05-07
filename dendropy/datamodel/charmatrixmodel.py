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
from dendropy.datamodel.charstatemodel import DNA_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import RNA_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import NUCLEOTIDE_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import PROTEIN_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import RESTRICTION_SITES_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import INFINITE_SITES_STATE_ALPHABET
from dendropy.datamodel import basemodel
from dendropy.datamodel import taxonmodel
from dendropy import dataio

###############################################################################
## ContinuousCharElement

class ContinuousCharElement(
        basemodel.DataObject,
        basemodel.Annotable):
    def __init__(self, value, column_def=None, label=None):
        basemodel.DataObject.__init__(self,
                label=label)
        self.value = value
        self.column_def = column_def

###############################################################################
## CharacterType

class CharacterType(
        basemodel.DataObject,
        basemodel.Annotable):
    """
    A character format or type of a particular column: i.e., maps
    a particular set of character state definitions to a column in a character matrix.
    """

    def __init__(self,
            label=None,
            state_alphabet=None):
        basemodel.DataObject.__init__(self, label=label)
        self._state_alphabet = None
        self.state_alphabet = state_alphabet

    def _get_state_alphabet(self):
        """
        The :class:`StateAlphabet` representing the state alphabet for this
        column: i.e., the collection of symbols and the state identities to
        which they map.
        """
        return self._state_alphabet
    def _set_state_alphabet(self, value):
        self._state_alphabet = value
    state_alphabet = property(_get_state_alphabet, _set_state_alphabet)

    def __copy__(self, memo=None):
        raise TypeError("Cannot directly copy {}".format(self.__class__.__name__))

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot directly copy {}".format(self.__class__.__name__))

    def __deepcopy__(self, memo=None):
        return basemodel.Annotable.__deepcopy__(self, memo=memo)

###############################################################################
## CharacterSequence

class CharacterSequence(
        basemodel.Annotable,
        list
        ):
    """
    A sequence of character states or values for a particular taxon or entry in
    a data matrix. Extends list by supporting metadata annotation, and richer
    suite of representing data in self (e.g., as list or string of symbols
    instead of just the raw data values).
    """

    def __init__(self, values=None):
        """
        Parameters
        ----------
        values : iterable of states
            A set of values for this sequence.
        """
        if values:
            super(CharacterSequence, self).__init__(values)
        else:
            super(CharacterSequence, self).__init__()

    def values(self):
        """
        Returns list of values of this vector.

        Returns
        -------
        v : list
            List of values making up this vector.
        """
        return list(self)

    def symbols_as_list(self):
        """
        Returns list of string representation of values of this vector.

        Returns
        -------
        v : list
            List of string representation of values making up this vector.
        """
        return [str(v) for v in self]

    def symbols_as_string(self, sep=""):
        """
        Returns values of this vector as a single string, with individual value
        elements separated by `sep`.
        Returns
        -------
        s : string
            String representation of values making up this vector.
        """
        return sep.join(self.symbols_as_list())

    def __str__(self):
        return self.symbols_as_string()

###############################################################################
## Subset of Character (Columns)

class CharacterSubset(
        basemodel.DataObject,
        basemodel.Annotable,
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
        basemodel.DataObject.__init__(self, label=label)
        if character_indices is None:
            self.character_indices = set()
        else:
            self.character_indices = set(character_indices)

    def __len__(self):
        return len(self.character_indices)

    def __iter__(self):
        return iter(self.character_indices)

    def __deepcopy__(self, memo):
        return basemodel.Annotable.__deepcopy__(self, memo=memo)

###############################################################################
## CharacterMatrix

class CharacterMatrix(
        taxonmodel.TaxonNamespaceAssociated,
        basemodel.Annotable,
        basemodel.Readable,
        basemodel.Writeable,
        basemodel.DataObject):
    """
    A data structure that manages assocation of operational taxononomic unit
    concepts to sequences of character state identities or values.
    """

    ###########################################################################
    ### Class Variables

    data_type = None
    character_sequence_type = CharacterSequence

    ###########################################################################
    ### Factory (Class) Methods

    def _parse_from_stream(cls,
            stream,
            schema,
            matrix_offset=0,
            **kwargs):
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        if taxon_namespace is None:
            taxon_namespace = taxonmodel.TaxonNamespace()
        def tns_factory(label):
            if label is not None and taxon_namespace.label is None:
                taxon_namespace.label = label
            return taxon_namespace
        label = kwargs.pop("label", None)
        kwargs["data_type"] = cls.data_type
        reader = dataio.get_reader(schema, **kwargs)
        char_matrices = reader.read_char_matrices(
                stream=stream,
                taxon_namespace_factory=tns_factory,
                char_matrix_factory=new_char_matrix,
                global_annotations_target=None)
        if len(char_matrices) == 0:
            raise ValueError("No character data in data source")
        char_matrix = char_matrices[matrix_offset]
        if char_matrix.data_type != cls.data_type:
            raise ValueError(
                "Data source (at offset {}) is of type '{}', "
                "but current CharacterMatrix is of type '{}'.".format(char_matrix.data_type,
                    cls.data_type))
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
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
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

    ###########################################################################
    ### Lifecycle and Identity

    def __init__(self, *args, **kwargs):
        basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
        taxonmodel.TaxonNamespaceAssociated.__init__(self,
                taxon_namespace=taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None))
        self._taxon_sequence_map = {}
        self.character_types = []
        self.character_subsets = container.OrderedCaselessDict()
        self.markup_as_sequences = True

        if len(args) > 0:
            raise NotImplementedError

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

        # if len(args) > 1:
        #     raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        # if len(args) == 1:
        #     if ("stream" in kwargs and kwargs["stream"] is not None) \
        #             or ("schema" in kwargs and kwargs["schema"] is not None):
        #         raise error.MultipleInitializationSourceError(class_name=self.__class__.__name__, arg=args[0])
        #     if isinstance(args[0], self.__class__):
        #         self.clone_from(args[0])
        #     else:
        #         raise error.InvalidArgumentValueError(func_name=self.__class__.__name__, arg=args[0])
        # if "label" in kwargs:
        #     self.label = kwargs["label"]

    ###########################################################################
    ### Data I/O

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

    ###########################################################################
    ### Taxon Management

    def reindex_subcomponent_taxa(self):
        """
        Synchronizes `Taxon` objects of map to `taxon_namespace` of self.
        """
        ti_mutable = self.taxon_namespace._is_mutable
        self.taxon_namespace._is_mutable = True
        new_map = TaxonCharacterSequenceMap()
        for taxon, seq in self._taxon_sequence_map.items():
            taxon = self.taxon_namespace.require_taxon(label=taxon.label)
            new_map[taxon] = seq
        self.taxon_namespace._is_mutable = ti_mutable
        self._taxon_sequence_map = new_map

    def update_taxon_namespace(self):
        """
        Updates local taxa block by adding taxa not already managed.
        Mainly for use after map extension
        """
        assert self.taxon_namespace is not None
        for taxon in self._taxon_sequence_map:
            if taxon not in self.taxon_namespace:
                self.taxon_namespace.add(taxon)

    def prune_taxa(self, taxa, update_taxon_namespace=False):
        """
        Removes given taxa from matrix. If `preserve_taxon_namespace` is
        `True`, then the taxa are removed from the associated TaxonNamespace
        object as well. Otherwise this is not modified (default).
        """
        for taxon in taxa:
            if taxon in self._taxon_sequence_map:
                del self._taxon_sequence_map[taxon]
                if update_taxon_namespace and taxon in self.taxon_namespace:
                    self.taxon_namespace.remove(taxon)

    ###########################################################################
    ### Sequence CRUD

    def _resolve_key(self, key):
        """
        Resolves map access key into :class:`Taxon` instance.

        If `key` is integer, assumed to be taxon index.
        If `key` string, assumed to be taxon label.
        Otherwise, assumed to be :class:`Taxon` instance directly.
        """
        if isinstance(key, int):
            if key >= 0 and key < len(self.taxon_namespace):
                taxon = self.taxon_namespace[key]
            else:
                raise IndexError(key)
        elif isinstance(key, str):
            taxon = self.taxon_namespace.get_taxon(label=key)
            if taxon is None:
                raise KeyError(key)
        else:
            taxon = key
        return taxon

    def __getitem__(self, key):
        "Dictionary interface implementation for direct access to character map."
        taxon = self._resolve_key(key)
        try:
            return self._taxon_sequence_map[taxon]
        except KeyError:
            return self.new_sequence(taxon)

    def __setitem__(self, key, values):
        "Dictionary interface implementation for direct access to character map."
        taxon = self._resolve_key(key)
        if taxon not in self.taxon_namespace:
            raise ValueError(key)
        if not isinstance(values, self.__class__.character_sequence_type):
            values = self.__class__.character_sequence_type(values)
        self._taxon_sequence_map[taxon] = values

    def new_sequence(self, taxon, data=None):
        if taxon in self._taxon_sequence_map:
            raise ValueError("Character data vector for taxon {} already exists".format(taxon))
        if taxon not in self.taxon_namespace:
            raise ValueError("Taxon {} is not in object taxon namespace".format(taxon))
        cv = self.__class__.character_sequence_type(data)
        self._taxon_sequence_map[taxon] = cv
        return cv

    def __contains__(self, key):
        "Returns true if character map has `key`"
        return self._taxon_sequence_map.__contains__(key)

    def __delitem__(self, key):
        "Remove item from character map with specified key."
        return self._taxon_sequence_map.__delitem__(key)

    def clear(self):
        "Deletes all items from the character map dictionary."
        self._taxon_sequence_map.clear()

    def has_key(self, key):
        "Returns true if character map has key, regardless of case."
        return key in self._taxon_sequence_map

    ###########################################################################
    ### Column-based Operations

    def remove_column(self, index):
        raise NotImplementedError()

    def append_column(self, index):
        raise NotImplementedError()

    def insert_column(self, index, values):
        raise NotImplementedError()

    ###########################################################################
    ### Sequence Iteration

    def __iter__(self):
        "Returns an iterator over character map's ordered keys."
        for t in self.taxon_namespace:
            if t in self._taxon_sequence_map:
                yield t

    # def iterkeys(self):
    #     "Dictionary interface implementation for direct access to character map."
    #     for t in self.taxon_namespace:
    #         if t in self._taxon_sequence_map:
    #             yield t

    # def itervalues(self):
    #     "Dictionary interface implementation for direct access to character map."
    #     for t in self.taxon_namespace:
    #         if t in self._taxon_sequence_map:
    #             yield self._taxon_sequence_map[t]

    # def items(self):
    #     "Returns character map key, value pairs in key-order."
    #     return [(t, self._taxon_sequence_map[t]) for t in self.taxon_namespace if t in self._taxon_seq_map]

    # def values(self):
    #     "Returns list of values."
    #     return [self._taxon_sequence_map[t] for t in self.taxon_namespace if t in self._taxon_seq_map]

    # def pop(self, key, alt_val=None):
    #     "a.pop(k[, x]):  a[k] if k in a, else x (and remove k)"
    #     return self._taxon_sequence_map.pop(key, alt_val)

    # def popitem(self):
    #     "a.popitem()  remove and last (key, value) pair"
    #     return self._taxon_sequence_map.popitem()

    # def keys(self):
    #     "Returns a copy of the ordered list of character map keys."
    #     return list(self._taxon_sequence_map.keys())

    ###########################################################################
    ### Metrics

    def __len__(self):
        "Dictionary interface implementation for direct access to character map."
        return len(self._taxon_sequence_map)

    def _get_vector_size(self):
        """
        Returns number of characters in *first* sequence.
        """
        if len(self):
            # yuck, but len(self.values())
            # means we have to create and populatre a list ...
            return len(self[next(iter(self._taxon_sequence_map))])
        else:
            return 0
    vector_size = property(_get_vector_size, None, None, __doc__)

    def _get_max_vector_size(self):
        """
        Returns maximum number of characters across all sequences.
        """
        max_len = 0
        for k in self:
            if len(self[k]) > max_len:
                max_len  = len(self._taxon_sequence_map[k])
        return max_len
    max_vector_size = property(_get_max_vector_size, None, None, __doc__)

    ###########################################################################
    ### Mass/Bulk Operations

    def fill(self, state):
        """
        Pads out all sequences in `self` by adding `state` to each sequence
        that is less then the length of the longest sequence.

        Parameters
        ----------
        state : object
            A valid state (e.g., a numeric value for continuous characters, or a :class:`StateIdentity`
            for discrete character).
        """
        max_size = self.max_vector_size
        for k in self:
            v = self[k]
            while len(v) < max_size:
                v.append(state)

    def extend_characters(self, other_matrix):
        """
        Extends this matrix by adding characters from sequences of taxa
        in given matrix to sequences of taxa with correspond labels in
        this one. Taxa in the second matrix that do not exist in the
        current one are ignored.
        """
        self._taxon_sequence_map.extend_characters(other_matrix.taxon_seq_map)

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
        self._taxon_sequence_map.extend(other_map,
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
        self._taxon_sequence_map.extend(other_matrix.taxon_seq_map,
            overwrite_existing=overwrite_existing,
            extend_existing=extend_existing)
        self.update_taxon_namespace()

    ###########################################################################
    ### Character Subset Management

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

    ###########################################################################
    ### Export

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

    def vectors(self):
        "Returns list of vectors.        "
        if self.taxon_namespace is not None and self._taxon_sequence_map is not None:
            if len(self._taxon_sequence_map) > 0:
                return [self._taxon_sequence_map[t] for t in self.taxon_namespace]
            return []
        return None

    def create_taxon_to_state_set_map(self, char_indices=None):
        """Returns a dictionary that maps taxon objects to lists of sets of state
        indices
        if `char_indices` is not None it should be a iterable collection of
        character indices to include.
        """
        taxon_to_state_indices = {}
        for t in self._taxon_sequence_map.keys():
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

    ###########################################################################
    ### Representation

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
                           len(self._taxon_sequence_map[t])))

        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    ###########################################################################
    ### Legacy

    def _get_taxon_seq_map(self):
        error.dump_stack(msg)
        warnings.warn(
                "All methods and features of 'CharacterMatrix.taxon_seq_map' have"
                " been integrated directly into 'CharacterMatrix', or otherwise"
                " replaced entirely")
        return self
    taxon_seq_map = property(_get_taxon_seq_map)

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
        self._default_state_alphabet = None
        if len(args) > 0:
            self.clone_from(*args)

    def _get_default_state_alphabet(self):
        if self._default_state_alphabet is not None:
            return self._default_state_alphabet
        elif len(self.state_alphabets) == 1:
            return self.state_alphabets[0]
        elif len(self.state_alphabets) > 1:
            raise TypeError("Multiple state alphabets defined for this matrix with no default specified")
        elif len(self.state_alphabets) == 0:
            raise TypeError("No state alphabets defined for this matrix")
        return None
    def _set_default_state_alphabet(self, s):
        if s not in self.state_alphabets:
            self.state_alphabets.append(s)
        self._default_state_alphabet = s
    default_state_alphabet = property(_get_default_state_alphabet, _set_default_state_alphabet)

    def append_taxon_sequence(self, taxon, state_symbols):
        if taxon not in self:
            self[taxon] = CharacterSequence()
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
        for vi, vec in enumerate(self._taxon_sequence_map.values()):
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
    matrix_type = get_char_matrix_type(data_type=data_type)
    m = matrix_type(**kwargs)
    return m
