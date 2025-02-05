#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
Character and character-sequence data structures.
"""

import warnings
import copy
from io import StringIO
import math
import collections
from dendropy.utility import textprocessing
from dendropy.utility import error
from dendropy.utility import deprecate
from dendropy.utility import container
from dendropy.datamodel import charstatemodel
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
    A character format or type of a particular column: i.e., maps a particular
    set of character state definitions to a column in a character matrix.
    """

    def __init__(self,
            label=None,
            state_alphabet=None):
        basemodel.DataObject.__init__(self, label=label)
        self._state_alphabet = None
        self.state_alphabet = state_alphabet

    def _get_state_alphabet(self):
        """
        The |StateAlphabet| representing the state alphabet for this
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
## CharacterDataSequence

class CharacterDataSequence(
        basemodel.Annotable,
        ):
    """
    A sequence of character values or values for a particular taxon or entry in
    a data matrix.

    Objects of this class can be (almost) treated as simple lists, where the
    elements are the values of characters (typically, real values in the case
    of continuous data, and special instances of |StateIdentity| objects in the
    case of discrete data.

    Character type data (represented by |CharacterType| instances) and metadata
    annotations (represented by |AnnotationSet| instances), if any, are
    maintained in a parallel list that need to be accessed separately using the
    index of the value to which the data correspond. So, for example, the
    |AnnotationSet| object containing the metadata annotations for the first
    value in a sequence, ``s[0]``, is available through
    ``s.annotations_at(0)``, while the character type information for that
    first element is available through ``s.character_type_at(0)`` and can be
    set through ``s.set_character_type_at(0, c)``.

    In most cases where metadata annotations and character type information are
    not needed, treating objects of this class as a simple list provides all
    the functionality needed. Where metadata annotations or character type
    information are required, all the standard list mutation methods (e.g.,
    ``CharacterDataSequence.insert``, ``CharacterDataSequence.append``,
    ``CharacterDataSequence.extend``) also take optional ``character_type``
    and ``character_annotations`` argument in addition to the primary
    ``character_value`` argument, thus allowing for setting of the value,
    character type, and annotation set simultaneously. While iteration over
    character values are available through the standard list iteration
    interface, the method ``CharacterDataSequence.cell_iter()`` provides for
    iterating over ``<character-value, character-type,
    character-annotation-set>`` triplets.

    """

    ###############################################################################
    ## Life-cycle

    def __init__(self,
            character_values=None,
            character_types=None,
            character_annotations=None):
        """
        Parameters
        ----------
        character_values : iterable of values
            A set of values for this sequence.
        """
        self._character_values = []
        self._character_types = []
        self._character_annotations = []
        if character_values:
            self.extend(
                    character_values=character_values,
                    character_types=character_types,
                    character_annotations=character_annotations)

    ###############################################################################
    ## Life-cycle

    # def values(self):
    #     return list(self._character_values)

    def values(self):
        """
        Returns list of values of this vector.

        Returns
        -------
        v : list
            List of values making up this vector.
        """
        return self._character_values

    def symbols_as_list(self):
        """
        Returns list of string representation of values of this vector.

        Returns
        -------
        v : list
            List of string representation of values making up this vector.
        """
        return list(str(cs) for cs in self._character_values)

    def symbols_as_string(self, sep=""):
        """
        Returns values of this vector as a single string, with individual value
        elements separated by ``sep``.

        Returns
        -------
        s : string
            String representation of values making up this vector.
        """
        return sep.join(str(cs) for cs in self._character_values)

    def __str__(self):
        return self.symbols_as_string()

    def append(self, character_value, character_type=None, character_annotations=None):
        """
        Adds a value to ``self``.

        Parameters
        ----------
        character_value : object
            Value to be stored.
        character_type : |CharacterType|
            Description of character value.
        character_annotations : |AnnotationSet|
            Metadata annotations associated with this character.
        """
        self._character_values.append(character_value)
        self._character_types.append(character_type)
        self._character_annotations.append(character_annotations)

    def extend(self, character_values, character_types=None, character_annotations=None):
        """
        Extends ``self`` with values.

        Parameters
        ----------
        character_values : iterable of objects
            Values to be stored.
        character_types : iterable of |CharacterType| objects
            Descriptions of character values.
        character_annotations : iterable |AnnotationSet| objects
            Metadata annotations associated with characters.
        """
        self._character_values.extend(character_values)
        if character_types is None:
            self._character_types.extend( [None] * len(character_values) )
        else:
            assert len(character_types) == len(character_values)
            self._character_types.extend(character_types)
        if character_annotations is None:
            self._character_annotations.extend( [None] * len(character_values) )
        else:
            assert len(character_annotations) == len(character_values)
            self._character_annotations.extend(character_annotations)

    def __len__(self):
        return len(self._character_values)

    def __getitem__(self, idx):
        return self._character_values[idx]

    def __setitem__(self, idx, value):
        self._character_values[idx] = value

    def __iter__(self):
        return self.__next__()

    def __next__(self):
        for v in self._character_values:
            yield v

    def cell_iter(self):
        """
        Iterate over triplets of character values and associated
        |CharacterType| and |AnnotationSet| instances.
        """
        for v, t, a in zip(self._character_values, self._character_types, self._character_annotations):
            yield v, t, a

    def __delitem__(self, idx):
        del self._character_values[idx]
        del self._character_types[idx]
        del self._character_annotations[idx]

    def set_at(self, idx, character_value, character_type=None, character_annotations=None):
        """
        Set value and associated character type and metadata annotations for
        element at ``idx``.

        Parameters
        ----------
        idx : integer
            Index of element to set.
        character_value : object
            Value to be stored.
        character_type : |CharacterType|
            Description of character value.
        character_annotations : |AnnotationSet|
            Metadata annotations associated with this character.
        """
        to_add = (idx+1) - len(self._character_values)
        while to_add > 0:
            self.append(None)
            to_add -= 1
        self._character_values[idx] = character_value
        self._character_types[idx] = character_type
        self._character_annotations[idx] = character_annotations

    def insert(self, idx, character_value, character_type=None, character_annotations=None):
        """
        Insert value and associated character type and metadata annotations for
        element at ``idx``.

        Parameters
        ----------
        idx : integer
            Index of element to set.
        character_value : object
            Value to be stored.
        character_type : |CharacterType|
            Description of character value.
        character_annotations : |AnnotationSet|
            Metadata annotations associated with this character.
        """
        self._character_values.insert(idx, character_value)
        self._character_types.insert(idx, character_type)
        self._character_annotations.insert(idx, character_annotations)

    def value_at(self, idx):
        """
        Return value of character at ``idx``.

        Parameters
        ----------
        idx : integer
            Index of element value to return.

        Returns
        -------
        c : object
            Value of character at index ``idx``.
        """
        return self._character_values[idx]

    def character_type_at(self, idx):
        """
        Return type of character at ``idx``.

        Parameters
        ----------
        idx : integer
            Index of element character type to return.

        Returns
        -------
        c : |CharacterType|
            |CharacterType| associated with character index ``idx``.
        """
        return self._character_types[idx]

    def annotations_at(self, idx):
        """
        Return metadata annotations of character at ``idx``.

        Parameters
        ----------
        idx : integer
            Index of element annotations to return.

        Returns
        -------
        c : |AnnotationSet|
            |AnnotationSet| representing metadata annotations of character at index ``idx``.
        """
        if self._character_annotations[idx] is None:
            self._character_annotations[idx] = basemodel.AnnotationSet(self._character_types[idx])
        return self._character_annotations[idx]

    def has_annotations_at(self, idx):
        """
        Return |True| if character at ``idx`` has metadata annotations.

        Parameters
        ----------
        idx : integer
            Index of element annotations to check.

        Returns
        -------
        b : bool
            |True| if character at ``idx`` has metadata annotations, |False|
            otherwise.
        """
        return not self._character_annotations[idx] is None

    def set_character_type_at(self, idx, character_type):
        """
        Set type of character at ``idx``.

        Parameters
        ----------
        idx : integer
            Index of element character type to set.
        """
        self._character_types[idx] = character_type

    def set_annotations_at(self, idx, annotations):
        """
        Set metadata annotations of character at ``idx``.

        Parameters
        ----------
        idx : integer
            Index of element annotations to set.
        """
        self._character_annotations[idx] = annotations

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
        Parameters
        ----------
        label: str
            Name of this subset.
        character_indices: iterable of ``int``
            Iterable of 0-based (integer) indices of column positions that
            constitute this subset.

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
        basemodel.Deserializable,
        basemodel.NonMultiReadable,
        basemodel.Serializable,
        basemodel.DataObject):
    """
    A data structure that manages assocation of operational taxononomic unit
    concepts to sequences of character state identities or values.

    This is a base class that provides general functionality; derived classes
    specialize for particular data types. You will not be using the class
    directly, but rather one of the derived classes below, specialized for data
    types such as DNA, RNA, continuous, etc.

    This class and derived classes behave like a dictionary where the keys are
    |Taxon| objects and the values are `CharacterDataSequence` objects. Access
    to sequences based on taxon labels as well as indexes are also provided.
    Numerous methods are provided to manipulate and iterate over sequences.
    Character partitions can be managed through `CharacterSubset` objects,
    while management of detailed metadata on character types are available
    through |CharacterType| objects.

    Objects can be instantiated by reading data from external sources through
    the usual ``get_from_stream()``, ``get_from_path()``, or
    ``get_from_string()`` functions. In addition, a single matrix object can be
    instantiated from multiple matrices (``concatenate()``) or data sources
    (``concatenate_from_paths``).

    A range of methods also exist for importing data from another matrix object.
    These vary depending on how "new" and "existing" are treated.  A "new"
    sequence is a sequence in the other matrix associated with a |Taxon|
    object for which there is no sequence defined in the current matrix.  An
    "existing" sequence is a sequence in the other matrix associated with a
    |Taxon| object for which there *is* a sequence defined in the
    current matrix.

    +---------------------------------+---------------------------------------------+--------------------------------------------+
    |                                 | New Sequences: IGNORED                      | New Sequences: ADDED                       |
    +=================================+=============================================+============================================+
    | Existing Sequences: IGNORED     | [NO-OP]                                     | :meth:`CharacterMatrix.add_sequences()`    |
    +---------------------------------+---------------------------------------------+--------------------------------------------+
    | Existing Sequences: OVERWRITTEN | :meth:`CharacterMatrix.replace_sequences()` | :meth:`CharacterMatrix.update_sequences()` |
    +---------------------------------+---------------------------------------------+--------------------------------------------+
    | Existing Sequences: EXTENDED    | :meth:`CharacterMatrix.extend_sequences()`  | :meth:`CharacterMatrix.extend_matrix()`    |
    +---------------------------------+---------------------------------------------+--------------------------------------------+

    If character subsets have been defined, these subsets can be exported to independent matrices.

    """

    ###########################################################################
    ### Class Variables

    data_type = None
    character_sequence_type = CharacterDataSequence

    ###########################################################################
    ### Factory (Class) Methods

    @classmethod
    def _parse_and_create_from_stream(cls,
            stream,
            schema,
            matrix_offset=0,
            **kwargs):
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        if taxon_namespace is None:
            taxon_namespace = taxonmodel.TaxonNamespace(is_case_sensitive=kwargs.get("case_sensitive_taxon_labels", False))
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
                state_alphabet_factory=charstatemodel.StateAlphabet,
                global_annotations_target=None)
        if len(char_matrices) == 0:
            raise ValueError("No character data in data source")
        char_matrix = char_matrices[matrix_offset]
        if char_matrix.data_type != cls.data_type:
            raise ValueError(
                "Data source (at offset {}) is of type '{}', "
                "but current CharacterMatrix is of type '{}'.".format(
                    matrix_offset,
                    char_matrix.data_type,
                    cls.data_type))
        return char_matrix

    @classmethod
    def get(cls, **kwargs):
        r"""
        Instantiate and return a *new* character matrix object from a data source.

        **Mandatory Source-Specification Keyword Argument (Exactly One of the Following Required):**

            - **file** (*file*) -- File or file-like object of data opened for reading.
            - **path** (*str*) -- Path to file of data.
            - **url** (*str*) -- URL of data.
            - **data** (*str*) -- Data given directly.

        **Mandatory Schema-Specification Keyword Argument:**

            - **schema** (*str*) -- Identifier of format of data given by the
              "``file``", "``path``", "``data``", or "``url``" argument
              specified above: ":doc:`fasta </schemas/fasta>`", ":doc:`nexus
              </schemas/nexus>`", or ":doc:`nexml </schemas/nexml>`",
              ":doc:`phylip </schemas/phylip>`", etc.
              See "|Schemas|" for more details.

        **Optional General Keyword Arguments:**

            - **label** (*str*) -- Name or identifier to be assigned to the new
              object; if not given, will be assigned the one specified in the
              data source, or |None| otherwise.
            - **taxon_namespace** (|TaxonNamespace|) -- The |TaxonNamespace|
              instance to use to :doc:`manage the taxon names </primer/taxa>`.
              If not specified, a new one will be created.
            - **matrix_offset** (*int*) -- 0-based index of character block or
              matrix in source to be parsed. If not specified then the
              first matrix (offset = 0) is assumed.
            - **ignore_unrecognized_keyword_arguments** (*bool*) -- If |True|,
              then unsupported or unrecognized keyword arguments will not
              result in an error. Default is |False|: unsupported keyword
              arguments will result in an error.

        **Optional Schema-Specific Keyword Arguments:**

        These provide control over how the data is interpreted and
        processed, and supported argument names and values depend on
        the schema as specified by the value passed as the "``schema``"
        argument. See "|Schemas|" for more details.

        **Examples:**

        ::

            dna1 = dendropy.DnaCharacterMatrix.get(
                    file=open("pythonidae.fasta"),
                    schema="fasta")
            dna2 = dendropy.DnaCharacterMatrix.get(
                    url="http://purl.org/phylo/treebase/phylows/matrix/TB2:M2610?format=nexus",
                    schema="nexus")
            aa1 = dendropy.ProteinCharacterMatrix.get(
                    file=open("pythonidae.dat"),
                    schema="phylip")
            std1 = dendropy.StandardCharacterMatrix.get(
                    path="python_morph.nex",
                    schema="nexus")
            std2 = dendropy.StandardCharacterMatrix.get(
                    data=">t1\n01011\n\n>t2\n11100",
                    schema="fasta")

        """
        return cls._get_from(**kwargs)

    @classmethod
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
                raise ValueError("Different ``taxon_namespace`` references in matrices to be merged")
            if len(cm) != len(taxon_namespace):
                raise ValueError("Number of sequences not equal to the number of taxa")
            if len(cm) != nseqs:
                raise ValueError("Different number of sequences across alignments: %d (expecting %d based on first matrix)" % (len(cm), nseqs))
            v1 = len(cm[0])
            for t, s in cm.items():
                if len(s) != v1:
                    raise ValueError("Unequal length sequences in character matrix %d".format(cidx+1))
            concatenated_chars.extend_matrix(cm)
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

    @classmethod
    def concatenate_from_streams(cls, streams, schema, **kwargs):
        """
        Read a character matrix from each file object given in ``streams``,
        assuming data format/schema ``schema``, and passing any keyword arguments
        down to the underlying specialized reader. Merge the character matrices
        and return the combined character matrix. Component parts will be
        recorded as character subsets.
        """
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        if taxon_namespace is None:
            taxon_namespace = taxonmodel.TaxonNamespace(is_case_sensitive=kwargs.get("case_sensitive_taxon_labels", False))
        kwargs["taxon_namespace"] = taxon_namespace
        char_matrices = []
        for stream in streams:
            char_matrices.append(cls.get_from_stream(stream,
                schema=schema, **kwargs))
        return cls.concatenate(char_matrices)

    @classmethod
    def concatenate_from_paths(cls, paths, schema, **kwargs):
        """
        Read a character matrix from each file path given in ``paths``, assuming
        data format/schema ``schema``, and passing any keyword arguments down to
        the underlying specialized reader. Merge the and return the combined
        character matrix. Component parts will be recorded as character
        subsets.
        """
        try:
            streams = [open(path, "rU") for path in paths]
        except ValueError:
            streams = [open(path, "r") for path in paths]
        return cls.concatenate_from_streams(streams, schema, **kwargs)

    @classmethod
    def from_dict(cls,
            source_dict,
            char_matrix=None,
            case_sensitive_taxon_labels=False,
            **kwargs):
        r"""
        Populates character matrix from dictionary (or similar mapping type),
        creating |Taxon| objects and sequences as needed.

        Keys must be strings representing labels |Taxon| objects or
        |Taxon| objects directly. If key is specified as string, then it
        will be dereferenced to the first existing |Taxon| object in the
        current taxon namespace with the same label. If no such |Taxon|
        object can be found, then a new |Taxon| object is created and
        added to the current namespace. If a key is specified as a
        |Taxon| object, then this is used directly. If it is not in the
        current taxon namespace, it will be added.

        Values are the sequences (more generally, iterable of values).  If
        values are of type `CharacterDataSequence`, then they are added
        as-is.  Otherwise `CharacterDataSequence` instances are
        created for them. Values may be coerced into types compatible with
        particular matrices. The classmethod ``coerce_values()`` will be
        called for this.

        Examples
        --------

        The following creates a |DnaCharacterMatrix| instance with three
        sequences::

            d = {
                    "s1" : "TCCAA",
                    "s2" : "TGCAA",
                    "s3" : "TG-AA",
            }
            dna = DnaCharacterMatrix.from_dict(d)

        Three |Taxon| objects will be created, corresponding to the
        labels 's1', 's2', 's3'. Each associated string sequence will be
        converted to a `CharacterDataSequence`, with each symbol ("A", "C",
        etc.) being replaced by the DNA state represented by the symbol.

        Parameters
        ----------
        source_dict : dict or other mapping type
            Keys must be strings representing labels |Taxon| objects or
            |Taxon| objects directly. Values are sequences. See above
            for details.
        char_matrix : |CharacterMatrix|
            Instance of |CharacterMatrix| to populate with data. If not
            specified, a new one will be created using keyword arguments
            specified by ``kwargs``.
        case_sensitive_taxon_labels : boolean
            If |True|, matching of string labels specified as keys in ``d`` will
            be matched to |Taxon| objects in current taxon namespace
            with case being respected. If |False|, then case will be ignored.
        \*\*kwargs : keyword arguments, optional
            Keyword arguments to be passed to constructor of
            |CharacterMatrix| when creating new instance to populate, if
            no target instance is provided via ``char_matrix``.

        Returns
        -------
        char_matrix : |CharacterMatrix|
            |CharacterMatrix| populated by data from ``d``.
        """
        if char_matrix is None:
            char_matrix = cls(**kwargs)
        for key in source_dict:
            if textprocessing.is_str_type(key):
                taxon = char_matrix.taxon_namespace.require_taxon(key,
                        is_case_sensitive=case_sensitive_taxon_labels)
            else:
                taxon = key
                if taxon not in char_matrix.taxon_namespace:
                    char_matrix.taxon_namespace.add_taxon(taxon)
            s = char_matrix.coerce_values(source_dict[key])
            char_matrix[taxon] = s
        return char_matrix

    ###########################################################################
    ### Lifecycle and Identity

    def __init__(self, *args, **kwargs):
        if len(args) > 1:
            # only allow 1 positional argument
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        elif len(args) == 1 and isinstance(args[0], CharacterMatrix):
            self._clone_from(args[0], kwargs)
        else:
            basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
            taxonmodel.TaxonNamespaceAssociated.__init__(self,
                    taxon_namespace=taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None))
            self._taxon_sequence_map = {}
            self.character_types = []
            self.comments = []
            self.character_subsets = container.OrderedCaselessDict()
            if len(args) == 1:
                # takes care of all possible initializations, including. e.g.,
                # tuples and so on
                d = collections.OrderedDict(args[0])
                self.__class__.from_dict(d, char_matrix=self)
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def _clone_from(self, src, kwargs_dict):
        # super(Tree, self).__init__()
        memo = {}
        # memo[id(tree)] = self
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs_dict, src.taxon_namespace)
        memo[id(src.taxon_namespace)] = taxon_namespace
        if taxon_namespace is not src.taxon_namespace:
            for t1 in src.taxon_namespace:
                t2 = taxon_namespace.require_taxon(label=t1.label)
                memo[id(t1)] = t2
        else:
            for t1 in src.taxon_namespace:
                memo[id(t1)] = t1
        t = copy.deepcopy(src, memo)
        self.__dict__ = t.__dict__
        self.label = kwargs_dict.pop("label", src.label)
        return self

    def __copy__(self):
        other = self.__class__(label=self.label,
            taxon_namespace=self.taxon_namespace)
        for taxon in self._taxon_sequence_map:
            # other._taxon_sequence_map[taxon] = self.__class__.character_sequence_type(self._taxon_sequence_map[taxon])
            other._taxon_sequence_map[taxon] = self._taxon_sequence_map[taxon]
        memo = {}
        memo[id(self)] = other
        other.deep_copy_annotations_from(self, memo)
        return other

    def taxon_namespace_scoped_copy(self, memo=None):
        if memo is None:
            memo = {}
        # this populates ``memo`` with references to the
        # the TaxonNamespace and Taxon objects
        self.taxon_namespace.populate_memo_for_taxon_namespace_scoped_copy(memo)
        return self.__deepcopy__(memo=memo)

    def __deepcopy__(self, memo=None):
        return basemodel.Annotable.__deepcopy__(self, memo=memo)

    ###########################################################################
    ### Data I/O

    # def _parse_and_add_from_stream(self, stream, schema, **kwargs):
    #     """
    #     Populates objects of this type from ``schema``-formatted
    #     data in the file-like object source ``stream``, *replacing*
    #     all current data. If multiple character matrices are in the data
    #     source, a 0-based index of the character matrix to use can
    #     be specified using the ``matrix_offset`` keyword (defaults to 0, i.e., first
    #     character matrix).
    #     """
    #     warnings.warn("Repopulating a CharacterMatrix is now deprecated. Instantiate a new instance from the source instead.",
    #             DeprecationWarning)
    #     m = self.__class__._parse_and_create_from_stream(stream=stream,
    #             schema=schema,
    #             **kwargs)
    #     return self.clone_from(m)

    def _format_and_write_to_stream(self, stream, schema, **kwargs):
        r"""
        Writes out ``self`` in ``schema`` format to a destination given by
        file-like object ``stream``.

        Parameters
        ----------
        stream : file or file-like object
            Destination for data.
        schema : string
            Must be a recognized character file schema, such as "nexus",
            "phylip", etc, for which a specialized writer is available. If this
            is not implemented for the schema specified, then a
            UnsupportedSchemaError is raised.

        \*\*kwargs : keyword arguments, optional
            Keyword arguments will be passed directly to the writer for the
            specified schema. See documentation for details on keyword
            arguments supported by writers of various schemas.

        """
        writer = dataio.get_writer(schema, **kwargs)
        writer.write_char_matrices([self],
                stream)

    ###########################################################################
    ### Taxon Management

    def reconstruct_taxon_namespace(self,
            unify_taxa_by_label=True,
            taxon_mapping_memo=None):
        """
        See `TaxonNamespaceAssociated.reconstruct_taxon_namespace`.
        """
        if taxon_mapping_memo is None:
            taxon_mapping_memo = {}
        original_taxa = list(self._taxon_sequence_map.keys())
        for original_taxon in original_taxa:
            if unify_taxa_by_label or original_taxon not in self.taxon_namespace:
                t = taxon_mapping_memo.get(original_taxon, None)
                if t is None:
                    # taxon to use not given and
                    # we have not yet created a counterpart
                    if unify_taxa_by_label:
                        # this will force usage of any taxon with
                        # a label that matches the current taxon
                        t = self.taxon_namespace.require_taxon(label=original_taxon.label)
                    else:
                        # this will unconditionally create a new taxon
                        t = self.taxon_namespace.new_taxon(label=original_taxon.label)
                    taxon_mapping_memo[original_taxon] = t
                else:
                    # taxon to use is given by mapping
                    self.taxon_namespace.add_taxon(t)
                if t in self._taxon_sequence_map:
                    raise error.TaxonNamespaceReconstructionError("Multiple sequences for taxon with label '{}'".format(t.label))
                self._taxon_sequence_map[t] = self._taxon_sequence_map[original_taxon]
                del self._taxon_sequence_map[original_taxon]

    def poll_taxa(self, taxa=None):
        """
        Returns a set populated with all of |Taxon| instances associated
        with ``self``.

        Parameters
        ----------
        taxa : set()
            Set to populate. If not specified, a new one will be created.

        Returns
        -------
        taxa : set[|Taxon|]
            Set of taxa associated with ``self``.
        """
        if taxa is None:
            taxa = set()
        for taxon in self._taxon_sequence_map:
            taxa.add(taxon)
        return taxa

    def update_taxon_namespace(self):
        """
        All |Taxon| objects in ``self`` that are not in
        ``self.taxon_namespace`` will be added.
        """
        assert self.taxon_namespace is not None
        for taxon in self._taxon_sequence_map:
            if taxon not in self.taxon_namespace:
                self.taxon_namespace.add_taxon(taxon)

    def reindex_subcomponent_taxa(self):
        """
        Synchronizes |Taxon| objects of map to ``taxon_namespace`` of self.
        """
        raise NotImplementedError("'reindex_subcomponent_taxa()' is no longer supported; use '{}.reconstruct_taxon_namespace()' instead".format(self.__class__.__name__))

    ###########################################################################
    ### Sequence CRUD

    def _resolve_key(self, key):
        """
        Resolves map access key into |Taxon| instance.

        If ``key`` is integer, assumed to be taxon index.
        If ``key`` string, assumed to be taxon label.
        Otherwise, assumed to be |Taxon| instance directly.
        """
        if isinstance(key, int):
            if abs(key) < len(self.taxon_namespace):
                taxon = self.taxon_namespace[key]
            else:
                raise IndexError(key)
        elif textprocessing.is_str_type(key):
            taxon = self.taxon_namespace.get_taxon(label=key)
            if taxon is None:
                raise KeyError(key)
        else:
            taxon = key
        return taxon

    def new_sequence(self, taxon, values=None):
        """
        Creates a new `CharacterDataSequence` associated with |Taxon|
        ``taxon``, and populates it with values in ``values``.

        Parameters
        ----------
        taxon : |Taxon|
            |Taxon| instance with which this sequence is associated.
        values : iterable or |None|
            An initial set of values with which to populate the new character
            sequence.

        Returns
        -------
        s : `CharacterDataSequence`
            A new `CharacterDataSequence` associated with |Taxon|
            ``taxon``.
        """
        if taxon in self._taxon_sequence_map:
            raise ValueError("Character values vector for taxon {} already exists".format(repr(taxon)))
        if taxon not in self.taxon_namespace:
            raise ValueError("Taxon {} is not in object taxon namespace".format(repr(taxon)))
        cv = self.__class__.character_sequence_type(values)
        self._taxon_sequence_map[taxon] = cv
        return cv

    def __getitem__(self, key):
        """
        Retrieves sequence for ``key``, which can be a index or a label of a
        |Taxon| instance in the current taxon namespace, or a
        |Taxon| instance directly.

        If no sequence is currently associated with specified |Taxon|, a
        new one will be created. Note that the |Taxon| object must have
        already been defined in the curent taxon namespace.

        Parameters
        ----------
        key : integer, string, or |Taxon|
            If an integer, assumed to be an index of a |Taxon| object in
            the current |TaxonNamespace| object of ``self.taxon_namespace``.
            If a string, assumed to be a label of a |Taxon| object in
            the current |TaxonNamespace| object of ``self.taxon_namespace``.
            Otherwise, assumed to be |Taxon| instance directly. In all
            cases, the |Taxon| object must be (already) defined in the
            current taxon namespace.

        Returns
        -------
        s : `CharacterDataSequence`
            A sequence associated with the |Taxon| instance referenced
            by ``key``.
        """
        taxon = self._resolve_key(key)
        try:
            return self._taxon_sequence_map[taxon]
        except KeyError:
            return self.new_sequence(taxon)

    def __setitem__(self, key, values):
        """
        Assigns sequence ``values`` to taxon specified by ``key``, which can be a
        index or a label of a |Taxon| instance in the current taxon
        namespace, or a |Taxon| instance directly.

        If no sequence is currently associated with specified |Taxon|, a
        new one will be created.  Note that the |Taxon| object must have
        already been defined in the curent taxon namespace.

        Parameters
        ----------
        key : integer, string, or |Taxon|
            If an integer, assumed to be an index of a |Taxon| object in
            the current |TaxonNamespace| object of ``self.taxon_namespace``.
            If a string, assumed to be a label of a |Taxon| object in
            the current |TaxonNamespace| object of ``self.taxon_namespace``.
            Otherwise, assumed to be |Taxon| instance directly. In all
            cases, the |Taxon| object must be (already) defined in the
            current taxon namespace.

        """
        taxon = self._resolve_key(key)
        if taxon not in self.taxon_namespace:
            raise ValueError(repr(key))
        if not isinstance(values, self.__class__.character_sequence_type):
            values = self.__class__.character_sequence_type(values)
        self._taxon_sequence_map[taxon] = values

    def __contains__(self, key):
        """
        Returns |True| if a sequence associated with ``key`` is in ``self``, or
        |False| otherwise.

        Parameters
        ----------
        key : integer, string, or |Taxon|
            If an integer, assumed to be an index of a |Taxon| object in
            the current |TaxonNamespace| object of ``self.taxon_namespace``.
            If a string, assumed to be a label of a |Taxon| object in
            the current |TaxonNamespace| object of ``self.taxon_namespace``.
            Otherwise, assumed to be |Taxon| instance directly. In all
            cases, the |Taxon| object must be (already) defined in the
            current taxon namespace.

        Returns
        -------
        b : boolean
            |True| if ``key`` is in ``self``; |False| otherwise.
        """
        return self._taxon_sequence_map.__contains__(key)

    def __delitem__(self, key):
        """
        Removes sequence for ``key``, which can be a index or a label of a
        |Taxon| instance in the current taxon namespace, or a
        |Taxon| instance directly.

        Parameters
        ----------
        key : integer, string, or |Taxon|
            If an integer, assumed to be an index of a |Taxon| object in
            the current |TaxonNamespace| object of ``self.taxon_namespace``.
            If a string, assumed to be a label of a |Taxon| object in
            the current |TaxonNamespace| object of ``self.taxon_namespace``.
            Otherwise, assumed to be |Taxon| instance directly. In all
            cases, the |Taxon| object must be (already) defined in the
            current taxon namespace.

        """
        return self._taxon_sequence_map.__delitem__(key)

    def clear(self):
        """
        Removes all sequences from matrix.
        """
        self._taxon_sequence_map.clear()

    def sequences(self):
        """
        List of all sequences in self.

        Returns
        -------
        s : list of `CharacterDataSequence` objects in self

        """
        s = [self[taxon] for taxon in self]
        return s

    def vectors(self):
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'vectors()' will no longer be supported in future releases; use 'sequences()' instead")
        return self.sequences()

    ###########################################################################
    ### Symbol/alphabet management

    def coerce_values(self, values):
        """
        Converts elements of ``values`` to type of matrix.

        This method is called by :meth:`CharacterMatrix.from_dict` to create
        sequences from iterables of values.  This method should be overridden
        by derived classes to ensure that ``values`` consists of types compatible
        with the particular type of matrix. For example, a CharacterMatrix type
        with a fixed state alphabet (such as |DnaCharacterMatrix|) would
        dereference the string elements of ``values`` to return a list of
        |StateIdentity| objects corresponding to the symbols represented
        by the strings.  If there is no value-type conversion done, then
        ``values`` should be returned as-is. If no value-type conversion is
        possible (e.g., when the type of a value is dependent on positionaly
        information), then a TypeError should be raised.

        Parameters
        ----------
        values : iterable
            Iterable of values to be converted.

        Returns
        -------
        v : list of values.
        """
        return values

    ###########################################################################
    ### Sequence Access Iteration

    def __iter__(self):
        "Returns an iterator over character map's ordered keys."
        for t in self.taxon_namespace:
            if t in self._taxon_sequence_map:
                yield t

    def __contains__(self, key):
        if isinstance(key, str):
            for t in self._taxon_sequence_map:
                if t.label == key:
                    return True
            return False
        return key in self._taxon_sequence_map

    def values(self):
        """
        Iterates values (i.e. sequences) in this matrix.
        """
        for t in self:
            yield self[t]

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

    def items(self):
        "Returns character map key, value pairs in key-order."
        for t in self.taxon_namespace:
            if t in self._taxon_sequence_map:
                yield t, self._taxon_sequence_map[t]

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
        """
        Number of sequences in matrix.

        Returns
        -------
        n : Number of sequences in matrix.
        """
        return len(self._taxon_sequence_map)

    def _get_sequence_size(self):
        """
        Number of characters in *first* sequence in matrix.

        Returns
        -------
        n : integer
            Number of sequences in matrix.
        """
        if len(self):
            # yuck, but len(self.values())
            # means we have to create and populate a list ...
            return len(self[next(iter(self._taxon_sequence_map))])
        else:
            return 0
    sequence_size = property(_get_sequence_size, None, None)
    vector_size = property(_get_sequence_size, None, None) # legacy

    def _get_max_sequence_size(self):
        """
        Maximum number of characters across all sequences in matrix.

        Returns
        -------
        n : integer
            Maximum number of characters across all sequences in matrix.
        """
        max_len = 0
        for k in self:
            if len(self[k]) > max_len:
                max_len  = len(self._taxon_sequence_map[k])
        return max_len
    max_sequence_size = property(_get_max_sequence_size, None, None)

    ###########################################################################
    ### Mass/Bulk Operations

    def fill(self, value, size=None, append=True):
        """
        Pads out all sequences in ``self`` by adding ``value`` to each sequence
        until its length is ``size`` long or equal to the length of the longest
        sequence if ``size`` is not specified.

        Parameters
        ----------
        value : object
            A valid value (e.g., a numeric value for continuous characters, or
            a |StateIdentity| for discrete character).
        size : integer or None
            The size (length) up to which the sequences will be padded. If |None|, then
            the maximum (longest) sequence size will be used.
        append : boolean
            If |True| (default), then new values will be added to the end of
            each sequence. If |False|, then new values will be inserted to the
            front of each sequence.
        """
        if size is None:
            size = self.max_sequence_size
        for k in self:
            v = self[k]
            while len(v) < size:
                if append:
                    v.append(value)
                else:
                    v.insert(0, value)
        return size

    def fill_taxa(self):
        """
        Adds a new (empty) sequence for each |Taxon| instance in
        current taxon namespace that does not have a sequence.
        """
        for taxon in self.taxon_namespace:
            if taxon not in self:
                self[taxon] = CharacterDataSequence()

    def pack(self, value=None, size=None, append=True):
        """
        Adds missing sequences for all |Taxon| instances in current
        namespace, and then pads out all sequences in ``self`` by adding ``value``
        to each sequence until its length is ``size`` long or equal to the length
        of the longest sequence if ``size`` is not specified. A combination of
        :meth:`CharacterMatrix.fill_taxa()` and
        :meth:`CharacterMatrix.fill()`.

        Parameters
        ----------
        value : object
            A valid value (e.g., a numeric value for continuous characters, or
            a |StateIdentity| for discrete character).
        size : integer or None
            The size (length) up to which the sequences will be padded. If |None|, then
            the maximum (longest) sequence size will be used.
        append : boolean
            If |True| (default), then new values will be added to the end of
            each sequence. If |False|, then new values will be inserted to the
            front of each sequence.
        """
        self.fill_taxa()
        self.fill(value=value, size=size, append=append)

    def add_sequences(self, other_matrix):
        """
        Adds sequences for |Taxon| objects that are in ``other_matrix`` but not in
        ``self``.

        Parameters
        ----------
        other_matrix : |CharacterMatrix|
            Matrix from which to add sequences.

        Notes
        -----
            1. ``other_matrix`` must be of same type as ``self``.
            2. ``other_matrix`` must have the same |TaxonNamespace| as ``self``.
            3. Each sequence associated with a |Taxon| reference in ``other_matrix``
               but not in ``self`` will be added to ``self`` as a shallow-copy.
            4. All other sequences will be ignored.

        """
        if other_matrix.taxon_namespace is not self.taxon_namespace:
            raise error.TaxonNamespaceIdentityError(self, other_matrix)
        for taxon in other_matrix._taxon_sequence_map:
            if taxon not in self._taxon_sequence_map:
                self._taxon_sequence_map[taxon] = self.__class__.character_sequence_type(other_matrix._taxon_sequence_map[taxon])

    def replace_sequences(self, other_matrix):
        """
        Replaces sequences for |Taxon| objects shared between ``self`` and
        ``other_matrix``.

        Parameters
        ----------
        other_matrix : |CharacterMatrix|
            Matrix from which to replace sequences.

        Notes
        -----
            1. ``other_matrix`` must be of same type as ``self``.
            2. ``other_matrix`` must have the same |TaxonNamespace| as ``self``.
            3. Each sequence in ``self`` associated with a |Taxon| that is
               also represented in ``other_matrix`` will be replaced with a
               shallow-copy of the corresponding sequence from ``other_matrix``.
            4. All other sequences will be ignored.
        """
        if other_matrix.taxon_namespace is not self.taxon_namespace:
            raise error.TaxonNamespaceIdentityError(self, other_matrix)
        for taxon in other_matrix._taxon_sequence_map:
            if taxon in self._taxon_sequence_map:
                self._taxon_sequence_map[taxon] = self.__class__.character_sequence_type(other_matrix._taxon_sequence_map[taxon])

    def update_sequences(self, other_matrix):
        """
        Replaces sequences for |Taxon| objects shared between ``self`` and
        ``other_matrix`` and adds sequences for |Taxon| objects that are
        in ``other_matrix`` but not in ``self``.

        Parameters
        ----------
        other_matrix : |CharacterMatrix|
            Matrix from which to update sequences.

        Notes
        -----
            1. ``other_matrix`` must be of same type as ``self``.
            2. ``other_matrix`` must have the same |TaxonNamespace| as ``self``.
            3. Each sequence associated with a |Taxon| reference in ``other_matrix``
               but not in ``self`` will be added to ``self``.
            4. Each sequence in ``self`` associated with a |Taxon| that is
               also represented in ``other_matrix`` will be replaced with a
               shallow-copy of the corresponding sequence from ``other_matrix``.
        """
        if other_matrix.taxon_namespace is not self.taxon_namespace:
            raise error.TaxonNamespaceIdentityError(self, other_matrix)
        for taxon in other_matrix._taxon_sequence_map:
            self._taxon_sequence_map[taxon] = self.__class__.character_sequence_type(other_matrix._taxon_sequence_map[taxon])

    def extend_sequences(self, other_matrix, is_add_new_sequences=False):
        """
        Extends sequences in ``self`` with characters associated with
        corresponding |Taxon| objects in ``other_matrix``.

        Parameters
        ----------
        other_matrix : |CharacterMatrix|
            Matrix from which to extend sequences.

        Notes
        -----
            1. ``other_matrix`` must be of same type as ``self``.
            2. ``other_matrix`` must have the same |TaxonNamespace| as ``self``.
            3. Each sequence associated with a |Taxon| reference in
               ``other_matrix`` that is also in ``self`` will be appended to the
               sequence currently associated with that |Taxon| reference
               in ``self``.
            4. All other sequences will be ignored.
        """
        if other_matrix.taxon_namespace is not self.taxon_namespace:
            raise error.TaxonNamespaceIdentityError(self, other_matrix)
        for taxon in other_matrix._taxon_sequence_map:
            if taxon not in self._taxon_sequence_map:
                if not is_add_new_sequences:
                    continue
                self._taxon_sequence_map[taxon] = self.__class__.character_sequence_type(other_matrix._taxon_sequence_map[taxon])
            else:
                self._taxon_sequence_map[taxon].extend(other_matrix._taxon_sequence_map[taxon])

    def extend_matrix(self, other_matrix):
        """
        Extends sequences in ``self`` with characters associated with
        corresponding |Taxon| objects in ``other_matrix`` and adds
        sequences for |Taxon| objects that are in ``other_matrix`` but not
        in ``self``.

        Parameters
        ----------
        other_matrix : |CharacterMatrix|
            Matrix from which to extend.

        Notes
        -----
            1. ``other_matrix`` must be of same type as ``self``.
            2. ``other_matrix`` must have the same |TaxonNamespace| as ``self``.
            3. Each sequence associated with a |Taxon| reference in ``other_matrix``
               that is also in ``self`` will be appending
               to the sequence currently associated with that |Taxon|
               reference in ``self``.
            4. Each sequence associated with a |Taxon| reference in
               ``other_matrix`` that is also in ``self`` will replace the sequence
               currently associated with that |Taxon| reference in ``self``.
        """
        if other_matrix.taxon_namespace is not self.taxon_namespace:
            raise error.TaxonNamespaceIdentityError(self, other_matrix)
        for taxon in other_matrix._taxon_sequence_map:
            if taxon in self._taxon_sequence_map:
                self._taxon_sequence_map[taxon].extend(other_matrix._taxon_sequence_map[taxon])
            else:
                self._taxon_sequence_map[taxon]= self.__class__.character_sequence_type(other_matrix._taxon_sequence_map[taxon])

    def remove_sequences(self, taxa):
        """
        Removes sequences associated with |Taxon| instances specified in
        ``taxa``. A KeyError is raised if a |Taxon| instance is
        specified for which there is no associated sequences.

        Parameters
        ----------
        taxa : iterable[|Taxon|]
            List or some other iterable of |Taxon| instances.
        """
        for taxon in taxa:
            del self._taxon_sequence_map[taxon]

    def discard_sequences(self, taxa):
        """
        Removes sequences associated with |Taxon| instances specified in
        ``taxa`` if they exist.

        Parameters
        ----------
        taxa : iterable[|Taxon|]
            List or some other iterable of |Taxon| instances.
        """
        for taxon in taxa:
            try:
                del self._taxon_sequence_map[taxon]
            except KeyError:
                pass

    def keep_sequences(self, taxa):
        """
        Discards all sequences *not* associated with any of the |Taxon| instances.

        Parameters
        ----------
        taxa : iterable[|Taxon|]
            List or some other iterable of |Taxon| instances.
        """
        to_keep = set(taxa)
        for taxon in tuple(self._taxon_sequence_map.keys()):
            if taxon not in to_keep:
                del self._taxon_sequence_map[taxon]

    # def extend_characters(self, other_matrix):
    #     """
    #     DEPRECATED
    #     Extends this matrix by adding characters from sequences of taxa
    #     in given matrix to sequences of taxa with correspond labels in
    #     this one. Taxa in the second matrix that do not exist in the
    #     current one are ignored.
    #     """
    #     self._taxon_sequence_map.extend_characters(other_matrix.taxon_seq_map)

    # def extend_map(self,
    #                   other_map,
    #                   overwrite_existing=False,
    #                   extend_existing=False):
    #     """
    #     DEPRECATED
    #     Extends this matrix by adding taxa and characters from the given
    #     map to this one.  If ``overwrite_existing`` is True and a taxon
    #     in the other map is already present in the current one, then
    #     the sequence associated with the taxon in the second map
    #     replaces the sequence in the current one. If ``extend_existing``
    #     is True and a taxon in the other matrix is already present in
    #     the current one, then the squence map with the taxon in
    #     the second map will be added to the sequence in the current
    #     one. If both are True, then an exception is raised. If neither
    #     are True,  and a taxon in the other map is already present in
    #     the current one, then the sequence is ignored.
    #     """
    #     self._taxon_sequence_map.extend(other_map,
    #         overwrite_existing=overwrite_existing,
    #         extend_existing=extend_existing)
    #     self.update_taxon_namespace()

    # def extend(self,
    #            other_matrix,
    #            overwrite_existing=False,
    #            extend_existing=False):
    #     """
    #     Extends this matrix by adding taxa and characters from the given
    #     matrix to this one.  If ``overwrite_existing`` is True and a taxon
    #     in the other matrix is already present in the current one, then
    #     the sequence associated with the taxon in the second matrix
    #     replaces the sequence in the current one. If ``extend_existing``
    #     is True and a taxon in the other matrix is already present in
    #     the current one, then the sequence associated with the taxon in
    #     the second matrix will be added to the sequence in the current
    #     one. If both are True, then an exception is raised. If neither
    #     are True, and a taxon in the other matrix is already present in
    #     the current one, then the sequence is ignored.
    #     """
    #     self._taxon_sequence_map.extend(other_matrix.taxon_seq_map,
    #         overwrite_existing=overwrite_existing,
    #         extend_existing=extend_existing)
    #     self.update_taxon_namespace()

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
    ### CharacterType Management

    def new_character_type(self, *args, **kwargs):
        return CharacterType(*args, **kwargs)

    ###########################################################################
    ### Export

    def export_character_subset(self, character_subset):
        """
        Returns a new CharacterMatrix (of the same type) consisting only
        of columns given by the CharacterSubset, ``character_subset``.
        Note that this new matrix will still reference the same taxon set.
        """
        if textprocessing.is_str_type(character_subset):
            if character_subset not in self.character_subsets:
                raise KeyError(character_subset)
            else:
                character_subset = self.character_subsets[character_subset]
        return self.export_character_indices(character_subset.character_indices)

    def export_character_indices(self, indices):
        """
        Returns a new CharacterMatrix (of the same type) consisting only
        of columns given by the 0-based indices in ``indices``.
        Note that this new matrix will still reference the same taxon set.
        """
        clone = self.__class__(self)
        # clear out character subsets; otherwise all indices will have to be
        # recalculated, which will require some careful and perhaps arbitrary
        # handling of corner cases
        clone.character_subsets = container.OrderedCaselessDict()
        indices = set(indices)
        for vec in clone.values():
            for cell_idx in range(len(vec)-1, -1, -1):
                if cell_idx not in indices:
                    del(vec[cell_idx])
        return clone

    ###########################################################################
    ### Representation

    def description(self, depth=1, indent=0, itemize="", output=None):
        """
        Returns description of object, up to level ``depth``.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        label = " (%s: '%s')" % (id(self), self.label)
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
        warnings.warn("All methods and features of 'CharacterMatrix.taxon_seq_map' have been integrated directly into 'CharacterMatrix', or otherwise replaced entirely",
                stacklevel=2)
        return self
    taxon_seq_map = property(_get_taxon_seq_map)

###############################################################################
## Specialized Matrices

### Continuous Characters ##################################################

class ContinuousCharacterDataSequence(CharacterDataSequence):
    """
    A sequence of continuous character values for a particular taxon or entry
    in a data matrix. Specializes `CharacterDataSequence` by assuming all
    values are primitive numerics (i.e., either floats or integers) when
    copying or representing self.
    """

    def symbols_as_list(self):
        """
        Returns list of string representation of values of this vector.

        Returns
        -------
        v : list
            List of string representation of values making up this vector.
        """
        return [str(v) for v in self]

    def symbols_as_string(self, sep=" "):
        # different default
        return CharacterDataSequence.symbols_as_string(self, sep=sep)

class ContinuousCharacterMatrix(CharacterMatrix):
    """
    Specializes |CharacterMatrix| for continuous data.

    Sequences stored using |ContinuousCharacterDataSequence|, with values of
    elements assumed to be ``float`` .
    """

    character_sequence_type = ContinuousCharacterDataSequence
    data_type = "continuous"

    def __init__(self, *args, **kwargs):
        CharacterMatrix.__init__(self, *args, **kwargs)

### Discrete Characters ##################################################

class DiscreteCharacterDataSequence(CharacterDataSequence):
    pass

class DiscreteCharacterMatrix(CharacterMatrix):

    character_sequence_type = DiscreteCharacterDataSequence

    data_type = "discrete"

    def __init__(self, *args, **kwargs):
        CharacterMatrix.__init__(self, *args, **kwargs)
        self.state_alphabets = []
        self._default_state_alphabet = None

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
            self[taxon] = CharacterDataSequence()
        for value in state_symbols:
            if textprocessing.is_str_type(value):
                symbol = value
            else:
                symbol = str(value)
            self[taxon].append(self.default_symbol_state_map[symbol])

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
        for vi, vec in enumerate(self._taxon_sequence_map.values()):
            for ci, cell in enumerate(vec):
                vec[ci] = state_alphabet[cell.symbol]
        for ct in self.character_types:
            if ct is not None:
                ct.state_alphabet = state_alphabet
        if purge_other_state_alphabets:
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

    def taxon_state_sets_map(self,
            char_indices=None,
            gaps_as_missing=True,
            gap_state=None,
            no_data_state=None):
        """
        Returns a dictionary that maps taxon objects to lists of sets of
        fundamental state indices.

        Parameters
        ----------

        char_indices : iterable of ints
            An iterable of indexes of characters to include (by column). If not
            given or |None| [default], then all characters are included.

        gaps_as_missing : boolean
            If |True| [default] then gap characters will be treated as missing
            data values. If |False|, then they will be treated as an additional
            (fundamental) state.`

        Returns
        -------
        d : dict
            A dictionary with class:|Taxon| objects as keys and a list of sets
            of fundamental state indexes as values.

            E.g., Given the following matrix of DNA characters:

                T1 AGN
                T2 C-T
                T3 GC?

            Return with ``gaps_as_missing==True`` ::

                {
                    <T1> : [ set([0]), set([2]),        set([0,1,2,3]) ],
                    <T2> : [ set([1]), set([0,1,2,3]),  set([3]) ],
                    <T3> : [ set([2]), set([1]),        set([0,1,2,3]) ],
                }

            Return with ``gaps_as_missing==False`` ::

                {
                    <T1> : [ set([0]), set([2]),        set([0,1,2,3]) ],
                    <T2> : [ set([1]), set([4]),        set([3]) ],
                    <T3> : [ set([2]), set([1]),        set([0,1,2,3,4]) ],
                }

            Note that when gaps are treated as a fundamental state, not only
            does '-' map to a distinct and unique state (4), but '?' (missing
            data) maps to set consisting of all bases *and* the gap
            state, whereas 'N' maps to a set of all bases but not including the
            gap state.

            When gaps are treated as missing, on the other hand, then '?' and
            'N' and '-' all map to the same set, i.e. of all the bases.

        """
        taxon_to_state_indices = {}
        for t in self:
            cdv = self[t]
            if char_indices is None:
                ci = range(len(cdv))
            else:
                ci = char_indices
            v = []
            for char_index in ci:
                state = cdv[char_index]
                if gaps_as_missing:
                    v.append(set(state.fundamental_indexes_with_gaps_as_missing))
                else:
                    v.append(set(state.fundamental_indexes))
            taxon_to_state_indices[t] = v
        return taxon_to_state_indices

    def folded_site_frequency_spectrum(self, is_pad_vector_to_unfolded_length=False):
        r"""
        Returns the folded or minor site/allele frequency spectrum.

        Given $N$ chromosomes, the site frequency spectrum is a vector $(f_0,
        f_1, f_2, ..., f_N)$, where the value $f_i$ is the number of
        sites where $i$ derived alleles are segregating in the sample: 0
        alleles, 1 allele, 2 alleles, etc.

        The *folded* site frequency spectrum is a vector $(f_0, f_1, f_2, ...,
        f_m), m = \ceil{\frac{N}{2}}$, where the values are the number of minor
        alleles in the site.

        Parameters
        ----------
        is_pad_vector_to_unfolded_length: bool
            If False, then the vector length will be $\ceil{\frac{N}{2}}$,
            where $N$ is the number of taxa. Otherwise, by default,
            True, length of vector will be number of taxa + 1, with the
            first element the number of monomorphic sites not contributing to
            the site frequency spectrum.

        Returns
        -------
        v : list[int]
            A vector of integers representing the folded site frequency
            spectrum.
        """
        site_columns = zip(*self.sequences())
        nsites = 0
        if is_pad_vector_to_unfolded_length:
            sfs = [0 for idx in range(len(self._taxon_sequence_map)+1)]
        else:
            sfs = [0 for idx in range(int(math.ceil(len(self._taxon_sequence_map)/2.0))+1)]
        for site in site_columns:
            counter = collections.Counter(site)
            nsites += 1
            if len(counter) == 1:
                sfs[0] += 1
                continue
            del counter[counter.most_common(1)[0][0]]
            sfs[sum(counter.values())] += 1
        assert sum(sfs) == nsites
        return sfs

### Fixed Alphabet Characters ##################################################

class FixedAlphabetCharacterDataSequence(CharacterDataSequence):
    pass

class FixedAlphabetCharacterMatrix(DiscreteCharacterMatrix):

    character_sequence_type = FixedAlphabetCharacterDataSequence
    data_type = "fixed"
    datatype_alphabet = None

    def __init__(self, *args, **kwargs):
        DiscreteCharacterMatrix.__init__(self, *args, **kwargs)
        self.state_alphabets.append(self.__class__.datatype_alphabet)
        self._default_state_alphabet = self.__class__.datatype_alphabet

    def coerce_values(self, values):
        if self.datatype_alphabet is None:
            raise ValueError("'datatype_alphabet' not set")
        return charstatemodel.coerce_to_state_identities(
                state_alphabet=self.datatype_alphabet,
                values=values)

### DNA Characters ##################################################

class DnaCharacterDataSequence(FixedAlphabetCharacterDataSequence):
    pass

class DnaCharacterMatrix(FixedAlphabetCharacterMatrix):
    """
    Specializes |CharacterMatrix| for DNA data.
    """
    character_sequence_type = DnaCharacterDataSequence
    data_type = "dna"
    datatype_alphabet = DNA_STATE_ALPHABET

### RNA Characters ##################################################

class RnaCharacterDataSequence(FixedAlphabetCharacterDataSequence):
    pass

class RnaCharacterMatrix(FixedAlphabetCharacterMatrix):
    """
    Specializes |CharacterMatrix| for DNA data.
    """
    character_sequence_type = RnaCharacterDataSequence
    data_type = "rna"
    datatype_alphabet = RNA_STATE_ALPHABET

### Nucleotide Characters ##################################################

class NucleotideCharacterDataSequence(FixedAlphabetCharacterDataSequence):
    pass

class NucleotideCharacterMatrix(FixedAlphabetCharacterMatrix):
    """
    Specializes |CharacterMatrix| for RNA data.
    """
    character_sequence_type = NucleotideCharacterDataSequence
    data_type = "nucleotide"
    datatype_alphabet = NUCLEOTIDE_STATE_ALPHABET

### Protein Characters ##################################################

class ProteinCharacterDataSequence(FixedAlphabetCharacterDataSequence):
    pass

class ProteinCharacterMatrix(FixedAlphabetCharacterMatrix):
    """
    Specializes |CharacterMatrix| for protein or amino acid data.
    """
    character_sequence_type = ProteinCharacterDataSequence
    data_type = "protein"
    datatype_alphabet = PROTEIN_STATE_ALPHABET

### Restricted Site Characters ##################################################

class RestrictionSitesCharacterDataSequence(FixedAlphabetCharacterDataSequence):
    pass

class RestrictionSitesCharacterMatrix(FixedAlphabetCharacterMatrix):
    """
    Specializes |CharacterMatrix| for restriction site data.
    """
    character_sequence_type = RestrictionSitesCharacterDataSequence
    data_type = "restriction"
    datatype_alphabet = RESTRICTION_SITES_STATE_ALPHABET

### Infinite Sites Characters ##################################################

class InfiniteSitesCharacterDataSequence(FixedAlphabetCharacterDataSequence):
    pass

class InfiniteSitesCharacterMatrix(FixedAlphabetCharacterMatrix):
    """
    Specializes |CharacterMatrix| for infinite sites data.
    """
    character_sequence_type = InfiniteSitesCharacterDataSequence
    data_type = "infinite"
    datatype_alphabet = INFINITE_SITES_STATE_ALPHABET

### Standard Characters ##################################################

class StandardCharacterDataSequence(DiscreteCharacterDataSequence):
    pass

class StandardCharacterMatrix(DiscreteCharacterMatrix):
    """
    Specializes |CharacterMatrix| for "standard" data (i.e., generic discrete
    character data).

    """
    character_sequence_type = StandardCharacterDataSequence

    data_type = "standard"

    def __init__(self, *args, **kwargs):
        """
        A default state alphabet consisting of state symbols of 0-9 will
        automatically be created unless the ``default_state_alphabet=None`` is
        passed in. To specify a different default state alphabet::

            default_state_alphabet=dendropy.new_standard_state_alphabet("abc")
            default_state_alphabet=dendropy.new_standard_state_alphabet("ij")

        """
        if "default_state_alphabet" in kwargs:
            default_state_alphabet = kwargs.pop("default_state_alphabet")
        else:
            default_state_alphabet = charstatemodel.new_standard_state_alphabet()
        DiscreteCharacterMatrix.__init__(self, *args, **kwargs)
        if default_state_alphabet is not None:
            self.default_state_alphabet = default_state_alphabet

    def coerce_values(self, values):
        if self.default_state_alphabet is None:
            raise ValueError("'default_state_alphabet' not set")
        return charstatemodel.coerce_to_state_identities(
                state_alphabet=self.default_state_alphabet,
                values=values)


###############################################################################
## Main Character Matrix Factory Function

data_type_matrix_map = {
    'continuous' : ContinuousCharacterMatrix,
    'dna' : DnaCharacterMatrix,
    'rna' : RnaCharacterMatrix,
    'nucleotide' : NucleotideCharacterMatrix,
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
