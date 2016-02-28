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
This module defines the |DataSet|: a top-level data container object
that manages collections of |TaxonNamespace|, |TreeList|, and
(various kinds of) |CharacterMatrix| objects.
"""

import warnings
import copy
import sys
from dendropy.utility import container
from dendropy.utility import error
from dendropy.utility import deprecate
from dendropy.utility import textprocessing
from dendropy.datamodel import basemodel
from dendropy.datamodel import taxonmodel
from dendropy.datamodel import treecollectionmodel
from dendropy.datamodel import charmatrixmodel
from dendropy.datamodel import charstatemodel
from dendropy import dataio

###############################################################################
## DataSet

class DataSet(
        basemodel.Annotable,
        basemodel.Deserializable,
        basemodel.MultiReadable,
        basemodel.Serializable,
        basemodel.DataObject):
    """
    A phylogenetic data object that coordinates collections of
    |TaxonNamespace|, |TreeList|, and (various kinds of)
    |CharacterMatrix| objects.

    A |DataSet| has three attributes:

        ``taxon_namespaces``
            A list of |TaxonNamespace| objects, each representing
            a distinct namespace for operational taxononomic unit concept
            definitions.

        ``tree_lists``
            A list of |TreeList| objects, each representing a
            collection of |Tree| objects.

        ``char_matrices``
            A list of |CharacterMatrix|-derived objects (e.g.
            |DnaCharacterMatrix|).

    Multiple |TaxonNamespace| objects within a |DataSet| are
    allowed so as to support reading/loading of data from external sources that
    have multiple independent taxon namespaces defined within the same source
    or document (e.g., a Mesquite file with multiple taxa blocks, or a NeXML
    file with multiple OTU sections). Ideally, however, this would not
    be how data is managed. Recommended idiomatic usage would be to use a
    |DataSet| to manage multiple types of data that all share and
    reference the same, single taxon namespace.

    This convention can be enforced by setting the DataSet instance to
    "attached taxon namespace" mode::

        ds = dendropy.DataSet()
        tns = dendropy.TaxonNamespace()
        ds.attach_taxon_namespace(tns)

    After setting this mode, all subsequent data read or created will be
    coerced to use the same, common operational taxonomic unit concept
    namespace.

    Note that unless there is a need to collect and serialize a collection of
    data to the same file or external source, it is probably better
    semantically to use more specific data structures (e.g., a
    |TreeList| object for trees or a |DnaCharacterMatrix|
    object for an alignment). Similarly, when deserializing an external
    data source, if just a single type or collection of data is needed (e.g.,
    the collection of trees from a file that includes both trees and an
    alignment), then it is semantically cleaner to deserialize the data
    into a more specific structure (e.g., a |TreeList| to get all the
    trees). However, when deserializing a mixed external data source
    with, e.g. multiple alignments or trees and one or more alignments, and you
    need to access/use more than a single collection, it is more efficient to
    read the entire data source at once into a |DataSet| object and then
    independently extract the data objects as you need them from the various
    collections.

    """

    def _parse_and_create_from_stream(cls,
            stream,
            schema,
            **kwargs):
        """
        Constructs a new |DataSet| object and populates it with data
        from file-like object ``stream``.
        """
        exclude_trees = kwargs.pop("exclude_trees", False)
        exclude_chars = kwargs.pop("exclude_chars", False)
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        label = kwargs.pop("label", None)
        dataset = DataSet(label=label)
        if taxon_namespace is not None:
            dataset.attached_taxon_namespace = taxon_namespace
        reader = dataio.get_reader(schema, **kwargs)
        reader.read_dataset(
                stream=stream,
                dataset=dataset,
                taxon_namespace=taxon_namespace,
                exclude_trees=exclude_trees,
                exclude_chars=exclude_chars,
                state_alphabet_factory=charstatemodel.StateAlphabet,
                )
        return dataset
    _parse_and_create_from_stream = classmethod(_parse_and_create_from_stream)

    @classmethod
    def get(cls, **kwargs):
        """
        Instantiate and return a *new* |TreeList| object from a data source.

        **Mandatory Source-Specification Keyword Argument (Exactly One Required):**

            - **file** (*file*) -- File or file-like object of data opened for reading.
            - **path** (*str*) -- Path to file of data.
            - **url** (*str*) -- URL of data.
            - **data** (*str*) -- Data given directly.

        **Mandatory Schema-Specification Keyword Argument:**

            - **schema** (*str*) -- Identifier of format of data given by the
              "``file``", "``path``", "``data``", or "``url``" argument
              specified above: ":doc:`fasta </schemas/fasta>`", ":doc:`newick
              </schemas/newick>`", ":doc:`nexus </schemas/nexus>`", or
              ":doc:`nexml </schemas/nexml>`", ":doc:`phylip
              </schemas/phylip>`", etc. See "|Schemas|" for more details.

        **Optional General Keyword Arguments:**

            - **exclude_trees** (*bool*) -- If |True|, then all tree data in the data
              source will be skipped.
            - **exclude_chars** (*bool*) -- If |True|, then all character
              data in the data source will be skipped.
            - **taxon_namespace** (|TaxonNamespace|) -- The |TaxonNamespace|
              instance to use to :doc:`manage the taxon names </primer/taxa>`.
              If not specified, a new one will be created.
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

            dataset1 = dendropy.DataSet.get(
                    path="pythonidae.chars_and_trees.nex",
                    schema="nexus")
            dataset2 = dendropy.DataSet.get(
                    url="http://purl.org/phylo/treebase/phylows/study/TB2:S1925?format=nexml",
                    schema="nexml")

        """
        return cls._get_from(**kwargs)

    ###########################################################################
    ### Lifecycle and Identity

    def __init__(self, *args, **kwargs):
        """
        The constructor can take one argument. This can either be another
        |DataSet| instance or an iterable of |TaxonNamespace|,
        |TreeList|, or |CharacterMatrix|-derived instances.

        In the former case, the newly-constructed |DataSet| will be a
        shallow-copy clone of the argument.

        In the latter case, the newly-constructed |DataSet| will have
        the elements of the iterable added to the respective collections
        (``taxon_namespaces``, ``tree_lists``, or ``char_matrices``, as
        appropriate). This is essentially like calling :meth:`DataSet.add()`
        on each element separately.
        """
        if len(args) > 1:
            # only allow 1 positional argument
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        if "stream" in kwargs or "schema" in kwargs:
            raise TypeError("Constructing from an external stream is no longer supported: use the factory method 'DataSet.get_from_stream()'")
        elif len(args) == 1 and isinstance(args[0], DataSet):
            self._clone_from(args[0], kwargs)
        else:
            basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
            self.taxon_namespaces = container.OrderedSet()
            self.tree_lists = container.OrderedSet()
            self.char_matrices = container.OrderedSet()
            self.attached_taxon_namespace = None
            self._process_taxon_namespace_directives(kwargs)
            self.comments = []
            if len(args) == 1 and not isinstance(args[0], DataSet):
                for item in args[0]:
                    self.add(item)
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def _clone_from(self, dataset, kwargs_dict):
        raise NotImplementedError

    def __copy__(self):
        raise NotImplementedError

    def taxon_namespace_scoped_copy(self, memo=None):
        raise NotImplementedError

    def __deepcopy__(self, memo=None):
        raise NotImplementedError

    ###########################################################################
    ### Data I/O

    def _parse_and_add_from_stream(self,
            stream,
            schema,
            exclude_trees=False,
            exclude_chars=False,
            **kwargs):
        # exclude_trees = kwargs.pop("exclude_trees", False)
        # exclude_chars = kwargs.pop("exclude_chars", False)
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        if (self.attached_taxon_namespace is not None
                and taxon_namespace is not None
                and self.attached_taxon_namespace is not taxon_namespace):
            raise ValueError("DataSet has attached TaxonNamespace that is not the same as ``taxon_namespace``")
        if self.attached_taxon_namespace is not None and taxon_namespace is None:
            taxon_namespace = self.attached_taxon_namespace
        label = kwargs.pop("label", None)
        reader = dataio.get_reader(schema, **kwargs)
        n_tns = len(self.taxon_namespaces)
        n_tree_lists = len(self.tree_lists)
        n_char_matrices = len(self.char_matrices)
        reader.read_dataset(
                stream=stream,
                dataset=self,
                taxon_namespace=taxon_namespace,
                exclude_trees=exclude_trees,
                exclude_chars=exclude_chars,
                state_alphabet_factory=charstatemodel.StateAlphabet,
                )
        n_tns2 = len(self.taxon_namespaces)
        n_tree_lists2 = len(self.tree_lists)
        n_char_matrices2 = len(self.char_matrices)
        return (n_tns2-n_tns,
                n_tree_lists2-n_tree_lists,
                n_char_matrices2-n_char_matrices)

    def read(self, **kwargs):
        """
        Add data to ``self`` from data source.

        **Mandatory Source-Specification Keyword Argument (Exactly One Required):**

            - **file** (*file*) -- File or file-like object of data opened for reading.
            - **path** (*str*) -- Path to file of data.
            - **url** (*str*) -- URL of data.
            - **data** (*str*) -- Data given directly.

        **Mandatory Schema-Specification Keyword Argument:**

            - **schema** (*str*) -- Identifier of format of data given by the
              "``file``", "``path``", "``data``", or "``url``" argument
              specified above: ":doc:`newick </schemas/newick>`", ":doc:`nexus
              </schemas/nexus>`", or ":doc:`nexml </schemas/nexml>`". See
              "|Schemas|" for more details.

        **Optional General Keyword Arguments:**

            - **exclude_trees** (*bool*) -- If |True|, then all tree data in the data
              source will be skipped.
            - **exclude_chars** (*bool*) -- If |True|, then all character
              data in the data source will be skipped.
            - **taxon_namespace** (|TaxonNamespace|) -- The |TaxonNamespace|
              instance to use to :doc:`manage the taxon names </primer/taxa>`.
              If not specified, a new one will be created unless the DataSet
              object is in attached taxon namespace mode
              (``self.attached_taxon_namespace`` is not |None| but assigned
              to a specific |TaxonNamespace| instance).
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

            ds = dendropy.DataSet()
            ds.read(
                    path="pythonidae.chars_and_trees.nex",
                    schema="nexus")
            ds.read(
                    url="http://purl.org/phylo/treebase/phylows/study/TB2:S1925?format=nexml",
                    schema="nexml")

        """
        return basemodel.MultiReadable._read_from(self, **kwargs)

    def _format_and_write_to_stream(self,
            stream,
            schema,
            exclude_trees=False,
            exclude_chars=False,
            **kwargs):
        """
        Writes out ``self`` in ``schema`` format to a destination given by
        file-like object ``stream``.

        Parameters
        ----------
        stream : file or file-like object
            Destination for data.
        schema : string
            Must be a recognized character file schema, such as "nexus",
            "phylip", etc, for which a specialized tree list writer is
            available. If this is not implemented for the schema specified, then
            a UnsupportedSchemaError is raised.

        \*\*kwargs : keyword arguments, optional
            Keyword arguments will be passed directly to the writer for the
            specified schema. See documentation for details on keyword
            arguments supported by writers of various schemas.

        """
        writer = dataio.get_writer(schema, **kwargs)
        writer.write_dataset(self, stream, exclude_trees, exclude_chars)

    ###########################################################################
    ### Domain Data Management

    ### General ###

    def add(self, data_object, **kwargs):
        """
        Generic add for TaxonNamespace, TreeList or CharacterMatrix objects.
        """
        if isinstance(data_object, taxonmodel.TaxonNamespace):
            self.add_taxon_namespace(data_object)
        elif isinstance(data_object, treecollectionmodel.TreeList):
            self.add_tree_list(data_object)
        elif isinstance(data_object, charmatrixmodel.CharacterMatrix):
            self.add_char_matrix(data_object)
        else:
            raise error.InvalidArgumentValueError("Cannot add object of type {} to DataSet" .format(type(data_object)))

    ### TaxonNamespace ###

    def add_taxon_namespace(self, taxon_namespace):
        """
        Adds a taxonomic unit concept namespace represented by a
        |TaxonNamespace| instance to this dataset if it is not already
        there.

        Parameters
        ----------
        taxon_namespace : |TaxonNamespace|
            The |TaxonNamespace| object to be added.
        """
        self.taxon_namespaces.add(taxon_namespace)
        return taxon_namespace

    def new_taxon_namespace(self, *args, **kwargs):
        """
        Creates a new |TaxonNamespace| object, according to the arguments given
        (passed to `TaxonNamespace()`), and adds it to this |DataSet|.
        """
        t = taxonmodel.TaxonNamespace(*args, **kwargs)
        self.add_taxon_namespace(t)
        return t

    def attach_taxon_namespace(self, taxon_namespace=None):
        """
        Forces all read() calls of this DataSet to use the same |TaxonNamespace|. If
        ``taxon_namespace`` If ``taxon_namespace`` is None, then a new |TaxonNamespace| will be
        created, added to ``self.taxon_namespaces``, and that is the |TaxonNamespace| that will be
        attached.
        """
        if taxon_namespace is None:
            raise TypeError("Automatic creation of a new TaxonNamespace is no longer supported: ``taxon_namespace`` argument required to be passed a valid 'TaxonNamespace' instance. E.g.,:\n\n    taxon_namespace = dendropy.TaxonNamespace()\n    dataset.attach_taxon_namespace(taxon_namespace)\n")
            taxon_namespace = self.new_taxon_namespace()
        if type(taxon_namespace) == type(True):
            raise ValueError("attach_taxon_namespace() no longer accepts a bool argument: a valid 'TaxonNamespace' instance is required")
        if taxon_namespace not in self.taxon_namespaces:
            self.add_taxon_namespace(taxon_namespace)
        self.attached_taxon_namespace = taxon_namespace
        return self.attached_taxon_namespace

    def detach_taxon_namespace(self):
        """
        Relaxes constraint forcing all
        """
        t = self.attached_taxon_namespace
        self.attached_taxon_namespace = None
        return t

    def _process_taxon_namespace_directives(self, kwargs_dict):
        """
        The following idioms are supported:

            `taxon_namespace=tns`
                Attach ``tns`` as the bound (single, unified) taxonomic namespace
                reference for all objects.
            `attached_taxon_namespace=tns`
                Attach ``tns`` as the bound (single, unified) taxonomic namespace
                reference for all objects.
            `attach_taxon_namespace=True, attached_taxon_namespace=tns`
                Attach ``tns`` as the bound (single, unified) taxonomic namespace
                reference for all objects.
            `attach_taxon_namespace=True`
                Create a *new* |TaxonNamespace| and set it as the bound
                (single, unified) taxonomic namespace reference for all
                objects.
        """
        attach_taxon_namespace, taxon_namespace = taxonmodel.process_attached_taxon_namespace_directives(kwargs_dict)
        if attach_taxon_namespace or (taxon_namespace is not None):
            self.attach_taxon_namespace(taxon_namespace)
        return attach_taxon_namespace, taxon_namespace

    def unify_taxon_namespaces(self,
            taxon_namespace=None,
            case_sensitive_label_mapping=True,
            attach_taxon_namespace=True):
        """
        Reindices taxa across all subcomponents, mapping to single taxon set.
        """
        if len(self.taxon_namespaces) or len(self.tree_lists) or len(self.char_matrices):
            self.taxon_namespaces.clear()
            if taxon_namespace is None:
                taxon_namespace = self.new_taxon_namespace()
            taxon_mapping_memo = {}
            for tree_list in self.tree_lists:
                tree_list.migrate_taxon_namespace(
                        taxon_namespace=taxon_namespace,
                        unify_taxa_by_label=True,
                        # case_sensitive_label_mapping=case_sensitive_label_mapping,
                        taxon_mapping_memo=taxon_mapping_memo)
            for char_matrix in self.char_matrices:
                char_matrix.migrate_taxon_namespace(
                        taxon_namespace=taxon_namespace,
                        unify_taxa_by_label=True,
                        # case_sensitive_label_mapping=case_sensitive_label_mapping,
                        taxon_mapping_memo=taxon_mapping_memo)
        if attach_taxon_namespace:
            self.attach_taxon_namespace(taxon_namespace)

    def unify_taxa(self, taxon_set=None, bind=None):
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: '{class_name}.unify_taxa()' will no longer be supported in future releases; use '{class_name}.unify_taxon_namespaces()' instead".format(class_name=self.__class__.__name__))
        self.unify_taxon_namespaces(taxon_namespace=taxon_set,
                attach_taxon_namespace=bind)

    ### **Legacy** ###

    def _get_taxon_sets(self):
        self.taxon_sets_deprecation_warning()
        return self.taxon_namespaces
    def _set_taxon_sets(self, v):
        self.taxon_sets_deprecation_warning()
        self.taxon_namespaces = v
    def _del_taxon_sets(self):
        self.taxon_sets_deprecation_warning()
    taxon_sets = property(_get_taxon_sets, _set_taxon_sets, _del_taxon_sets)

    def taxon_sets_deprecation_warning(self, stacklevel=4):
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'taxon_sets' will no longer be supported in future releases; use 'taxon_namespaces' instead",
                stacklevel=stacklevel)

    def _get_attached_taxon_set(self):
        self.attached_taxon_set_deprecation_warning()
        return self.attached_taxon_namespace
    def _set_attached_taxon_set(self, v):
        self.attached_taxon_set_deprecation_warning()
        self.attached_taxon_namespace = v
    def _del_attached_taxon_set(self):
        self.attached_taxon_set_deprecation_warning()
    attached_taxon_set = property(_get_attached_taxon_set, _set_attached_taxon_set, _del_attached_taxon_set)

    def attached_taxon_set_deprecation_warning(self, stacklevel=4):
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'attached_taxon_set' will no longer be supported in future releases; use 'attached_taxon_namespace' instead",
                stacklevel=3)

    def add_taxon_set(self, taxon_set):
        """
        DEPRECATED: Use `add_taxon_namespace()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'add_taxon_set' will no longer be supported in future releases; use 'add_taxon_namespace' instead",
                stacklevel=3)
        return self.add_taxon_namespace(taxon_namespace=taxon_set)

    def new_taxon_set(self, *args, **kwargs):
        """
        DEPRECATED: Use `new_taxon_namespace()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'new_taxon_set' will no longer be supported in future releases; use 'new_taxon_namespace' instead",
                stacklevel=3)
        return self.new_taxon_namespace(*args, **kwargs)

    def attach_taxon_set(self, taxon_set=None):
        """
        DEPRECATED: Use `attach_taxon_namespace()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'attach_taxon_set' will no longer be supported in future releases; use 'attach_taxon_namespace' instead",
                stacklevel=3)
        return self.attach_taxon_namespace(taxon_namespace=taxon_set)

    def detach_taxon_set(self):
        """
        DEPRECATED: Use `detach_taxon_namespace()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'detach_taxon_set' will no longer be supported in future releases; use 'detach_taxon_namespace' instead",
                stacklevel=3)
        self.detach_taxon_namespace()

    ### TreeList ###

    def add_tree_list(self, tree_list):
        """
        Adds a |TreeList| instance to this dataset if it is not already
        there.

        Parameters
        ----------
        tree_list : |TreeList|
            The |TreeList| object to be added.
        """
        if tree_list.taxon_namespace not in self.taxon_namespaces:
            self.taxon_namespaces.add(tree_list.taxon_namespace)
        self.tree_lists.add(tree_list)
        return tree_list

    def new_tree_list(self, *args, **kwargs):
        """
        Creates a new |TreeList| instance, adds it to this DataSet.

        Parameters
        ----------
        \*args : positional arguments
            Passed directly to |TreeList| constructor.
        \*\*kwargs : keyword arguments, optional
            Passed directly to |TreeList| constructor.

        Returns
        -------
        t : |TreeList|
            The new |TreeList| instance created.
        """
        if self.attached_taxon_namespace is not None:
            if "taxon_namespace" in kwargs and kwargs["taxon_namespace"] is not self.attached_taxon_namespace:
                raise TypeError("DataSet object is attached to TaxonNamespace {}, but 'taxon_namespace' argument specifies different TaxonNamespace {}" .format(
                    repr(self.attached_taxon_namespace), repr(kwargs["taxon_namespace"])))
            else:
                kwargs["taxon_namespace"] = self.attached_taxon_namespace
        tree_list = treecollectionmodel.TreeList(*args, **kwargs)
        return self.add_tree_list(tree_list)

    def get_tree_list(self, label):
        """
        Returns a TreeList object specified by label.
        """
        for t in self.tree_lists:
            if t.label == label:
                return t
        return None

    ### CharacterMatrix ###

    def add_char_matrix(self, char_matrix):
        """
        Adds a |CharacterMatrix| or |CharacterMatrix|-derived
        instance to this dataset if it is not already there.

        Parameters
        ----------
        char_matrix : |CharacterMatrix|
            The |CharacterMatrix| object to be added.
        """
        if char_matrix.taxon_namespace not in self.taxon_namespaces:
            self.taxon_namespaces.add(char_matrix.taxon_namespace)
        self.char_matrices.add(char_matrix)
        return char_matrix

    def new_char_matrix(self, char_matrix_type, *args, **kwargs):
        """
        Creation and accession of new |CharacterMatrix| (of class
        ``char_matrix_type``) into ``chars`` of self."
        """
        if self.attached_taxon_namespace is not None:
            if "taxon_namespace" in kwargs and kwargs["taxon_namespace"] is not self.attached_taxon_namespace:
                raise TypeError("DataSet object is attached to TaxonNamespace %s, but 'taxon_namespace' argument specifies different TaxonNamespace %s" % (
                    repr(self.attached_taxon_namespace), repr(kwargs["taxon_namespace"])))
            else:
                kwargs["taxon_namespace"] = self.attached_taxon_namespace
        if textprocessing.is_str_type(char_matrix_type):
            char_matrix = charmatrixmodel.new_char_matrix(
                    data_type=char_matrix_type,
                    *args,
                    **kwargs)
        else:
            char_matrix = char_matrix_type(*args, **kwargs)
        return self.add_char_matrix(char_matrix)
