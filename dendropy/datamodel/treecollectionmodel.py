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
This module handles the core definition of classes that model collections of
trees.
"""

import collections
import math
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import copy
import sys
from dendropy.utility import container
from dendropy.utility import error
from dendropy.utility import bitprocessing
from dendropy.utility import deprecate
from dendropy.utility import constants
from dendropy.mathlib import statistics
from dendropy.datamodel import basemodel
from dendropy.datamodel import taxonmodel
from dendropy.datamodel import treemodel
from dendropy import dataio

##############################################################################
### TreeList

class TreeList(
        taxonmodel.TaxonNamespaceAssociated,
        basemodel.Annotable,
        basemodel.Readable,
        basemodel.Writeable,
        basemodel.DataObject):
    """
    A collection of :class:`Tree` objects, all referencing the same "universe" of
    opeational taxonomic unit concepts through the same :class:`TaxonNamespace`
    object reference.
    """

    def _parse_from_stream(cls,
            stream,
            schema,
            collection_offset=None,
            tree_offset=None,
            **kwargs):
        """
        Constructs a new :class:`TreeList` object and populates it with trees from
        file-like object `stream`.

        Notes
        -----
        *All* operational taxonomic unit concepts in the data source will be included
        in the :class:`TaxonNamespace` object associated with the new
        :class:`TreeList` object and its contained :class:`Tree` objects, even those
        not associated with trees or the particular trees being retrieved.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`

        collection_offset : integer or None
            0-based index indicating collection of trees to parse. If `None`,
            then all tree collections are retrieved, with each distinct
            collection parsed into a separate :class:`TreeList` object. If the
            tree colleciton offset index is equal or greater than the number of
            tree collections in the data source, then IndexError is raised.
            Negative offsets work like negative list indexes; e.g., a
            `collection_offset` of -1 means to read the last collection of
            trees in the data source. For data formats that do not support the
            concept of distinct tree collections (e.g. NEWICK) are considered
            single-collection data source (i.e, the only acceptable
            `collection_offset` values are -1 or 0).

        tree_offset : integer or None
            0-based index indicating particular tree within a particular
            collection of trees at which to begin reading.  If not specified or
            `None` (default), then all trees are parsed.  Otherwise, must be an
            integer value up the length of the collection minus 1.  A positive
            offset indicates the number of trees in the collection to skip;
            e.g. a `tree_offset` of 20 means to skip the first 20 trees in the
            collection.  Negative offsets work like negative list indexes;
            e.g., a `tree_offset` value of -10 means to retrieve the last 10
            trees in the collection.  If the tree offset index is equal or
            greater than the number of trees in the collection, then IndexError
            is raised. Requires that a particular tree collection has been
            identified using the `tree_collection_offset` parameter: if
            `tree_collection_offset` is not specified, a `TypeError` is raised.

        \*\*kwargs : keyword arguments
            Arguments to customize parsing, instantiation, processing, and
            accession of :class:`Tree` objects read from the data source, including
            schema- or format-specific handling.

            The following optional keyword arguments are recognized and handled
            by this function:

                * `label` Specifies the label or description of the new
                  :class:`TreeList`.
                * `taxon_namespace` specifies the :class:`TaxonNamespace`
                   object to be attached to the new :class:`TreeList` object.
                   Note that *all* operational taxonomic unit concepts in the
                   data source will be accessioned into the specified
                   :class:`TaxonNamespace` instance. This includes the
                   operation taxonomic unit definitions associated with all
                   tree collections and character matrices in the data source.
                * `tree_list` : **SPECIAL** If passed a :class:`TreeList` using
                  this keyword, then this instance is populated and returned
                  (instead of a new instance being created).

            All other keyword arguments are passed directly to `TreeList.read()`.
            Other keyword arguments may be available, depending on the implementation
            of the reader specialized to handle `schema` formats.

        Notes
        -----
        Note that in most cases, even if `collection_offset` and `tree_offset`
        are specified to restrict the trees returned, the *entire* data source
        is still parsed and processed. So this is not more efficient than
        reading all the trees and then manually-extracting them later; just
        more convenient. If you need just a single subset of trees from a data
        source, there is no gain in efficiency. If you need multiple trees or
        subsets of trees from the same data source, it would be much more
        efficient to read the entire data source, and extract trees as needed.

        Returns
        -------
        A :class:`TreeList` object.

        """
        # these must be pulled before passing the kwargs
        # down to the reader
        tree_list = kwargs.pop("tree_list", None)
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        label = kwargs.pop("label", None)

        # get the reader
        reader = dataio.get_reader(schema, **kwargs)

        # Accommodate an existing TreeList object being passed
        if tree_list is None:
            tree_list = cls(label=label, taxon_namespace=taxon_namespace)

        if collection_offset is None:
            if tree_offset is not None:
                raise TypeError("Cannot specify `tree_offset` without specifying `collection_offset`")
            # coerce all tree products into this list
            reader.read_tree_lists(
                        stream=stream,
                        taxon_namespace_factory=tree_list._taxon_namespace_pseudofactory,
                        tree_list_factory=tree_list._tree_list_pseudofactory,
                        global_annotations_target=None)
        else:
            tree_lists = reader.read_tree_lists(
                        stream=stream,
                        taxon_namespace_factory=tree_list._taxon_namespace_pseudofactory,
                        tree_list_factory=tree_list.__class__,
                        global_annotations_target=None)
            # if collection_offset < 0:
            #     raise IndexError("Collection offset out of range: {} (minimum valid tree offset = 0)".format(collection_offset))
            if collection_offset >= len(tree_lists):
                raise IndexError("Collection offset out of range: {} (number of collections = {}, maximum valid collection offset = {})".format(collection_offset, len(tree_lists), len(tree_lists)-1))
            target_tree_list = tree_lists[collection_offset]
            tree_list.copy_annotations_from(target_tree_list)
            if tree_offset is not None:
                # if tree_offset < 0:
                #     raise IndexError("Tree offset out of range: {} (minimum offset = 0)".format(tree_offset))
                if tree_offset >= len(target_tree_list):
                    raise IndexError("Tree offset out of range: {} (number of trees in source = {}, maximum valid tree offset = {})".format(tree_offset, len(target_tree_list), len(target_tree_list)-1))
                for tree in target_tree_list[tree_offset:]:
                    tree_list._trees.append(tree)
            else:
                for tree in target_tree_list:
                    tree_list._trees.append(tree)
        return tree_list
        # taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        # label = kwargs.pop("label", None)
        # tree_list = cls(label=label,
        #         taxon_namespace=taxon_namespace)
        # tree_list.read(
        #         stream=stream,
        #         schema=schema,
        #         collection_offset=collection_offset,
        #         tree_offset=tree_offset,
        #         **kwargs)
        # return tree_list
    _parse_from_stream = classmethod(_parse_from_stream)

    DEFAULT_TREE_TYPE = treemodel.Tree

    def tree_factory(cls, *args, **kwargs):
        """
        Creates and returns a :class:`Tree` of a type that this list understands how to
        manage.

        Deriving classes can override this to provide for custom Tree-type
        object lists. You can simple override the class-level variable
        `DEFAULT_TREE_TYPE` in your derived class if the constructor signature
        of the alternate tree type is the same as `:class:Tree`.
        If you want to have a TreeList *instance* that generates
        custom trees (i.e., as opposed to a TreeList-ish *class* of instances),
        set the `tree_type` attribute of the TreeList instance.

        Parameters
        ----------
        \*args : positional arguments
            Passed directly to constructor of :class:`Tree`.

        \*\*kwargs : keyword arguments
            Passed directly to constructor of :class:`Tree`.

        Returns
        -------
        A :class:`Tree` object.

        """
        tree = cls.DEFAULT_TREE_TYPE(*args, **kwargs)
        return tree
    tree_factory = classmethod(tree_factory)

    ###########################################################################
    ### Lifecycle and Identity

    def __init__(self, *args, **kwargs):
        """
        Constructs a new :class:`TreeList` object, populating it with any iterable
        container with Tree object members passed as unnamed argument, or from
        a data source if `stream` and `schema` are passed.

        If passed an iterable container, the objects in that container must be
        of type :class:`Tree` (or derived). If the container is of type :class:`TreeList`,
        then, because each :class:`Tree` object must have the same :class:`TaxonNamespace`
        reference as the containing :class:`TreeList`, the trees in the container
        passed as an initialization argument will be **deep**-copied (except
        for associated :class:`TaxonNamespace` and :class:`Taxon` objects, which will
        be shallow-copied). If the container is any other type of
        iterable, then the :class:`Tree` objects will be **shallow**-copied.

        :class:`TreeList` objects can directly thus be instantiated in the
        following ways::

            # /usr/bin/env python

            import StringIO
            from dendropy import TaxonNamespace, Tree, TreeList

            # instantiate an empty tree
            tlst1 = TreeList()

            # the canonical way to instantiate a TreeList from a data source
            # is `get_from_*` family of static factory methods
            tlst2 = TreeList.get_from_stream(open('treefile.tre', 'rU'), "newick")
            tlst3 = TreeList.get_from_path('sometrees.nexus', "nexus")
            tlst4 = TreeList.get_from_string("((A,B),(C,D));((A,C),(B,D));", "newick")

            # can also call `read()` on a TreeList object; each read adds
            # (appends) the tree(s) found to the TreeList
            tlst5 = TreeList()
            tlst5.read(open('boot1.tre', 'rU'), "newick")
            tlst5.read_from_stream(open('boot2.tre', 'rU'), "newick") # same as above
            tlst5.read_from_string("((A,B),(C,D));((A,C),(B,D));", "newick")
            tlst5.read_from_path("boot3.tre", "newick")

            # populated from list of Tree objects
            tlist6_1 = Tree.get_from_string("((A,B),(C,D))",
                    schema="newick")
            tlist6_2 = Tree.get_from_string("((A,C),(B,D))",
                    schema="newick")
            tlist6 = TreeList([tlist5_1, tlist5_2])

            # passing keywords to underlying tree parser
            tlst8 = TreeList.get_from_string("((A,B),(C,D));((A,C),(B,D));",
                             schema="newick",
                             taxon_namespace=tlst3.taxon_namespace,
                             rooting="force-rooted",
                             extract_comment_metadata=True,
                             store_tree_weights=False,
                             preserve_underscores=True)

            # Subsets of trees can be read. Note that in most cases, the entire
            # data source is parsed, so this is not more efficient than reading
            # all the trees and then manually-extracting them later; just more
            # convenient

            # skip the first 100 trees in the first collection of trees
            trees = TreeList.get_from_path("mcmc.tre", "newick",
                        collection_offset=0, tree_offset=100)

            # get the last 10 trees in the first collection of trees
            trees = TreeList.get_from_path("mcmc.tre", "newick",
                        collection_offset=0, tree_offset=-10)

            # get the last 10 trees in the second-to-last collection of trees
            trees = TreeList.get_from_path("mcmc.xml", "nexml",
                        collection_offset=-2, tree_offset=-10)

            # Slices give shallow-copy: trees are references
            tlst4copy0a = t4[:]
            assert tlst4copy0a[0] is t4[0]
            tlst4copy0b = t4[:4]
            assert tlst4copy0b[0] is t4[0]

            # 'Taxon-namespace-scoped' copy:
            # I.e., Deep-copied objects but taxa and taxon namespace
            # are copied as references
            tlst4copy1a = TreeList(t4)
            tlst4copy1b = TreeList([Tree(t) for t in tlst5])
            assert tlst4copy1a[0] is not tlst4[0] # True
            assert tlst4copy1a.taxon_namespace is tlst4.taxon_namespace # True
            assert tlst4copy1b[0] is not tlst4[0] # True
            assert tlst4copy1b.taxon_namespace is tlst4.taxon_namespace # True


        """
        if len(args) > 1:
            # only allow 1 positional argument
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        elif len(args) == 1 and isinstance(args[0], TreeList):
            self._clone_from(args[0], kwargs)
        else:
            basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
            taxonmodel.TaxonNamespaceAssociated.__init__(self,
                    taxon_namespace=taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None))
            self.tree_type = kwargs.pop("tree_type", self.__class__.DEFAULT_TREE_TYPE)
            self._trees = []
            self.comments = []
            if len(args) == 1:
                for aidx, a in enumerate(args[0]):
                    if not isinstance(a, self.tree_type):
                        raise ValueError("Cannot add object not of 'Tree' type to 'TreeList'")
                    self.append(a)
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return (
            isinstance(other, TreeList)
            and (self.taxon_namespace is other.taxon_namespace)
            and (self._trees == other._trees)
        )

    def _clone_from(self, tree_list, kwargs_dict):
        memo = {}
        # memo[id(tree)] = self
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs_dict, tree_list.taxon_namespace)
        memo[id(tree_list.taxon_namespace)] = taxon_namespace
        if taxon_namespace is not tree_list.taxon_namespace:
            for t1 in tree_list.taxon_namespace:
                t2 = taxon_namespace.require_taxon(label=t1.label)
                memo[id(t1)] = t2
        else:
            for t1 in tree_list.taxon_namespace:
                memo[id(t1)] = t1
        t = copy.deepcopy(tree_list, memo)
        self.__dict__ = t.__dict__
        self.label = kwargs_dict.pop("label", tree_list.label)
        return self

    def __copy__(self):
        other = TreeList(label=self.label, taxon_namespace=self.taxon_namespace)
        other._trees = list(self._trees)
        memo = {}
        memo[id(self)] = other
        other.deep_copy_annotations_from(self, memo)
        return other

    def taxon_namespace_scoped_copy(self, memo=None):
        if memo is None:
            memo = {}
        # this populates `memo` with references to the
        # the TaxonNamespace and Taxon objects
        self.taxon_namespace.populate_memo_for_taxon_namespace_scoped_copy(memo)
        return self.__deepcopy__(memo=memo)

    def __deepcopy__(self, memo=None):
        return basemodel.Annotable.__deepcopy__(self, memo=memo)

    ###########################################################################
    ### Representation

    def __str__(self):
        return "<TreeList {} '{}': [{}]>".format(hex(id(self)), self.label, ", ".join(repr(i) for i in self._trees))

    ###########################################################################
    ### Data I/O

    def _taxon_namespace_pseudofactory(self, **kwargs):
        """
        Dummy factory to coerce all :class:`TaxonNamespace` objects required when
        parsing a data source to reference `self.taxon_namespace`.
        """
        if "label" in kwargs and kwargs["label"] is not None and self.taxon_namespace.label is None:
            self.taxon_namespace.label = kwargs["label"]
        return self.taxon_namespace

    def _tree_list_pseudofactory(self, **kwargs):
        """
        Dummy factory to coerce all :class:`TreeList` objects required when
        parsing a data source to reference `self`.
        """
        if "label" in kwargs and kwargs["label"] is not None and self.label is None:
            self.label = kwargs["label"]
        return self

    def read(self,
            stream,
            schema,
            collection_offset=None,
            tree_offset=None,
            **kwargs):
        """
        Parses :class:`Tree` objects from data source and adds to this collection.

        Notes
        -----
        *All* operational taxonomic unit concepts in the data source will be included
        in the :class:`TaxonNamespace` object associated with the new
        :class:`TreeList` object and its contained :class:`Tree` objects, even those
        not associated with trees or the particular trees being retrieved.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`.

        collection_offset : integer or None
            0-based index indicating collection of trees to parse. If `None`,
            then all tree collections are retrieved, with each distinct
            collection parsed into a separate :class:`TreeList` object. If the
            tree colleciton offset index is equal or greater than the number of
            tree collections in the data source, then IndexError is raised.
            Negative offsets work like negative list indexes; e.g., a
            `collection_offset` of -1 means to read the last collection of
            trees in the data source. For data formats that do not support the
            concept of distinct tree collections (e.g. NEWICK) are considered
            single-collection data source (i.e, the only acceptable
            `collection_offset` values are -1 or 0).

        tree_offset : integer or None
            0-based index indicating particular tree within a particular
            collection of trees at which to begin reading.  If not specified or
            `None` (default), then all trees are parsed.  Otherwise, must be an
            integer value up the length of the collection minus 1.  A positive
            offset indicates the number of trees in the collection to skip;
            e.g. a `tree_offset` of 20 means to skip the first 20 trees in the
            collection.  Negative offsets work like negative list indexes;
            e.g., a `tree_offset` value of -10 means to retrieve the last 10
            trees in the collection.  If the tree offset index is equal or
            greater than the number of trees in the collection, then IndexError
            is raised. Requires that a particular tree collection has been
            identified using the `tree_collection_offset` parameter: if
            `tree_collection_offset` is not specified, a `TypeError` is raised.

        \*\*kwargs : keyword arguments

            Arguments to customize parsing, instantiation, processing, and
            accession of :class:`Tree` objects read from the data source, including
            schema- or format-specific handling. These will be passed to the
            underlying schema-specific reader for handling.

            General (schema-agnostic) keyword arguments are:

                * `rooted` specifies the default rooting interpretation of the tree.
                * `edge_len_type` specifies the type of the edge lengths (int or
                  float; defaults to 'float')

            Other keyword arguments are available depending on the schema. See
            specific schema handlers (e.g., :class:`NewickReader`, :class:`NexusReader`,
            :class:`NexmlReader`) for more details.

        Notes
        -----
        Note that in most cases, even if `collection_offset` and `tree_offset`
        are specified to restrict the trees read, the *entire* data source
        is still parsed and processed. So this is not more efficient than
        reading all the trees and then manually-extracting them later; just
        more convenient. If you need just a single subset of trees from a data
        source, there is no gain in efficiency. If you need multiple trees or
        subsets of trees from the same data source, it would be much more
        efficient to read the entire data source, and extract trees as needed.

        Returns
        -------
        n : `int`
            The number of :class:`Tree` objects read.

        """
        if "taxon_namespace" in kwargs or "taxon_namespace" in kwargs:
            raise TypeError("Cannot change `taxon_namespace` when reading into an existing TreeList")
        kwargs["taxon_namespace"] = self.taxon_namespace
        kwargs["tree_list"] = self
        cur_size = len(self._trees)
        TreeList._parse_from_stream(
                stream=stream,
                schema=schema,
                collection_offset=collection_offset,
                tree_offset=tree_offset,
                **kwargs)
        new_size = len(self._trees)
        return new_size - cur_size

    def write(self, stream, schema, **kwargs):
        """
        Writes out `self` in `schema` format to a destination given by
        file-like object `stream`.

        Parameters
        ----------
        stream : file or file-like object
            Destination for data.
        schema : string
            Must be a recognized and tree file schema, such as "nexus",
            "newick", etc, for which a specialized tree list writer is
            available. If this is not implemented for the schema specified, then
            a UnsupportedSchemaError is raised.

        \*\*kwargs : keyword arguments, optional
            Keyword arguments will be passed directly to the writer for the
            specified schema. See documentation for details on keyword
            arguments supported by writers of various schemas.

        """
        writer = dataio.get_writer(schema, **kwargs)
        writer.write_tree_list(self, stream)

    ###########################################################################
    ### List Interface

    def _import_tree_to_taxon_namespace(self,
            tree,
            taxon_import_strategy="migrate",
            **kwargs):
        if tree.taxon_namespace is not self.taxon_namespace:
            if taxon_import_strategy == "migrate":
                tree.migrate_taxon_namespace(taxon_namespace=self.taxon_namespace,
                        **kwargs)
            elif taxon_import_strategy == "add":
                tree._taxon_namespace = self.taxon_namespace
                tree.update_taxon_namespace()
            else:
                raise ValueError("Unrecognized taxon import strategy: '{}'".format(taxon_import_strategy))
        # assert tree.taxon_namespace is self.taxon_namespace
        return tree

    def insert(self,
            index,
            tree,
            taxon_import_strategy="migrate",
            **kwargs):
        """
        Inserts a :class:`Tree` object, `tree`, into the collection before
        `index`.

        The :class:`TaxonNamespace` reference of `tree` will be set to that of
        `self`.  Any :class:`Taxon` objects associated with nodes in `tree`
        that are not already in `self.taxon_namespace` will be handled
        according to `taxon_import_strategy`:

            - 'migrate'
                :class:`Taxon` objects associated with `tree` that are not
                already in `self.taxon_nameaspace` will be remapped based on
                their labels, with new :class`Taxon` objects being
                reconstructed if none with matching labels are found.
                Specifically, '`tree.migrate_taxon_namespace(**kwargs)`' will
                be called on `tree`, where `kwargs` is as passed to this
                function.
            - 'add'
                :class:`Taxon` objects associated with `tree` that are not already
                in `self.taxon_namespace` will be added. Note that this might
                result in :class:`Taxon` objects with duplicate labels as no
                attempt at mapping to existing :class:`Taxon` objects based on
                label-matching is done.

        Parameters
        ----------
        index : integer
            Position before which to insert `tree`.
        tree : A :class:`Tree` instance
            The :class:`Tree` object to be added.
        taxon_import_strategy : string
            If `tree` is associated with a different :class:`TaxonNamespace`,
            this argument determines how new :class:`Taxon` objects in `tree`
            are handled: 'migrate' or 'add'. See above for details.
        \*\*kwargs : keyword arguments
            These arguments will be passed directly to
            'migrate_taxon_namespace()' method call on `tree`.

        See Also
        --------

        :meth:`Tree.migrate_taxon_namespace`

        """
        self._import_tree_to_taxon_namespace(
                tree=tree,
                taxon_import_strategy=taxon_import_strategy,
                **kwargs)
        self._trees.insert(index, tree)

    def append(self,
            tree,
            taxon_import_strategy="migrate",
            **kwargs):
        """
        Adds a :class:`Tree` object, `tree`, to the collection.

        The :class:`TaxonNamespace` reference of `tree` will be set to that of
        `self`.  Any :class:`Taxon` objects associated with nodes in `tree`
        that are not already in `self.taxon_namespace` will be handled
        according to `taxon_import_strategy`:

            - 'migrate'
                :class:`Taxon` objects associated with `tree` that are not
                already in `self.taxon_nameaspace` will be remapped based on
                their labels, with new :class`Taxon` objects being
                reconstructed if none with matching labels are found.
                Specifically, '`tree.migrate_taxon_namespace(**kwargs)`' will
                be called on `tree`, where `kwargs` is as passed to this
                function.
            - 'add'
                :class:`Taxon` objects associated with `tree` that are not already
                in `self.taxon_namespace` will be added. Note that this might
                result in :class:`Taxon` objects with duplicate labels as no
                attempt at mapping to existing :class:`Taxon` objects based on
                label-matching is done.

        Parameters
        ----------
        tree : A :class:`Tree` instance
            The :class:`Tree` object to be added.
        taxon_import_strategy : string
            If `tree` is associated with a different :class:`TaxonNamespace`,
            this argument determines how new :class:`Taxon` objects in `tree`
            are handled: 'migrate' or 'add'. See above for details.
        \*\*kwargs : keyword arguments
            These arguments will be passed directly to
            'migrate_taxon_namespace()' method call on `tree`.

        See Also
        --------

        :meth:`Tree.migrate_taxon_namespace`

        """
        self._import_tree_to_taxon_namespace(
                tree=tree,
                taxon_import_strategy=taxon_import_strategy,
                **kwargs)
        self._trees.append(tree)

    def extend(self, other):
        """
        In-place addition of :class:`Tree` objects in `other` to `self`.

        If `other` is a :class:`TreeList`, then the :class:`Trees` are *copied*
        and migrated into `self.taxon_namespace`; otherwise, the original
        objects are migrated into `self.taxon_namespace` and added directly.

        Parameters
        ----------
        other : iterable of :class:`Tree` objects

        Returns
        -------
        `self` : :class:`TreeList`
        """
        if isinstance(other, TreeList):
            for t0 in other:
                t1 = self.tree_type(t0, taxon_namespace=self.taxon_namespace)
                self._trees.append(t1)
        else:
            for t0 in other:
                self.append(t0)
        return self

    def __iadd__(self, other):
        """
        In-place addition of :class:`Tree` objects in `other` to `self`.

        If `other` is a :class:`TreeList`, then the :class:`Trees` are *copied*
        and migrated into `self.taxon_namespace`; otherwise, the original
        objects are migrated into `self.taxon_namespace` and added directly.

        Parameters
        ----------
        other : iterable of :class:`Tree` objects

        Returns
        -------
        `self` : :class:`TreeList`
        """
        return self.extend(other)

    def __add__(self, other):
        """
        Creates and returns new :class:`TreeList` with clones of all :class:`Trees` in `self`
        as well as all :class:`Tree` objects in `other`.  If `other` is a
        :class:`TreeList`, then the :class:`Trees` are *cloned* and migrated into
        `self.taxon_namespace`; otherwise, the original objects are migrated into
        `self.taxon_namespace` and added directly.

        Parameters
        ----------
        other : iterable of :class:`Tree` objects

        Returns
        -------
        tlist : :class:`TreeList` object
            :class:`TreeList` object containing clones of :class:`Tree` objects
            in `self` and `other`.
        """
        tlist = TreeList(taxon_namespace=self.taxon_namespace)
        tlist += self
        tlist += other
        return tlist

    def __contains__(self, tree):
        return tree in self._trees

    def __delitem__(self, tree):
        del self._trees[tree]

    def __iter__(self):
        return iter(self._trees)

    def __reversed__(self):
        return reversed(self._trees)

    def __len__(self):
        return len(self._trees)

    def __getitem__(self, index):
        """
        If `index` is an integer, then :class:`Tree` object at position `index`
        is returned. If `index` is a slice, then a :class:`TreeList` is returned
        with references (i.e., not copies or clones, but the actual original
        instances themselves) to :class:`Tree` objects in the positions given
        by the slice. The :class:`TaxonNamespace` is the same as `self`.

        Parameters
        ----------
        index : integer or slice
            Index or slice.

        Returns
        -------
        t : :class:`Tree` object or :class:`TreeList` object

        """
        if isinstance(index, slice):
            r = self._trees[index]
            return TreeList(r,
                    taxon_namespace=self.taxon_namespace)
        else:
            return self._trees[index]

    def __setitem__(self, index, value):
        if isinstance(index, slice):
            if isinstance(value, TreeList):
                tt = []
                for t0 in value:
                    t1 = self.tree_type(t0,
                            taxon_namespace=self.taxon_namespace)
                    tt.append(t1)
                value = tt
            else:
                for t in value:
                    self._import_tree_to_taxon_namespace(t)
            self._trees[index] = value
        else:
            self._trees[index] = self._import_tree_to_taxon_namespace(value)

    def clear(self):
        # list.clear() only with 3.4 or so ...
        self._trees = []

    def index(self, tree):
        return self._trees.index(tree)

    def pop(self, index=-1):
        return self._trees.pop(index)

    def remove(self, tree):
        self._trees.remove(tree)

    def reverse(self):
        self._trees.reverse()

    def sort(self, key=None, reverse=False):
        self._trees.sort(key=key, reverse=reverse)

    def new_tree(self, *args, **kwargs):
        tns = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, self.taxon_namespace)
        if tns is not self.taxon_namespace:
            raise TypeError("Cannot create new Tree with different TaxonNamespace")
        kwargs["taxon_namespace"] = self.taxon_namespace
        if self.tree_type is not None:
            tree = self.tree_type(*args, **kwargs)
        else:
            tree = self.tree_factory(*args, **kwargs)
        self._trees.append(tree)
        return tree

   ##############################################################################
   ## Taxon Handling

    def reconstruct_taxon_namespace(self,
            unify_taxa_by_label=True,
            taxon_mapping_memo=None):
        if taxon_mapping_memo is None:
            taxon_mapping_memo = {}
        for tree in self._trees:
            tree._taxon_namespace = self.taxon_namespace
            tree.reconstruct_taxon_namespace(
                unify_taxa_by_label=unify_taxa_by_label,
                taxon_mapping_memo=taxon_mapping_memo,
            )

    def update_taxon_namespace(self):
        for tree in self._trees:
            tree._taxon_namespace = self.taxon_namespace
            tree.update_taxon_namespace()

    def poll_taxa(self, taxa=None):
        """
        Returns a set populated with all of :class:`Taxon` instances associated
        with `self`.

        Parameters
        ----------
        taxa : set()
            Set to populate. If not specified, a new one will be created.

        Returns
        -------
        taxa : set[:class:`Taxon`]
            Set of taxa associated with `self`.
        """
        if taxa is None:
            taxa = set()
        for tree in self:
            tree.poll_taxa(taxa)
        return taxa

    def reindex_subcomponent_taxa():
        raise NotImplementedError()

   ##############################################################################
   ## Special Calculations and Operations on Entire Collection

    def consensus(self, min_freq=0.5, is_bipartitions_updated=False, **kwargs):
        """
        Returns a consensus tree of all trees in self, with minumum frequency
        of bipartition to be added to the consensus tree given by `min_freq`.
        """
        from dendropy.calculate import treesum
        return treesum.consensus_tree(
                self,
                min_freq=min_freq,
                is_bipartitions_updated=is_bipartitions_updated,
                **kwargs)

    def frequency_of_bipartition(self, **kwargs):
        """
        Given a bipartition specified as:

            - a :class:`Bipartition` instance given the keyword 'bipartition'
            - a split bitmask given the keyword 'split_bitmask'
            - a list of `Taxon` objects given with the keyword `taxa`
            - a list of taxon labels given with the keyword `labels`

        this function returns the proportion of trees in self
        in which the split is found.
        """
        split = None
        if "split_bitmask" in kwargs:
            split = kwargs["split_bitmask"]
        elif "bipartition" in kwargs:
            split = kwargs["bipartition"].split_bitmask
        elif "taxa" in kwargs or "labels" in kwargs:
            split = self.taxon_namespace.taxa_bitmask(**kwargs)
            if "taxa" in kwargs:
                k = len(kwargs["taxa"])
            else:
                k = len(kwargs["labels"])
            if bitprocessing.num_set_bits(split) != k:
                raise IndexError('Not all taxa could be mapped to bipartition (%s): %s' \
                    % (self.taxon_namespace.bitmask_as_bitstring(split), k))
        found = 0
        total = 0
        is_bipartitions_updated = kwargs.get("is_bipartitions_updated", False)
        for tree in self:
            if not is_bipartitions_updated or not tree.bipartitions:
                tree.encode_bipartitions()
            bipartition_encoding = set(b.split_bitmask for b in tree.bipartition_encoding)
            total += 1
            if split in bipartition_encoding:
                found += 1
        try:
            return float(found)/total
        except ZeroDivisionError:
            return 0

###############################################################################
### SplitDistribution

class SplitDistribution(taxonmodel.TaxonNamespaceAssociated):
    """
    Collects information regarding splits over multiple trees.
    """

    def __init__(self,
            taxon_namespace=None,
            ignore_edge_lengths=False,
            ignore_node_ages=True,
            use_tree_weights=True,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION):

        # Taxon Namespace
        taxonmodel.TaxonNamespaceAssociated.__init__(self,
                taxon_namespace=taxon_namespace)

        # configuration
        self.ignore_edge_lengths = ignore_edge_lengths
        self.ignore_node_ages = ignore_node_ages
        self.use_tree_weights = use_tree_weights
        self.ultrametricity_precision = ultrametricity_precision

        # storage
        self.total_trees_counted = 0
        self.sum_of_tree_weights = 0.0
        self.tree_rooting_types_counted = set()
        self.split_counts = collections.defaultdict(float)
        self.split_edge_lengths = collections.defaultdict(list)
        self.split_node_ages = collections.defaultdict(list)

        # secondary/derived/generated/collected data
        self._is_rooted = False
        self._split_freqs = None
        self._trees_counted_for_freqs = 0
        self._split_edge_length_summaries = None
        self._split_node_age_summaries = None
        self._trees_counted_for_summaries = 0

        # services
        self.tree_decorator = None

    ###########################################################################
    ### Utility

    def normalize_bitmask(self, bitmask):
        """
        "Normalizes" split, by ensuring that the least-significant bit is
        always 1 (used on unrooted trees to establish split identity
        independent of rotation).

        Parameters
        ----------
        bitmask : integer
            Split bitmask hash to be normalized.

        Returns
        -------
        h : integer
            Normalized split bitmask.
        """
        return treemodel.Bipartition.normalize_bitmask(
                bitmask=bitmask,
                fill_bitmask=self.taxon_namespace.all_taxa_bitmask(),
                lowest_relevant_bit=1)

    ###########################################################################
    ### Configuration

    def _is_rooted_deprecation_warning(self):
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'SplitDistribution.is_rooted' and 'SplitDistribution.is_unrooted' are no longer valid attributes; rooting state tracking and management is now the responsibility of client code.",
                stacklevel=4,
                )
    def _get_is_rooted(self):
        self._is_rooted_deprecation_warning()
        return self._is_rooted
    def _set_is_rooted(self, val):
        self._is_rooted_deprecation_warning()
        self._is_rooted = val
    is_rooted = property(_get_is_rooted, _set_is_rooted)
    def _get_is_unrooted(self):
        self._is_rooted_deprecation_warning()
        return not self._is_rooted
    def _set_is_unrooted(self, val):
        self._is_rooted_deprecation_warning()
        self._is_rooted = not val
    is_unrooted = property(_get_is_unrooted, _set_is_unrooted)

    ###########################################################################
    ### Split Counting and Book-Keeping

    def add_split_count(self, split, count=1):
        self.split_counts[split] += count

    def count_splits_on_tree(self,
            tree,
            is_bipartitions_updated=False,
            default_edge_length_value=None):
        """
        Counts splits in this tree and add to totals. `tree` must be decorated
        with splits, and no attempt is made to normalize taxa.

        Parameters
        ----------
        tree : a :class:`Tree` object.
            The tree on which to count the splits.
        is_bipartitions_updated : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.

        Returns
        --------
        s : iterable of splits
            A list of split bitmasks from `tree`.
        e :
            A list of edge length values from `tree`.
        a :
            A list of node age values from `tree`.
        """
        assert tree.taxon_namespace is self.taxon_namespace
        self.total_trees_counted += 1
        if not self.ignore_node_ages:
            tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
        if tree.weight is not None and self.use_tree_weights:
            weight_to_use = float(tree.weight)
        else:
            weight_to_use = 1.0
        self.sum_of_tree_weights += weight_to_use
        if tree.is_rooted:
            self.tree_rooting_types_counted.add(True)
        else:
            self.tree_rooting_types_counted.add(False)
        if not is_bipartitions_updated:
            tree.encode_bipartitions()
        splits = []
        edge_lengths = []
        node_ages = []
        for bipartition in tree.bipartition_encoding:
            split = bipartition.split_bitmask

            ## if edge is stored as an attribute, might be faster to:
            # edge = bipartition.edge
            edge = tree.bipartition_edge_map[bipartition]

            splits.append(split)
            self.split_counts[split] += weight_to_use
            if not self.ignore_edge_lengths:
                sel = self.split_edge_lengths.setdefault(split,[])
                if edge.length is None:
                    elen = default_edge_length_value
                else:
                    elen = edge.length
                sel.append(elen)
                edge_lengths.append(elen)
            else:
                sel = None
            if not self.ignore_node_ages:
                sna = self.split_node_ages.setdefault(split, [])
                if edge.head_node is not None:
                    nage = edge.head_node.age
                else:
                    nage = None
                sna.append(nage)
                node_ages.append(nage)
            else:
                sna = None
        return splits, edge_lengths, node_ages

    def splits_considered(self):
        """
        Returns 4 values:
            total number of splits counted
            total *weighted* number of unique splits counted
            total number of non-trivial splits counted
            total *weighted* number of unique non-trivial splits counted
        """
        if not self.split_counts:
            return 0, 0, 0, 0
        num_splits = 0
        num_unique_splits = 0
        num_nt_splits = 0
        num_nt_unique_splits = 0
        taxa_mask = self.taxon_namespace.all_taxa_bitmask()
        for s in self.split_counts:
            num_unique_splits += 1
            num_splits += self.split_counts[s]
            if not treemodel.Bipartition.is_trivial_bitmask(s, taxa_mask):
                num_nt_unique_splits += 1
                num_nt_splits += self.split_counts[s]
        return num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits

    def calc_freqs(self):
        "Forces recalculation of frequencies."
        self._split_freqs = {}
        if self.total_trees_counted == 0:
            for split in self.split_counts:
                self._split_freqs[split] = 1.0
        else:
            # total = self.total_trees_counted
            if not self.sum_of_tree_weights:
                total_weight = self.total_trees_counted
            else:
                total_weight = float(self.sum_of_tree_weights)
            for split in self.split_counts:
                count = self.split_counts[split]
                self._split_freqs[split] = float(self.split_counts[split]) / total_weight
        self._trees_counted_for_freqs = self.total_trees_counted
        self._split_edge_length_summaries = None
        self._split_node_age_summaries = None
        return self._split_freqs

    def update(self, split_dist):
        self.total_trees_counted += split_dist.total_trees_counted
        self.sum_of_tree_weights += split_dist.sum_of_tree_weights
        self._split_edge_length_summaries = None
        self._split_node_age_summaries = None
        self._trees_counted_for_summaries = 0
        self.tree_rooting_types_counted.update(split_dist.tree_rooting_types_counted)
        for split in split_dist.split_counts:
            self.split_counts[split] += split_dist.split_counts[split]
            self.split_edge_lengths[split] += split_dist.split_edge_lengths[split]
            self.split_node_ages[split] += split_dist.split_node_ages[split]

    ###########################################################################
    ### Basic Information Access

    def __len__(self):
        return len(self.split_counts)

    def __getitem__(self, split_bitmask):
        """
        Returns freqency of split_bitmask.
        """
        return self._get_split_frequencies().get(split_bitmask, 0.0)

    def _get_split_frequencies(self):
        if self._split_freqs is None or self._trees_counted_for_freqs != self.total_trees_counted:
            self.calc_freqs()
        return self._split_freqs
    split_frequencies = property(_get_split_frequencies)

    def is_mixed_rootings_counted(self):
        return ( (True in self.tree_rooting_types_counted)
                and (False in self.tree_rooting_types_counted or None in self.tree_rooting_types_counted) )

    def is_all_counted_trees_rooted(self):
        return (True in self.tree_rooting_types_counted) and (len(self.tree_rooting_types_counted) == 1)

    def is_all_counted_trees_strictly_unrooted(self):
        return (False in self.tree_rooting_types_counted) and (len(self.tree_rooting_types_counted) == 1)

    def is_all_counted_trees_treated_as_unrooted(self):
        return True not in self.tree_rooting_types_counted

    ###########################################################################
    ### Summarization

    def split_support_iter(self,
            tree,
            is_bipartitions_updated=False,
            include_external_splits=False,
            traversal_strategy="preorder",
            node_support_attr_name=None,
            edge_support_attr_name=None,
            ):
        """
        Returns iterator over support values for the splits of a given tree,
        where the support value is given by the proportional frequency of the
        split in the current split distribution.

        Parameters
        ----------
        tree : :class:`Tree`
            The :class:`Tree` which will be scored.
        is_bipartitions_updated : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.
        include_external_splits : bool
            If `True`, then non-internal split posteriors will be included.
            If `False`, then these are skipped. This should only make a
            difference when dealing with splits collected from trees of
            different leaf sets.
        traversal_strategy : str
            One of: "preorder" or "postorder". Specfies order in which splits
            are visited.

        Returns
        -------
        s : list of floats
            List of values for splits in the tree corresponding to the
            proportional frequency that the split is found in the current
            distribution.
        """
        if traversal_strategy == "preorder":
            if include_external_splits:
                iter_func = tree.preorder_node_iter
            else:
                iter_func = tree.preorder_internal_node_iter
        elif traversal_strategy == "postorder":
            if include_external_splits:
                iter_func = tree.postorder_node_iter
            else:
                iter_func = tree.postorder_internal_node_iter
        else:
            raise ValueError("Traversal strategy not supported: '{}'".format(traversal_strategy))
        if not is_bipartitions_updated:
            tree.encode_bipartitions()
        split_frequencies = self._get_split_frequencies()
        for nd in iter_func():
            split = nd.edge.split_bitmask
            support = split_frequencies.get(split, 0.0)
            yield support

    def summarize_edge_lengths(self):
        self._split_edge_length_summaries = {}
        for split, elens in self.split_edge_lengths.items():
            if not elens:
                continue
            try:
                self._split_edge_length_summaries[split] = statistics.summarize(elens)
            except ValueError:
                pass
        return self._split_edge_length_summaries

    def summarize_node_ages(self):
        self._split_node_age_summaries = {}
        for split, ages in self.split_node_ages.items():
            if not ages:
                continue
            try:
                self._split_node_age_summaries[split] = statistics.summarize(ages)
            except ValueError:
                pass
        return self._split_node_age_summaries

    def _get_split_edge_length_summaries(self):
        if self._split_edge_length_summaries is None \
                or self._trees_counted_for_summaries != self.total_trees_counted:
            self.summarize_edge_lengths()
        return dict(self._split_edge_length_summaries)
    split_edge_length_summaries = property(_get_split_edge_length_summaries)

    def _get_split_node_age_summaries(self):
        if self._split_node_age_summaries is None \
                or self._trees_counted_for_summaries != self.total_trees_counted:
            self.summarize_node_ages()
        return dict(self._split_node_age_summaries)
    split_node_age_summaries = property(_get_split_node_age_summaries)

    def log_product_of_split_support_on_tree(self,
            tree,
            is_bipartitions_updated=False,
            include_external_splits=False,
            ):
        """
        Calculates the (log) product of the support of the splits of the
        tree, where the support is given by the proportional frequency of the
        split in the current split distribution.

        The tree that has the highest product of split support out of a sample
        of trees corresponds to the "maximum credibility tree" for that sample.
        This can also be referred to as the "maximum clade credibility tree",
        though this latter term is sometimes use for the tree that has the
        highest *sum* of split support (see
        :meth:`SplitDistribution.sum_of_split_support_on_tree()`).

        Parameters
        ----------
        tree : :class:`Tree`
            The tree for which the score should be calculated.
        is_bipartitions_updated : bool
            If `True`, then the splits are assumed to have already been encoded
            and will not be updated on the trees.
        include_external_splits : bool
            If `True`, then non-internal split posteriors will be included in
            the score. Defaults to `False`: these are skipped. This should only
            make a difference when dealing with splits collected from trees of
            different leaf sets.

        Returns
        -------
        s : numeric
            The log product of the support of the splits of the tree.
        """
        log_product_of_split_support = 0.0
        for split_support in self.split_support_iter(
                tree=tree,
                is_bipartitions_updated=is_bipartitions_updated,
                include_external_splits=include_external_splits,
                traversal_strategy="preorder",
                ):
            if split_support:
                log_product_of_split_support += math.log(split_support)
        return log_product_of_split_support

    def sum_of_split_support_on_tree(self,
            tree,
            is_bipartitions_updated=False,
            include_external_splits=False,
            ):
        """
        Calculates the sum of the support of the splits of the tree, where the
        support is given by the proportional frequency of the split in the
        current distribtion.

        Parameters
        ----------
        tree : :class:`Tree`
            The tree for which the score should be calculated.
        is_bipartitions_updated : bool
            If `True`, then the splits are assumed to have already been encoded
            and will not be updated on the trees.
        include_external_splits : bool
            If `True`, then non-internal split posteriors will be included in
            the score. Defaults to `False`: these are skipped. This should only
            make a difference when dealing with splits collected from trees of
            different leaf sets.

        Returns
        -------
        s : numeric
            The sum of the support of the splits of the tree.
        """
        sum_of_split_support = 0.0
        for split_support in self.split_support_iter(
                tree=tree,
                is_bipartitions_updated=is_bipartitions_updated,
                include_external_splits=include_external_splits,
                traversal_strategy="preorder",
                ):
            sum_of_split_support += split_support
        return sum_of_split_support

    def consensus_tree(self,
            min_freq=constants.GREATER_THAN_HALF,
            is_rooted=None,
            **split_summarization_kwargs
            ):
        """
        Returns a consensus tree from splits in `self`.

        Parameters
        ----------

        min_freq : real
            The minimum frequency of a split in this distribution for it to be
            added to the tree.

        is_rooted : bool
            Should tree be rooted or not? If *all* trees counted for splits are
            explicitly rooted or unrooted, then this will default to `True` or
            `False`, respectively. Otherwise it defaults to `None`.

        \*\*split_summarization_kwargs : keyword arguments
            These will be passed directly to the underlying
            :class:`SplitDistributionSummarizer` object. See
            :meth:`SplitDistributionSummarizer.configure` for options.

        Returns
        -------
        t : consensus tree

        """
        if is_rooted is None:
            if self.is_all_counted_trees_rooted():
                is_rooted = True
            elif self.is_all_counted_trees_strictly_unrooted:
                is_rooted = False
        split_frequencies = self._get_split_frequencies()
        to_try_to_add = []
        _almost_one = lambda x: abs(x - 1.0) <= 0.0000001
        for s in split_frequencies:
            freq = split_frequencies[s]
            if (min_freq is None) or (freq >= min_freq) or (_almost_one(min_freq) and _almost_one(freq)):
                to_try_to_add.append((freq, s))
        to_try_to_add.sort(reverse=True)
        splits_for_tree = [i[1] for i in to_try_to_add]
        con_tree = treemodel.Tree.from_split_bitmasks(
                split_bitmasks=splits_for_tree,
                taxon_namespace=self.taxon_namespace,
                is_rooted=is_rooted)
        self.summarize_splits_on_tree(
            tree=con_tree,
            is_bipartitions_updated=False,
            **split_summarization_kwargs
            )
        return con_tree

    def summarize_splits_on_tree(self,
            tree,
            is_bipartitions_updated=False,
            **split_summarization_kwargs
            ):
        """
        Summarizes support of splits/edges/node on tree.

        Parameters
        ----------

        tree: :class:`Tree` instance
            Tree to be decorated with support values.

        is_bipartitions_updated: bool
            If `True`, then bipartitions will not be recalculated.

        \*\*split_summarization_kwargs : keyword arguments
            These will be passed directly to the underlying
            :class:`SplitDistributionSummarizer` object. See
            :meth:`SplitDistributionSummarizer.configure` for options.

        """
        if self.tree_decorator is None:
            self.tree_decorator = SplitDistributionSummarizer()
        self.tree_decorator.configure(**split_summarization_kwargs)
        self.tree_decorator.summarize_splits_on_tree(
                split_distribution=self,
                tree=tree,
                is_bipartitions_updated=is_bipartitions_updated)
        return tree

    ###########################################################################
    ### legacy

    def _get_taxon_set(self):
        from dendropy import taxonmodel
        taxon_model.taxon_set_deprecation_warning()
        return self.taxon_namespace

    def _set_taxon_set(self, v):
        from dendropy import taxonmodel
        taxon_model.taxon_set_deprecation_warning()
        self.taxon_namespace = v

    def _del_taxon_set(self):
        from dendropy import taxonmodel
        taxon_model.taxon_set_deprecation_warning()

    taxon_set = property(_get_taxon_set, _set_taxon_set, _del_taxon_set)

###############################################################################
### SplitDistributionSummarizer

class SplitDistributionSummarizer(object):

    def __init__(self, **kwargs):
        """
        See :meth:`SplitDistributionSummarizer.configure` for configuration
        options.
        """
        self.configure(**kwargs)

    def configure(self, **kwargs):
        """
        Configure rendition/mark-up.

        Parameters
        ----------

        set_edge_lengths : string
            For each edge, set the length based on:

                - "support": use support values split corresponding to edge
                - "mean-length": mean of edge lengths for split
                - "median-length": median of edge lengths for split
                - "mean-age": such that split age is equal to mean of ages
                - "median-age": such that split age is equal to mean of ages
                - `None`: do not set edge lengths

        add_support_as_node_attribute: bool
            Adds each node's support value as an attribute of the node,
            "`support`".

        add_support_as_node_annotation: bool
            Adds support as a metadata annotation, "`support". If
            `add_support_as_node_attribute` is `True`, then the value will be
            dynamically-bound to the value of the node's "`support`" attribute.

        set_support_as_node_label : bool
            Sets the `label` attribute of each node to the support value.

        add_node_age_summaries_as_node_attributes: bool
            Summarizes the distribution of the ages of each node in the
            following attributes:

                - `age_mean`
                - `age_median`
                - `age_sd`
                - `age_hpd95`
                - `age_range`

        add_node_age_summaries_as_node_annotations: bool
            Summarizes the distribution of the ages of each node in the
            following metadata annotations:

                - `age_mean`
                - `age_median`
                - `age_sd`
                - `age_hpd95`
                - `age_range`

            If `add_node_age_summaries_as_node_attributes` is `True`, then the
            values will be dynamically-bound to the corresponding node
            attributes.

        add_edge_length_summaries_as_edge_attributes: bool
            Summarizes the distribution of the lengths of each edge in the
            following attribtutes:

                - `length_mean`
                - `length_median`
                - `length_sd`
                - `length_hpd95`
                - `length_range`

        add_edge_length_summaries_as_edge_annotations: bool
            Summarizes the distribution of the lengths of each edge in the
            following metadata annotations:

                - `length_mean`
                - `length_median`
                - `length_sd`
                - `length_hpd95`
                - `length_range`

            If `add_edge_length_summaries_as_edge_attributes` is `True`, then the
            values will be dynamically-bound to the corresponding edge
            attributes.

        support_label_decimals: int
            Number of decimal places to express when rendering the support
            value as a string for the node label.

        support_as_percentages: bool
            Whether or not to express the support value as percentages (default
            is probability or proportion).

        minimum_edge_length : numeric
            All edge lengths calculated to have a value less than this will be
            set to this.

        error_on_negative_edge_lengths : bool
            If `True`, an inferred edge length that is less than 0 will result
            in a `ValueError`.

        """
        self.set_edge_lengths = kwargs.pop("set_edge_lengths", None)
        self.add_support_as_node_attribute = kwargs.pop("add_support_as_node_attribute", True)
        self.add_support_as_node_annotation = kwargs.pop("add_support_as_node_annotation", True)
        self.set_support_as_node_label = kwargs.pop("set_support_as_node_label", None)
        self.add_node_age_summaries_as_node_attributes = kwargs.pop("add_node_age_summaries_as_node_attributes", True)
        self.add_node_age_summaries_as_node_annotations = kwargs.pop("add_node_age_summaries_as_node_annotations", True)
        self.add_edge_length_summaries_as_edge_attributes = kwargs.pop("add_edge_length_summaries_as_edge_attributes", True)
        self.add_edge_length_summaries_as_edge_annotations = kwargs.pop("add_edge_length_summaries_as_edge_annotations", True)
        self.support_label_decimals = kwargs.pop("support_label_decimals", 4)
        self.support_as_percentages = kwargs.pop("support_as_percentages", False)
        self.support_label_compose_func = kwargs.pop("support_label_compose_func", None)
        self.primary_fieldnames = ["support",]
        self.summary_stats_fieldnames = ['mean', 'median', 'sd', 'hpd95', 'quant_5_95', 'range']
        self.node_age_summaries_fieldnames = list("age_{}".format(f) for f in self.summary_stats_fieldnames)
        self.edge_length_summaries_fieldnames = list("length_{}".format(f) for f in self.summary_stats_fieldnames)
        self.fieldnames = self.primary_fieldnames + self.node_age_summaries_fieldnames + self.edge_length_summaries_fieldnames
        for fieldname in self.fieldnames:
            setattr(self, "{}_attr_name".format(fieldname), kwargs.pop("{}_attr_name".format(fieldname), fieldname))
            setattr(self, "{}_annotation_name".format(fieldname), kwargs.pop("{}_annotation_name".format(fieldname), fieldname))
            setattr(self, "is_{}_annotation_dynamic".format(fieldname), kwargs.pop("is_{}_annotation_dynamic".format(fieldname), True))
        self.minimum_edge_length = kwargs.pop("minimum_edge_length", None)
        self.error_on_negative_edge_lengths = kwargs.pop("error_on_negative_edge_lengths", False)
        if kwargs:
            TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def _decorate(self,
            target,
            fieldname,
            value,
            set_attribute,
            set_annotation,
            ):
        attr_name = getattr(self, "{}_attr_name".format(fieldname))
        annotation_name = getattr(self, "{}_annotation_name".format(fieldname))
        if set_attribute:
            setattr(target, attr_name, value)
            if set_annotation:
                target.annotations.drop(name=annotation_name)
                if getattr(self, "is_{}_annotation_dynamic".format(fieldname)):
                    target.annotations.add_bound_attribute(
                        attr_name=attr_name,
                        annotation_name=annotation_name,
                        )
                else:
                    target.annotations.add_new(
                            name=annotation_name,
                            value=value,
                            )
        elif set_annotation:
            target.annotations.drop(name=annotation_name)
            target.annotations.add_new(
                    name=annotation_name,
                    value=value,
                    )

    def summarize_splits_on_tree(self,
            split_distribution,
            tree,
            is_bipartitions_updated=False):
        if not is_bipartitions_updated:
            tree.encode_bipartitions()
        if self.support_label_compose_func is not None:
            support_label_fn = lambda freq: self.support_label_compose_func(freq)
        else:
            support_label_fn = lambda freq: "{:.{places}f}".format(freq, places=self.support_label_decimals)
        node_age_summaries = split_distribution.split_node_age_summaries
        edge_length_summaries = split_distribution.split_edge_length_summaries
        split_freqs = split_distribution.split_frequencies
        assert len(self.node_age_summaries_fieldnames) == len(self.summary_stats_fieldnames)
        for node in tree:
            split_bitmask = node.edge.bipartition.split_bitmask
            split_support = split_freqs.get(split_bitmask, 0.0)
            if self.support_as_percentages:
                split_support = split_support * 100
            self._decorate(
                target=node,
                fieldname="support",
                value=split_support,
                set_attribute=self.add_support_as_node_attribute,
                set_annotation=self.add_support_as_node_annotation,
                )
            if self.set_support_as_node_label:
                node.label = support_label_fn(split_support)
            if (self.add_node_age_summaries_as_node_attributes or self.add_node_age_summaries_as_node_annotations) and node_age_summaries:
                for fieldname, stats_fieldname in zip(self.node_age_summaries_fieldnames, self.summary_stats_fieldnames):
                    self._decorate(
                        target=node,
                        fieldname=fieldname,
                        value=node_age_summaries[split_bitmask].get(stats_fieldname, 0.0),
                        set_attribute=self.add_node_age_summaries_as_node_attributes,
                        set_annotation=self.add_node_age_summaries_as_node_annotations,
                        )
            if (self.add_edge_length_summaries_as_edge_attributes or self.add_edge_length_summaries_as_edge_annotations) and edge_length_summaries:
                for fieldname, stats_fieldname in zip(self.edge_length_summaries_fieldnames, self.summary_stats_fieldnames):
                    self._decorate(
                        target=node.edge,
                        fieldname=fieldname,
                        value=edge_length_summaries[split_bitmask].get(stats_fieldname, 0.0),
                        set_attribute=self.add_edge_length_summaries_as_edge_attributes,
                        set_annotation=self.add_edge_length_summaries_as_edge_annotations,
                        )
            if self.set_edge_lengths is None:
                pass
            elif self.set_edge_lengths == "support":
                node.edge.length = split_support
            else:
                if self.set_edge_lengths in ("mean-age", "median-age"):
                    if not node_age_summaries:
                        raise ValueError("Node ages not available")
                    if self.set_edge_lengths == "mean-age":
                        node.age = node_age_summaries[split_bitmask]["mean"]
                    elif self.set_edge_lengths == "median-age":
                        node.age = node_age_summaries[split_bitmask]["median"]
                    else:
                        raise ValueError(self.set_edge_lengths)
                elif self.set_edge_lengths in ("mean-length", "median-length"):
                    if not edge_length_summaries:
                        raise ValueError("Edge lengths not available")
                    if self.set_edge_lengths == "mean-length":
                        node.edge.length = edge_length_summaries[split_bitmask]["mean"]
                    elif self.set_edge_lengths == "median-length":
                        node.edge.length = edge_length_summaries[split_bitmask]["median"]
                    else:
                        raise ValueError(self.set_edge_lengths)
                    if self.minimum_edge_length is not None and edge.length < self.minimum_edge_length:
                        edge.length = self.minimum_edge_length
                else:
                    raise ValueError(self.set_edge_lengths)
        if self.set_edge_lengths in ("mean-age", "median-age"):
            tree.set_edge_lengths_from_node_ages(
                    minimum_edge_length=self.minimum_edge_length,
                    error_on_negative_edge_lengths=self.error_on_negative_edge_lengths)
        return tree

###############################################################################
### TreeArray

class TreeArray(taxonmodel.TaxonNamespaceAssociated):
    """
    Storage of minimal tree structural information as represented by toplogy
    and edge lengths.

    This class stores trees as collections of splits and edge lengths. All
    other information, such as labels, metadata annotations, etc. will be
    discarded. A full :class:`Tree` instance can be reconstructed as needed
    from the structural information stored by this class.
    """

    class IncompatibleTreeArrayUpdate(Exception):
        pass
    class IncompatibleRootingTreeArrayUpdate(IncompatibleTreeArrayUpdate):
        pass
    class IncompatibleEdgeLengthsTreeArrayUpdate(IncompatibleTreeArrayUpdate):
        pass
    class IncompatibleNodeAgesTreeArrayUpdate(IncompatibleTreeArrayUpdate):
        pass
    class IncompatibleTreeWeightsTreeArrayUpdate(IncompatibleTreeArrayUpdate):
        pass

    ##############################################################################
    ## Life-Cycle

    def __init__(self,
            taxon_namespace=None,
            is_rooted_trees=None,
            ignore_edge_lengths=False,
            ignore_node_ages=True,
            use_tree_weights=True,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            ):
        """
        Parameters
        ----------
        taxon_namespace : :class:`TaxonNamespace`
            The operational taxonomic unit concept namespace to manage taxon
            references.
        is_rooted_trees : bool
            If not set, then it will be set based on the rooting state of the
            first tree added. If `True`, then trying to add an unrooted tree
            will result in an error. If `False`, then trying to add a rooted
            tree will result in an error.
        ignore_edge_lengths : bool
            If `True`, then edge lengths of splits will not be stored. If
            `False`, then edge lengths will be stored.
        ignore_node_ages : bool
            If `True`, then node ages of splits will not be stored. If
            `False`, then node ages will be stored.
        use_tree_weights : bool
            If `False`, then tree weights will not be used to weight splits.
        """
        taxonmodel.TaxonNamespaceAssociated.__init__(self,
                taxon_namespace=taxon_namespace)

        # Configuration
        self._is_rooted_trees = is_rooted_trees
        self.ignore_edge_lengths = ignore_edge_lengths
        self.ignore_node_ages = ignore_node_ages
        self.use_tree_weights = use_tree_weights
        self.default_edge_length_value = 0 # edge.length of `None` gets this value
        self.tree_type = treemodel.Tree

        # Storage
        self._tree_split_bitmasks = []
        self._tree_edge_lengths = []
        self._tree_leafset_bitmasks = []
        self._split_distribution = SplitDistribution(
                taxon_namespace=self.taxon_namespace,
                ignore_edge_lengths=self.ignore_edge_lengths,
                ignore_node_ages=self.ignore_node_ages,
                ultrametricity_precision=ultrametricity_precision,
                )

    ##############################################################################
    ## Book-Keeping

    def _get_is_rooted_trees(self):
        return self._is_rooted_trees
    is_rooted_trees = property(_get_is_rooted_trees)

    def _get_split_distribution(self):
        return self._split_distribution
    split_distribution = property(_get_split_distribution)

    def validate_rooting(self, rooting_of_other):
        if self._is_rooted_trees is None:
            self._is_rooted_trees = rooting_of_other
        elif self._is_rooted_trees != rooting_of_other:
            if self._is_rooted_trees:
                ta = "rooted"
                t = "unrooted"
            else:
                ta = "unrooted"
                t = "rooted"
            raise error.MixedRootingError("Cannot add {tree_rooting} tree to TreeArray with {tree_array_rooting} trees".format(
                tree_rooting=t,
                tree_array_rooting=ta))

    ##############################################################################
    ## Updating from Another TreeArray

    def update(self, other):
        if len(self) > 0:
            # self.validate_rooting(other._is_rooted_trees)
            if self._is_rooted_trees is not other._is_rooted_trees:
                raise TreeArray.IncompatibleRootingTreeArrayUpdate("Updating from incompatible TreeArray: 'is_rooted_trees' should be '{}', but is instead '{}'".format(other._is_rooted_trees, self._is_rooted_trees, ))
            if self.ignore_edge_lengths is not other.ignore_edge_lengths:
                raise TreeArray.IncompatibleEdgeLengthsTreeArrayUpdate("Updating from incompatible TreeArray: 'ignore_edge_lengths' is not: {} ".format(other.ignore_edge_lengths, self.ignore_edge_lengths, ))
            if self.ignore_node_ages is not other.ignore_node_ages:
                raise TreeArray.IncompatibleNodeAgesTreeArrayUpdate("Updating from incompatible TreeArray: 'ignore_node_ages' should be '{}', but is instead '{}'".format(other.ignore_node_ages, self.ignore_node_ages))
            if self.use_tree_weights is not other.use_tree_weights:
                raise TreeArray.IncompatibleTreeWeightsTreeArrayUpdate("Updating from incompatible TreeArray: 'use_tree_weights' should be '{}', but is instead '{}'".format(other.use_tree_weights, self.use_tree_weights))
        else:
            self._is_rooted_trees = other._is_rooted_trees
            self.ignore_edge_lengths = other.ignore_edge_lengths
            self.ignore_node_ages = other.ignore_node_ages
            self.use_tree_weights = other.use_tree_weights
        self._tree_split_bitmasks.extend(other._tree_split_bitmasks)
        self._tree_edge_lengths.extend(other._tree_edge_lengths)
        self._tree_leafset_bitmasks.extend(other._tree_leafset_bitmasks)
        self._split_distribution.update(other._split_distribution)

    ##############################################################################
    ## Fundamental Tree Accession

    def add_tree(self,
            tree,
            is_bipartitions_updated=False,
            index=None):
        """
        Adds the structure represented by a :class:`Tree` instance to the
        collection.

        Parameters
        ----------
        tree : :class:`Tree`
            A :class:`Tree` instance. This must have the same rooting state as
            all the other trees accessioned into this collection as well as
            that of `self.is_rooted_trees`.
        is_bipartitions_updated : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.
        index : integer
            Insert before index.

        Returns
        -------
        index : int
            The index of the accession.
        s : iterable of splits
            A list of split bitmasks from `tree`.
        e :
            A list of edge length values from `tree`.
        """
        if self.taxon_namespace is not tree.taxon_namespace:
            raise error.TaxonNamespaceIdentityError(self, tree)
        self.validate_rooting(tree.is_rooted)
        # if self._is_rooted_trees is None:
        #     self._is_rooted_trees = tree.is_rooted
        # else:
        #     if self._is_rooted_trees != tree.is_rooted:
        #         if self._is_rooted_trees:
        #             ta = "rooted"
        #             t = "unrooted"
        #         else:
        #             ta = "unrooted"
        #             t = "rooted"
        #         raise error.MixedRootingError("Cannot add {tree_rooting} tree to TreeArray of {tree_array_rooting} trees '{}'".format(
        #             tree_rooting=t,
        #             tree_array_rootin=ta))
        splits, edge_lengths, node_ages = self._split_distribution.count_splits_on_tree(
                tree=tree,
                is_bipartitions_updated=is_bipartitions_updated,
                default_edge_length_value=self.default_edge_length_value)
        if index is None:
            self._tree_leafset_bitmasks.append(tree.seed_node.edge.bipartition.leafset_bitmask)
        else:
            self._tree_leafset_bitmasks.insert(index,
                    tree.seed_node.edge.bipartition.leafset_bitmask)
        return self.add_splits(splits=splits,
                edge_lengths=edge_lengths,
                index=index)

    def add_splits(self, splits, edge_lengths, index=None):
        """
        Adds a "tree" as represented by a list of splits and edge lengths.

        Parameters
        ----------
        splits : iterable of ints
            An iterable of split bitmasks.
        edge_lengths : iterable of values
            An iterable of (usually numeric) values for each split listed in ``splits``.
        index : integer
            Insert before index.

        """
        splits = tuple(splits)
        if self.ignore_edge_lengths:
            # edge_lengths = tuple( [None] * len(splits) )
            edge_lengths = tuple( None for x in range(len(splits)) )
        else:
            assert len(splits) == len(edge_lengths), "Unequal vectors:\n    Splits: {}\n    Edges: {}\n".format(splits, edge_lengths)
            edge_lengths = tuple(edge_lengths)
        if index is None:
            index = len(self._tree_split_bitmasks)
            self._tree_split_bitmasks.append(splits)
            self._tree_edge_lengths.append(edge_lengths)
        else:
            self._tree_split_bitmasks.insert(index, splits)
            self._tree_edge_lengths.insert(index, edge_lengths)
        return index, splits, edge_lengths

    def add_trees(self, trees, is_bipartitions_updated=False):
        """
        Adds multiple structures represneted by an iterator over or iterable of
        :class:`Tree` instances to the collection.

        Parameters
        ----------
        trees : iterator over or iterable of :class:`Tree` instances
            An iterator over or iterable of :class:`Tree` instances. Thess must
            have the same rooting state as all the other trees accessioned into
            this collection as well as that of `self.is_rooted_trees`.
        is_bipartitions_updated : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.

        """
        for tree in trees:
            self.add_tree(tree,
                    is_bipartitions_updated=is_bipartitions_updated)

    def read_from_files(self,
            files,
            schema,
            **kwargs):
        """
        Adds multiple structures from one or more external file sources to the
        collection.

        Parameters
        ----------
        files : iterable of strings and/or file objects
            A list or some other iterable of file paths or file-like objects
            (string elements will be assumed to be paths to files, while all
            other types of elements will be assumed to be file-like
            objects opened for reading).
        schema : string
            The data format of the source. E.g., "nexus", "newick", "nexml".
        \*\*kwargs : keyword arguments
            These will be passed directly to the underlying schema-specific
            reader implementation.
        """
        if "taxon_namespace" in kwargs:
            if kwargs["taxon_namespace"] is not self.taxon_namespace:
                raise ValueError("TaxonNamespace object passed as keyword argument is not the same as self's TaxonNamespace reference")
            kwargs.pop("taxon_namespace")
        target_tree_offset = kwargs.pop("tree_offset", 0)
        tree_yielder = self.tree_type.yield_from_files(
                files=files,
                schema=schema,
                taxon_namespace=self.taxon_namespace,
                **kwargs)
        current_source_index = None
        current_tree_offset = None
        for tree_idx, tree in enumerate(tree_yielder):
            current_yielder_index = tree_yielder.current_file_index
            if current_source_index != current_yielder_index:
                current_source_index = current_yielder_index
                current_tree_offset = 0
            if current_tree_offset >= target_tree_offset:
                self.add_tree(tree=tree, is_bipartitions_updated=False)
            current_tree_offset += 1

    ##############################################################################
    ## Convenient I/O

    def read_from_stream(self, fileobj, schema, **kwargs):
        """
        Reads trees from a file. See :meth:`TreeList.read_from_stream()`.
        """
        return self.read_from_files(
                files=[fileobj],
                schema=schema,
                **kwargs)

    def read_from_path(self, filepath, schema, **kwargs):
        """
        Reads trees from a path. See :meth:`TreeList.read_from_path()`.
        """
        return self.read_from_files(
                files=[filepath],
                schema=schema,
                **kwargs)

    def read_from_string(self, src_str, schema, **kwargs):
        """
        Reads trees from a string. See :meth:`TreeList.read_from_string()`.
        """
        return self.read_from_files(
                files=[StringIO(src_str)],
                schema=schema,
                **kwargs)

    ##############################################################################
    ## Container (List) Interface

    def append(tree, is_bipartitions_updated=False):
        """
        Adds a :class:`Tree` instance to the collection before position given
        by `index`.

        Parameters
        ----------
        tree : :class:`Tree`
            A :class:`Tree` instance. This must have the same rooting state as
            all the other trees accessioned into this collection as well as
            that of `self.is_rooted_trees`.
        is_bipartitions_updated : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.

        """
        return self.add_tree(tree=tree,
                is_bipartitions_updated=is_bipartitions_updated)

    def insert(index, tree, is_bipartitions_updated=False):
        """
        Adds a :class:`Tree` instance to the collection before position given
        by `index`.

        Parameters
        ----------
        index : integer
            Insert before index.
        tree : :class:`Tree`
            A :class:`Tree` instance. This must have the same rooting state as
            all the other trees accessioned into this collection as well as
            that of `self.is_rooted_trees`.
        is_bipartitions_updated : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.

        Returns
        -------
        index : int
            The index of the accession.
        s : iterable of splits
            A list of split bitmasks from `tree`.
        e :
            A list of edge length values `tree`.
        """
        return self.add_tree(tree=tree,
                is_bipartitions_updated=is_bipartitions_updated,
                index=index)

    def extend(self, tree_array):
        """
        Accession of data from `tree_array` to self.

        Parameters
        ----------
        tree_array : :class:`TreeArray`
            A :class:`TreeArray` instance from which to add data.

        """
        assert self.taxon_namespace is tree_array.taxon_namespace
        assert self._is_rooted_trees is tree_array._is_rooted_trees
        assert self.ignore_edge_lengths is tree_array.ignore_edge_lengths
        assert self.ignore_node_ages is tree_array.ignore_node_ages
        assert self.use_tree_weights is tree_array.use_tree_weights
        self._tree_split_bitmasks.extend(tree_array._tree_split_bitmasks)
        self._tree_edge_lengths.extend(tree_array._tree_edge_lengths)
        self._split_distribution.update(tree_array._split_distribution)
        return self

    def __iadd__(self, tree_array):
        """
        Accession of data from `tree_array` to self.

        Parameters
        ----------
        tree_array : :class:`TreeArray`
            A :class:`TreeArray` instance from which to add data.

        """
        return self.extend(tree_array)

    def __add__(self, other):
        """
        Creates and returns new :class:`TreeArray`.

        Parameters
        ----------
        other : iterable of :class:`Tree` objects

        Returns
        -------
        tlist : :class:`TreeArray` object
            :class:`TreeArray` object containing clones of :class:`Tree` objects
            in `self` and `other`.
        """
        ta = TreeArray(
                taxon_namespace=self.taxon_namespace,
                is_rooted_trees=self._is_rooted_trees,
                ignore_edge_lengths=self.ignore_edge_lengths,
                ignore_node_ages=self.ignore_node_ages,
                use_tree_weights=self.use_tree_weights)
        ta += self
        ta += other
        return ta

    def __contains__(self, splits):
        # expensive!!
        return tuple(splits) in self._tree_split_bitmasks

    def __delitem__(self, index):
        raise NotImplementedError
        # expensive!!
        # tree_split_bitmasks = self._trees_splits[index]
        ### TODO: remove this "tree" from underlying splits distribution
        # for split in tree_split_bitmasks:
        #   self._split_distribution.split_counts[split] -= 1
        # etc.
        # becomes complicated because tree weights need to be updated etc.
        # del self._tree_split_bitmasks[index]
        # del self._tree_edge_lengths[index]
        # return

    def __iter__(self):
        """
        Yields pairs of (split, edge_length) from the store.
        """
        for split, edge_length in zip(self._tree_split_bitmasks, self._tree_edge_lengths):
            yield split, edge_length

    def __reversed__(self):
        raise NotImplementedError

    def __len__(self):
        return len(self._tree_split_bitmasks)

    def __getitem__(self, index):
        raise NotImplementedError
        # """
        # Returns a pair of tuples, ( (splits...), (lengths...) ), corresponding
        # to the "tree" at `index`.
        # """
        # return self._tree_split_bitmasks[index], self._tree_edge_lengths[index]

    def __setitem__(self, index, value):
        raise NotImplementedError

    def clear(self):
        raise NotImplementedError
        self._tree_split_bitmasks = []
        self._tree_edge_lengths = []
        self._tree_leafset_bitmasks = []
        self._split_distribution.clear()

    def index(self, splits):
        raise NotImplementedError
        return self._tree_split_bitmasks.index(splits)

    def pop(self, index=-1):
        raise NotImplementedError

    def remove(self, tree):
        raise NotImplementedError

    def reverse(self):
        raise NotImplementedError

    def sort(self, key=None, reverse=False):
        raise NotImplementedError

    ##############################################################################
    ## Accessors/Settors

    def get_split_bitmask_and_edge_tuple(self, index):
        """
        Returns a pair of tuples, ( (splits...), (lengths...) ), corresponding
        to the "tree" at `index`.
        """
        return self._tree_split_bitmasks[index], self._tree_edge_lengths[index]

    ##############################################################################
    ## Calculations

    def calculate_log_product_of_split_supports(self,
            include_external_splits=False,
            ):
        """
        Calculates the log product of split support for all trees in the
        collection. The tree with the maximum log product of the split support
        is the 'maximum credibility tree' (MCT) or the 'maximum clade
        credibility tree' (MCCT), though sometimes one or both these terms
        are used to refer to the tree with the highest *sum* of split support.

        Parameters
        ----------
        include_external_splits : bool
            If `True`, then non-internal split posteriors will be included in
            the score. Defaults to `False`: these are skipped. This should only
            make a difference when dealing with splits collected from trees of
            different leaf sets.

        Returns
        -------
        s : tuple(list[numeric], integer)
            Returns a tuple, with the first element being the list of scores
            and the second being the index of the highest score. The element order
            corresponds to the trees accessioned in the collection.
        """
        assert len(self._tree_leafset_bitmasks) == len(self._tree_split_bitmasks)
        scores = []
        max_score = None
        max_score_tree_idx = None
        split_frequencies = self._split_distribution.split_frequencies
        for tree_idx, (leafset_bitmask, split_bitmasks) in enumerate(zip(self._tree_leafset_bitmasks, self._tree_split_bitmasks)):
            log_product_of_split_support = 0.0
            for split_bitmask in split_bitmasks:
                if include_external_splits or not treemodel.Bipartition.is_trivial_bitmask(split_bitmask, leafset_bitmask):
                    split_support = split_frequencies.get(split_bitmask, 0.0)
                    if split_support:
                        log_product_of_split_support += math.log(split_support)
            if max_score is None or max_score < log_product_of_split_support:
                max_score = log_product_of_split_support
                max_score_tree_idx = tree_idx
            scores.append(log_product_of_split_support)
        return scores, max_score_tree_idx

    def maximum_product_of_split_support_tree(self,
            include_external_splits=False,
            tree_class=None,
            **split_summarization_kwargs
            ):
        """
        Return the tree with that maximizes the product of split supports.
        The tree with the maximum log product of the split support is the
        'maximum credibility tree' (MCT) or the 'maximum clade credibility
        tree' (MCCT), though sometimes one or both these terms are used to
        refer to the tree with the highest *sum* of split support.

        Parameters
        ----------
        include_external_splits : bool
            If `True`, then non-internal split posteriors will be included in
            the score. Defaults to `False`: these are skipped. This should only
            make a difference when dealing with splits collected from trees of
            different leaf sets.

        Returns
        -------
        mct_tree : Tree
            Tree that maximizes the product of split supports.
        """
        scores, max_score_tree_idx = self.calculate_log_product_of_split_supports(
                include_external_splits=include_external_splits,
                )
        tree = self.restore_tree(
                index=max_score_tree_idx,
                tree_class=tree_class,
                **split_summarization_kwargs)
        tree.log_product_of_split_support = scores[max_score_tree_idx]
        self._split_distribution.summarize_splits_on_tree(
            tree=tree,
            is_bipartitions_updated=True,
            **split_summarization_kwargs
            )
        return tree

    def calculate_sum_of_split_supports(self,
            include_external_splits=False,
            ):
        """
        Calculates the *sum* of split support for all trees in the
        collection.

        Parameters
        ----------
        include_external_splits : bool
            If `True`, then non-internal split posteriors will be included in
            the score. Defaults to `False`: these are skipped. This should only
            make a difference when dealing with splits collected from trees of
            different leaf sets.

        Returns
        -------
        s : tuple(list[numeric], integer)
            Returns a tuple, with the first element being the list of scores
            and the second being the index of the highest score. The element order
            corresponds to the trees accessioned in the collection.
        """
        assert len(self._tree_leafset_bitmasks) == len(self._tree_split_bitmasks)
        scores = []
        max_score = None
        max_score_tree_idx = None
        split_frequencies = self._split_distribution.split_frequencies
        for tree_idx, (leafset_bitmask, split_bitmasks) in enumerate(zip(self._tree_leafset_bitmasks, self._tree_split_bitmasks)):
            sum_of_support = 0.0
            for split_bitmask in split_bitmasks:
                if include_external_splits or not treemodel.Bipartition.is_trivial_bitmask(split_bitmask, leafset_bitmask):
                    split_support = split_frequencies.get(split_bitmask, 0.0)
                    sum_of_support += split_support
            if max_score is None or max_score < sum_of_support:
                max_score = sum_of_support
                max_score_tree_idx = tree_idx
            scores.append(sum_of_support)
        return scores, max_score_tree_idx

    def maximum_sum_of_split_support_tree(self,
            include_external_splits=False,
            tree_class=None,
            **split_summarization_kwargs
            ):
        """
        Return the tree with that maximizes the sum of split supports.
        The tree with the maximum log sum of the split support is the
        'maximum credibility tree' (MCT) or the 'maximum clade credibility
        tree' (MCCT), though sometimes one or both these terms are used to
        refer to the tree with the highest *sum* of split support.

        Parameters
        ----------
        include_external_splits : bool
            If `True`, then non-internal split posteriors will be included in
            the score. Defaults to `False`: these are skipped. This should only
            make a difference when dealing with splits collected from trees of
            different leaf sets.

        Returns
        -------
        mct_tree : Tree
            Tree that maximizes the sum of split supports.
        """
        scores, max_score_tree_idx = self.calculate_sum_of_split_supports(
                include_external_splits=include_external_splits,
                )
        tree = self.restore_tree(
                index=max_score_tree_idx,
                tree_class=tree_class,
                **split_summarization_kwargs
                )
        tree.sum_of_split_support = scores[max_score_tree_idx]
        self._split_distribution.summarize_splits_on_tree(
            tree=tree,
            is_bipartitions_updated=True,
            **split_summarization_kwargs
            )
        return tree

    def consensus_tree(self,
            min_freq=constants.GREATER_THAN_HALF,
            **split_summarization_kwargs
            ):
        """
        Returns a consensus tree from splits in `self`.

        Parameters
        ----------

        min_freq : real
            The minimum frequency of a split in this distribution for it to be
            added to the tree.

        is_rooted : bool
            Should tree be rooted or not? If *all* trees counted for splits are
            explicitly rooted or unrooted, then this will default to `True` or
            `False`, respectively. Otherwise it defaults to `None`.

        \*\*split_summarization_kwargs : keyword arguments
            These will be passed directly to the underlying
            :class:`SplitDistributionSummarizer` object. See
            :meth:`SplitDistributionSummarizer.configure` for options.

        Returns
        -------
        t : consensus tree

        """
        return self._split_distribution.consensus_tree(
                min_freq=min_freq,
                is_rooted=self.is_rooted,
                **split_summarization_kwargs
                )

    ##############################################################################
    ## Tree Reconstructions

    def restore_tree(self,
            index,
            tree_class=None,
            summarize_splits_on_tree=False,
            **split_summarization_kwargs
            ):
        split_bitmasks = self._tree_split_bitmasks[index]
        if self.ignore_edge_lengths:
            split_edge_lengths = None
        else:
            assert len(self._tree_split_bitmasks) == len(self._tree_edge_lengths)
            edge_lengths = self._tree_edge_lengths[index]
            split_edge_lengths = dict(zip(split_bitmasks, edge_lengths))
        if tree_class is None:
            tree_class = self.tree_type
        tree = tree_class.from_split_bitmasks(
                split_bitmasks=split_bitmasks,
                taxon_namespace=self.taxon_namespace,
                is_rooted=self._is_rooted_trees,
                split_edge_lengths=split_edge_lengths,
                )
        # if update_bipartitions:
        #     tree.encode_bipartitions()
        if summarize_splits_on_tree:
            split_summarization_kwargs["is_bipartitions_updated"] = True
            self._split_distribution.summarize_splits_on_tree(
                    tree=tree,
                    **split_summarization_kwargs)
        return tree

    def consensus_tree(self, *args, **kwargs):
        return self._split_distribution.consensus_tree(*args, **kwargs)


