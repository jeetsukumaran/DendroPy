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

import collections
import warnings

###############################################################################
## IOService

class IOService(object):
    """
    Base class for all readers/writers.
    """

    def __init__(self, **kwargs):
        pass

###############################################################################
## DataReader

class DataReader(IOService):
    """
    Base class for all readers.

    Consumes a stream and builds or composes corresponding phylogenetic data
    object or objects. Abstract class, to be implemented by derived classes
    specializing in particular data formats.
    """

    Product = collections.namedtuple(
            "product",
            ["taxon_namespaces", "tree_lists", "char_matrices"]
            )

    def __init__(self, **kwargs):
        """
        Constructs and configures a `DataReader` object by "harvesting" keyword
        arguments and setting state accordingly. Keyword arguments recognized
        and processed will be removed from the keyword argument dictionary.

        Parameters
        ----------

        **kwargs : schema- and implementation-specific keyword arguments

        """
        IOService.__init__(self, **kwargs)

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            global_annotations_target=None):
        """
        Deriving classes should implement this method to build a data product
        from the information in `stream` using the provided factory functions.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        taxon_namespace_factory : function object
            A function that takes one named argument, `label`, and returns a
            `TaxonNamespace` object to be used for each distinct block of
            operational taxonomic unit concept definitions encountered in the
            data source.

            The function will be called as::

                tns = taxon_namespace_factory(label="label")

            In the simplest case, a new `TaxonNamespace` object can be created
            for each block of taxon definitions in the data source by setting
            the factory function to::

                taxon_namespace_factory = TaxonNamespace

            If all data objects are to be organized into a DataSet object,
            then:

                taxon_namespace_factory = dataset.new_taxon_namespace

            If all data objects should reference the *same* `TaxonNamespace`
            object, then:

                taxon_namespace_factory = lambda label : taxon_namespace

            where `taxon_namespace` is an existing `TaxonNamespace` object that
            should be used.

            If `taxon_namespace_factor` is `None`, then no tree data will be
            parsed.

        tree_list_factory : function object
            A function that takes two named arguments, `label` and
            `taxon_namespace`, and returns a `TreeList` or equivalent object to
            be used to manage each distinct collection of trees in the data
            source.

            The function will be called as::

                tns = taxon_namespace_factory(label="label")
                tlist = tree_list_factory(label="label", taxon_namespace=tns)

            In the simplest case, a new `TreeList` object can be created for
            each block of tree definitions in the data source by setting the
            factory function to::

                tree_list_factory = TreeList

            If all data objects are to be organized into a DataSet object,
            then:

                tree_list = dataset.new_tree_list

            If all Tree data objects instantiated should be accessioned into
            the *same* `TreeList` object, then:

                taxon_namespace_factory = lambda label : tree_list.taxon_namespace
                tree_list_factory = lambda label, taxon_namespace : tree_list

            where `tree_list` is an existing `TreeList` object that should be
            used.

        char_matrix_factory : function object
            A function that takes two named arguments, `label` and
            `taxon_namespace`, and returns a `CharacterMatrix` or equivalent object to
            be used to manage each aligment or distinct set of sequences in the data
            source.

            The function will be called as::

                tns = taxon_namespace_factory(label="label")
                cm = char_matrix_factory(label="label", taxon_namespace=tns)

            In the simplest case, a new `CharacterMatrix` object can be created for
            each alignment or set of sequences in the data source by setting the
            factory function to, for e.g.::

                char_matrix_factory = DnaCharacterMatrix

            If all data objects are to be organized into a DataSet object,
            then:

                char_matrix = dataset.new_char_matrix

            If `char_matrix_factory` is `None`, then no character data will be
            parsed.

        global_annotations_target : `Annotable` object
            Any object that will be the target (or subject, in the grammatical
            sense) of general metadata or annotations in the data source. If
            `None`, then such metadata or annotations will not be stored.

        Returns
        -------

        A `Product` object : a `namedtuple` with the following attributes:
            "taxon_namespaces", "tree_lists", "char_matrices".

        """
        raise NotImplementedError

    def read_dataset(self,
            stream,
            dataset,
            exclude_trees=False,
            exclude_chars=False,
            global_annotations_target=None):
        """
        Populates the given :class:`DataSet` object from external data source.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.
        dataset : :class:`DataSet` object
            The target :class:`DataSet` to populate/build.
        exclude_trees : boolean, default: False
            If set to `True`, tree data will not be read from the source.
        exclude_chars : boolean, default: False
            If set to `True`, character data will not be read from the source.
        global_annotations_target : :class:`Annotable` object
            Any object that will be the target (or subject, in the grammatical
            sense) of general metadata or annotations in the data source. If
            `None`, then such metadata or annotations will not be stored.

        """
        raise NotImplementedError
        taxon_namespace_factory = dataset.new_taxon_namespace
        if exclude_chars:
            char_matrix_factory = None
        else:
            char_matrix_factory = dataset.new_char_matrix
        if exclude_trees:
            tree_list_factory = None
        else:
            tree_list_factory = dataset.new_tree_list
        global_annotations_target = dataset

    def read_trees(self,
            stream,
            taxon_namespace_factory,
            tree_list_factory,
            global_annotations_target=None):
        """
        Reads tree data from source into tree objects.

        With data schemas that support the concept of multiple distinct blocks
        or sets of trees (e.g. NEXUS or NeXML), each tree block will be
        accessioned into a separate `TreeList` instantiated by calling
        `tree_list_factory(label)`. If trees should be accessioned into the
        same `TreeList`, then this can be coerced by, e.g.::

            t = TreeList()
            reader.read_trees(
                stream=stream,
                taxon_namespace_factory=lambda x: t.taxon_namespace,
                tree_list_factory=lambda x : t)

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        taxon_namespace_factory : function object
            A function that takes one named argument, `label`, and returns a
            `TaxonNamespace` object to be used for each distinct block of
            operational taxonomic unit concept definitions encountered in the
            data source.

            The function will be called as::

                tns = taxon_namespace_factory(label="label")

            In the simplest case, a new `TaxonNamespace` object can be created
            for each block of taxon definitions in the data source by setting
            the factory function to::

                taxon_namespace_factory = TaxonNamespace

            If all data objects are to be organized into a DataSet object,
            then:

                taxon_namespace_factory = dataset.new_taxon_namespace

            If all data objects should reference the *same* `TaxonNamespace`
            object, then:

                taxon_namespace_factory = lambda label : taxon_namespace

            where `taxon_namespace` is an existing `TaxonNamespace` object that
            should be used.

            If `taxon_namespace_factor` is `None`, then no tree data will be
            parsed.

        tree_list_factory : function object
            A function that takes two named arguments, `label` and
            `taxon_namespace`, and returns a `TreeList` or equivalent object to
            be used to manage each distinct collection of trees in the data
            source.

            The function will be called as::

                tns = taxon_namespace_factory(label="label")
                tlist = tree_list_factory(label="label", taxon_namespace=tns)

            In the simplest case, a new `TreeList` object can be created for
            each block of tree definitions in the data source by setting the
            factory function to::

                tree_list_factory = TreeList

            If all data objects are to be organized into a DataSet object,
            then:

                tree_list = dataset.new_tree_list

            If all Tree data objects instantiated should be accessioned into
            the *same* `TreeList` object, then:

                taxon_namespace_factory = lambda label : tree_list.taxon_namespace
                tree_list_factory = lambda label, taxon_namespace : tree_list

            where `tree_list` is an existing `TreeList` object that should be
            used.

        global_annotations_target : `Annotable` object
            Any object that will be the target (or subject, in the grammatical
            sense) of general metadata or annotations in the data source. If
            `None`, then such metadata or annotations will not be stored.

        Returns
        -------
        List of `TreeList` objects.

        """
        # `product` is a namedtuple("DataReaderProducts", ["taxon_namespaces", "tree_lists", "char_matrices"])
        product = self._read(stream=stream,
                taxon_namespace_factory=taxon_namespace_factory,
                tree_list_factory=tree_list_factory,
                char_matrix_factory=None,
                global_annotations_target=global_annotations_target)
        return product.tree_lists

###############################################################################
## DataReader

class DataWriter(IOService):
    """
    Base class for all writers.

    Writes a DendroPy phylogenetic data object to a stream. Abstract class, to
    be implemented by derived classes specializing in particular data formats.
    """

    def __init__(self, **kwargs):
        """
        Constructs and configures a `DataWriter` object by "harvesting" keyword
        arguments and setting state accordingly. Keyword arguments recognized
        and processed will be removed from the keyword argument dictionary.

        Parameters
        ----------

        **kwargs : schema- and implementation-specific keyword arguments

        """
        IOService.__init__(self, **kwargs)

    def _write(self,
            stream,
            taxon_namespaces=None,
            tree_lists=None,
            char_matrices=None,
            global_annotations_target=None):
        """
        Deriving classes should implement this method to information to `stream`
        in schema-specific formatting.

        Parameters
        ----------

        stream : file or file-like object
            Destination for data.
        taxon_namespaces : Iterable of :class:`TaxonNamespace` objects
            Collection of :class:`TaxonNamespace` objects to be serialized.
        tree_lists : Iterable of :class:`TreeList` objects
            Collection of :class:`TreeList` objects to be serialized.
        char_matrices : Iterable of :class:`CharacterMatrix` objects
            Collection of :class:`CharacterMatrix` objects to be serialized.
        global_annotations_target : `Annotable` object
            Any object that will be the source (or subject, in the grammatical
            sense) of general metadata or annotations for the data. If
            `None`, then such metadata or annotations will not be stored.

        """
        raise NotImplementedError

    def write_dataset(self,
            dataset,
            stream,
            exclude_trees,
            exclude_chars,
            ):
        """
        Writes the given :class:`DataSet` object to the file-like object `stream`.

        stream : file or file-like object
            Destination for data.
        dataset : :class:`DataSet` object
            The :class:`DataSet` to write.
        exclude_trees : boolean, default: False
            If set to `True`, tree data will not be written to the destination.
        exclude_chars : boolean, default: False
            If set to `True`, character data will not be written to the destination.
        global_annotations_target : :class:`Annotable` object
            Any object that will be the source (or subject, in the grammatical
            sense) of general metadata or annotations for the data. If
            `None`, then such metadata or annotations will not be stored.
        """
        tree_lists = dataset.tree_lists if not exclude_trees else None
        char_matrices = dataset.char_matrices if not exclude_trees else None
        self.write(
                stream=stream,
                taxon_namespaces=dataset.taxon_namespaces,
                tree_lists=tree_lists,
                char_matrices=char_matrices,
                global_annotations_target=dataset)

    def write_tree_list(self, tree_list, stream):
        self._write(
                stream=stream,
                taxon_namespaces=None,
                tree_lists=[tree_list],
                char_matrices=None,
                global_annotations_target=None)

    def write_tree(self, tree, stream):
        self._write(
                stream=stream,
                taxon_namespaces=None,
                tree_lists=[[tree]],
                char_matrices=None,
                global_annotations_target=None)

    def write_char_matrices(self, char_matrix_list, stream):
        self._write(
                stream=stream,
                taxon_namespaces=None,
                tree_lists=None,
                char_matrices=char_matrix_list,
                global_annotations_target=None)

    def write_char_matrix(self, char_matrix, stream):
        self._write(
                stream=stream,
                taxon_namespaces=None,
                tree_lists=None,
                char_matrices=[char_matrix],
                global_annotations_target=None)

