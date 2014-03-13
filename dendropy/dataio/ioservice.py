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

        **kwargs : keyword arguments

        """
        IOService.__init__(self, **kwargs)

    def read(self,
            stream,
            dataset=None,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None):
        raise NotImplementedError

    def read_dataset(self,
            stream,
            dataset,
            exclude_trees=False,
            exclude_chars=False):
        """
        Populates the given `DataSet` object from external data source.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        dataset : DataSet object
            The target `DataSet` to populate/build.

        exclude_trees : boolean, default = False
            If set to `True`, tree data will not be read from the source.

        exclude_chars : boolean, default = False
            If set to `True`, character data will not be read from the source.

        """
        taxon_namespace_factory = dataset.new_taxon_namespace
        if exclude_chars:
            char_matrix_factory = None
        else:
            char_matrix_factory = dataset.new_char_matrix
        if exclude_trees:
            tree_list_factory = None
        else:
            tree_list_factory = dataset.new_tree_list


    def read_trees(self,
            stream,
            taxon_namespace_factory,
            tree_list_factory):
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

        Returns
        -------
        List of `TreeList` objects.

        """
        # `product` is a namedtuple("DataReaderProducts", ["taxon_namespaces", "tree_lists", "char_matrices"])
        product = self.read(stream=stream,
                dataset=None,
                taxon_namespace_factory=taxon_namespace_factory,
                tree_list_factory=tree_list_factory,
                char_matrix_factory=None)
        return product.tree_lists

    def read_char_matrix(self, stream, char_matrix):
        """
        Populates the given `CharacterMatrixList` object with
        sequences/alignments from an external data source.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        char_matrix_list : CharacterMatrixList object
            The collection of CharacterMatrix objects into which data will be read.

        """
        taxon_namespace_factory = lambda : char_matrix.taxon_namespace
        char_matrix_factory = lambda : char_matrix
        tree_list_factory = None

    # def get_default_taxon_namespace(self, dataset, **kwargs):
    #     """
    #     Returns an appropriate TaxonNamespace object, based on current settings.
    #     """
    #     if dataset is None:
    #         raise TypeError("'dataset' is not defined")
    #     if "taxon_namespace" in kwargs:
    #         self.attached_taxon_namespace = kwargs["taxon_namespace"]
    #     if dataset.attached_taxon_namespace is not None:
    #         if self.attached_taxon_namespace is not None \
    #                 and dataset.attached_taxon_namespace is not self.attached_taxon_namespace:
    #             raise TypeError("DataSet is attached to different TaxonNamespace than that specified by 'taxon_namespace'")
    #         self.attached_taxon_namespace = dataset.attached_taxon_namespace
    #     if self.attached_taxon_namespace is not None:
    #         if self.attached_taxon_namespace not in dataset.taxon_namespaces:
    #             dataset.add(self.attached_taxon_namespace)
    #         return self.attached_taxon_namespace
    #     else:
    #         return dataset.new_taxon_namespace(**kwargs)

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
        """
        IOService.__init__(self, **kwargs)

    def write_dataset(self, dataset, stream):
        raise NotImplementedError

    def write_tree_list(self, tree_list, stream):
        raise NotImplementedError

    def write_char_matrices(self, char_matrix_list, stream):
        raise NotImplementedError

