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
This module handles the core definition of tree data structure class,
as well as all the structural classes that make up a tree.
"""

from dendropy.datamodel import treemodel
from dendropy.datamodel import taxonmodel
from dendropy.utility import error
from dendropy.calculate import treesplit

###############################################################################
### TreeSplitDistribution

class TreeSplitDistribution(taxonmodel.TaxonNamespaceAssociated):
    """
    Storage of tree structural information as represented by toplogy and edge
    lengths.

    This class stores trees as collections of splits and edge lengths. All
    other information, such as labels, metadata annotations, etc. will be
    discarded. A full :class:`Tree` instance can be reconstructed as needed
    from the structural information stored by this class.
    """

    ##############################################################################
    ## Life-Cycle

    def __init__(self,
            taxon_namespace=None,
            is_rooted_trees=None,
            tree_traversal_order="preorder",
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
        tree_traversal_order : str
            One of 'preorder' or 'postorder'. When decomposing trees into
            splits and edge lengths, the default order that the edges should be
            visited.
        """
        taxonmodel.TaxonNamespaceAssociated.__init__(self,
                taxon_namespace=taxon_namespace)
        self._is_rooted_trees = is_rooted_trees
        self._splits = []
        self._edge_lengths = []
        self._tree_traversal_order = None
        self._tree_traversal_method_name = None
        self.tree_traversal_order = tree_traversal_order
        self.default_edge_length_value = 0.0 # edge.length of `None` gets this value

    ##############################################################################
    ## Book-Keeping

    def _get_tree_traversal_order(self):
        return self._tree_traversal_order
    def _set_tree_traversal_order(self, tree_traversal_order):
        if tree_traversal_order == "preorder":
            self._tree_traversal_order = "preorder"
            self._tree_traversal_method_name = "preorder_edge_iter"
        elif tree_traversal_order == "postorder":
            self._tree_traversal_order = "postorder"
            self._tree_traversal_method_name = "postorder_edge_iter"
        else:
            raise ValueError(tree_traversal_order)
    tree_traversal_order = property(_get_tree_traversal_order, _set_tree_traversal_order)

    def _get_is_rooted_trees(self):
        return self._rooted_trees
    is_rooted_trees = property(_get_is_rooted_trees)

    ##############################################################################
    ## Metrics

    def __len__(self):
        return len(self._splits)

    ##############################################################################
    ## Tree Accession

    def add_splits(self, splits, edge_lengths):
        """
        Adds a "tree" as represented by a list of splits and edge lengths.

        Parameters
        ----------
        splits : iterable of ints
            An iterable of split bitmasks.
        edge_lengths : iterable of values
            An iterable of (usually numeric) values for each split listed in ``splits``.

        """
        assert len(splits) == len(edge_lengths)
        self._splits.append(tuple(splits))
        self._edge_lengths.append(tuple(edge_lengths))

    def add_tree(self, tree, is_splits_encoded=False):
        """
        Adds a :class:`Tree` instance to the collection.

        Parameters
        ----------
        tree : :class:`Tree`
            A :class:`Tree` instance. This must have the same rooting state as
            all the other trees accessioned into this collection as well as
            that of `self.is_rooted_trees`.
        is_splits_encoded : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.

        Returns
        -------
        i : int
            The index of the accession.
        s : iterable of splits
            A list of split bitmasks from `tree`.
        e :
            A list of edge length values `tree`.
        """
        if self.taxon_namespace is not tree.taxon_namespace:
            raise error.TaxonNamespaceIdentityError(self, tree)
        if self._is_rooted_trees is None:
            self._is_rooted_trees = tree.is_rooted
        else:
            assert self._is_rooted_trees == tree.is_rooted
        if not is_splits_encoded:
            tree.encode_splits()
        splits = []
        edge_lengths = []
        tree_traversal_iter = getattr(tree, self._tree_traversal_method_name)
        for edge in tree_traversal_iter():
            splits.append(edge.split_bitmask)
            if edge.length is None:
                edge_lengths.append(self.default_edge_length_value)
            else:
                edge_lengths.append(edge.length)
        idx = len(self._splits)
        self._splits.append(tuple(splits))
        self._edge_lengths.append(tuple(edge_lengths))
        return idx, splits, edge_lengths

    def add_trees(self, trees, is_splits_encoded=False):
        """
        Adds multiple :class:`Tree` instances to the collection.

        Parameters
        ----------
        trees : iterator over or iterable of :class:`Tree` instances
            An iterator over or iterable of :class:`Tree` instances. Thess must
            have the same rooting state as all the other trees accessioned into
            this collection as well as that of `self.is_rooted_trees`.
        is_splits_encoded : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.

        """
        for tree in trees:
            self.add_tree(tree,
                    is_splits_encoded=is_splits_encoded)

    def add_trees_from_file(self,
            source,
            schema,
            is_splits_encoded=False,
            **kwargs):
        """
        Adds multiple :class:`Tree` instances from an external data source.

        Parameters
        ----------
        source : string or file object
            If string, will be treated as a path name. Otherwise, will be
            assumed to be a file or file-like object opened for reading.
        schema : string
            The data format of the source. E.g., "nexus", "newick", "nexml".
        is_splits_encoded : bool
            If `False` [default], then the tree will have its splits encoded or
            updated. Otherwise, if `True`, then the tree is assumed to have its
            splits already encoded and updated.
        \*\*kwargs : keyword arguments
            These will be passed directly to the underlying schema-specific
            reader implementation.
        """
        for tree in treemodel.Tree.yield_from_files(
                files=[source],
                schema=schema,
                taxon_namespace=self.taxon_namespace,
                **kwargs):
            self.add_tree(tree=tree,
                    is_splits_encoded=is_splits_encoded)

    ##############################################################################
    ## Structural Data Access

    def __iter__(self):
        """
        Yields pairs of (split, edge_length) from the store.
        """
        for split, edge_length in zip(self._splits, self._edge_lengths):
            yield split, edge_length

    def __getitem__(self, index):
        """
        Returns pair (split, edge_length) for tree stored at `index`.
        """
        return self._splits[index], self._edge_lengths[index]

    def __setitem__(self, index, splits, edge_lengths):
        """
        Sets split and edge length at index.
        """
        self._splits[index] = tuple(splits)
        self._edge_lengths[index] = tuple(edge_lengths)

    ##############################################################################
    ## (Full) Tree Reconstruction/Access

    def restore_tree(self, index):
        pass

    def restore_tree_from_splits(self, splits, edge_lengths):
        pass

