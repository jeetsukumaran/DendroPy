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
This module handles the core definition of tree data structure class,
as well as all the structural classes that make up a tree.
"""

import collections
import math
from dendropy.utility.textprocessing import StringIO
import copy
import sys
from dendropy.utility import GLOBAL_RNG
from dendropy.utility import container
from dendropy.utility import terminal
from dendropy.utility import error
from dendropy.utility import bitprocessing
from dendropy.utility import deprecate
from dendropy.utility import constants
from dendropy.utility import textprocessing
from dendropy.datamodel import basemodel
from dendropy.datamodel import taxonmodel
from dendropy import dataio

##############################################################################
### Bipartition

class Bipartition(object):
    """
    A bipartition on a tree.

    A bipartition of a tree is a division or sorting of the leaves/tips of a
    tree into two mutually-exclusive and collectively-comprehensive subsets,
    obtained by bisecting the tree at a particular edge. There is thus a
    one-to-one correspondence with an edge of a tree and a bipartition. The
    term "split" is often also used to refer to the same concept, though this
    is typically applied to unrooted trees.

    A bipartition is modeled using a bitmask. This is a a bit array
    representing the membership of taxa, with the least-significant bit
    corresponding to the first taxon, the next least-signficant bit
    corresponding to the second taxon, and so on, till the last taxon
    corresponding to the most-significant bit. Taxon membership in one of two
    arbitrary groups, '0' or '1', is indicated by its corresponding bit being
    unset or set, respectively.

    To allow comparisons and correct identification of the same bipartition
    across different rotational and orientiational representations of unrooted
    trees, we *normalize* the bipartition such that the first taxon is always
    assigned to group '0' for bipartition representations of unrooted trees.

    The normalization of the bitmask loses information about the actual
    descendents of a particular edge. Thus in addition to the
    :attr:`Bipartition.bitmask` attribute, each |Bipartition| object
    also maintains a :attr:`Bipartition.leafset_bitmask` attribute which is
    *unnormalized*. This is a bit array representing the presence or absence of
    taxa in the subtree descending from the child node of the edge of which
    this bipartition is associated. The least-significant bit corresponds to
    the first taxon, the next least-signficant bit corresponds to the second
    taxon, and so on, with the last taxon corresponding to the most-significant
    bit. For rooted trees, the value of :attr:`Bipartition.bitmask` and
    :attr:`Bipartition.leafset_bitmask` are identical. For unrooted trees, they
    may or may not be equal.

    In general, we use :attr:`Bipartition.bitmask` data to establish the *identity*
    of a split or bipartition across *different* trees: for example, when
    computing the Robinson-Foulds distances between trees, or in assessing the
    support for different bipartitions given an MCMC or bootstrap sample of trees.
    Here the normalization of the bitmask in unrooted trees allows for the
    (arbitrarily-labeled) group '0' to be consistent across different
    representations, rotations, and orientations of trees.

    On the other hand, we use :attr:`Bipartition.leafset_bitmask` data to work
    with various ancestor-descendent relationships *within* the *same* tree:
    for example, to quickly assess if a taxon descends from a particular
    node in a given tree, or if a particular node is a common ancestor of
    two taxa in a given tree.

    The |Bipartition| object might be used in keys in dictionaries and
    look-up tables implemented as sets to allow for, e.g., calculation of
    support in terms of the number times a particular bipartition is observed.
    The :attr:`Bipartition.bitmask` is used as hash value for this purpose. As
    such, it is crucial that this value does not change once a particular
    |Bipartition| object is stored in a dictionary or set. To this end,
    we impose the constraint that |Bipartition| objects are immutable
    unless the ``is_mutable`` attribute is explicitly set to |True| as a sort
    of waiver signed by the client code. Client code does this at its risk,
    with the warning that anything up to and including the implosion of the
    universe may occur if the |Bipartition| object is a member of an set
    of dictionary at the time (or, at the very least, the modified
    |Bipartition| object may not be accessible from dictionaries
    and sets in which it is stored, or may occlude other
    |Bipartition| objects in the container).

    Note
    ----

    There are two possible ways of mapping taxa to bits in a bitarray or bitstring.

    In the "Least-Signficiant-Bit" (LSB) scheme, the first taxon corresponds to the
    least-significant, or left-most bit. So, given four taxa, indexed from 1 to 4,
    taxon 1 would map to 0b0001, taxon 2 would map to 0b0010, taxon 3 would map
    to 0b0100, and taxon 4 would map to 0b1000.

    In the "Most-Significant-Bit" (MSB) scheme, on the other hand, the first taxon
    corresponds to the most-significant, or right-most bit. So, given four
    taxa, indexed from 1 to 4, taxon 1 would map to 0b1000, taxon 2 would map
    to 0b0100, taxon 3 would map to 0b0010, and taxon 4 would map to 0b0001.

    We selected the Least Significant Bit (LSB) approach because the MSB scheme
    requires the size of the taxon namespace to fixed before the index can be
    assigned to any taxa. For example, under the MSB scheme, if there are 4
    taxa, the bitmask for taxon 1 is 0b1000 == 8, but if another taxon is
    added, then the bitmask for taxon 1 will become 0b10000 == 16. On the other
    hand, under the LSB scheme, the bitmask for taxon 1 will be 0b0001 == 1 if
    there are 4 taxa, and 0b00001 == 1 if there 5 taxa, and so on. This
    stability of taxon indexes even as the taxon namespace grows is a strongly
    desirable property, and this the adoption of the LSB scheme.

    Constraining the first taxon to be in group 0 (LSB-0) rather than group 1
    (LSB-1) is motivated by the fact that, in the former, we can would combine
    the bitmasks of child nodes using OR (logical addition) operations when
    calculating the bitmask for a parent node, whereas, with the latter, we
    would need to use AND operations. The former strikes us as more intuitive.

    """

    def normalize_bitmask(bitmask, fill_bitmask, lowest_relevant_bit):
        if bitmask & lowest_relevant_bit:
            return (~bitmask) & fill_bitmask             # force least-significant bit to 0
        else:
            return bitmask & fill_bitmask                # keep least-significant bit as 0
    normalize_bitmask = staticmethod(normalize_bitmask)

    def is_trivial_bitmask(bitmask, fill_bitmask):
        """
        Returns True if the bitmask occurs in any tree of the taxa ``mask`` -- if
        there is only fewer than two 1's or fewer than two 0's in ``bitmask`` (among
        all of the that are 1 in mask).
        """
        masked_split = bitmask & fill_bitmask
        if bitmask == 0 or bitmask == fill_bitmask:
            return True
        if ((masked_split - 1) & masked_split) == 0:
            return True
        cm = (~bitmask) & fill_bitmask
        if ((cm - 1) & cm) == 0:
            return True
        return False
    is_trivial_bitmask = staticmethod(is_trivial_bitmask)

    def is_trivial_leafset(leafset_bitmask):
        return bitprocessing.num_set_bits(leafset_bitmask) == 1
    is_trivial_leafset = staticmethod(is_trivial_leafset)

    def is_compatible_bitmasks(m1, m2, fill_bitmask):
        """
        Returns |True| if ``m1`` is compatible with ``m2``

        Parameters
        ----------
        m1 : int
            A bitmask representing a split.
        m2 : int
            A bitmask representing a split.

        Returns
        -------
        bool
            |True| if ``m1`` is compatible with ``m2``. |False| otherwise.
        """
        if fill_bitmask != 0:
            m1 = fill_bitmask & m1
            m2 = fill_bitmask & m2
        if 0 == (m1 & m2):
            return True
        c2 = m1 ^ m2
        if 0 == (m1 & c2):
            return True
        c1 = fill_bitmask ^ m1
        if 0 == (c1 & m2):
            return True
        if 0 == (c1 & c2):
            return True
        return False
    is_compatible_bitmasks = staticmethod(is_compatible_bitmasks)

    ##############################################################################
    ## Life-cycle

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------
        bitmask : integer
            A bit array representing the membership of taxa, with the
            least-significant bit corresponding to the first taxon, the next
            least-signficant bit correspodning to the second taxon, and so on,
            till the last taxon corresponding to the most-significant bit.
            Taxon membership in one of two arbitrary groups, '0' or '1', is
            indicated by its correspondign bit being unset or set,
            respectively.
        leafset_bitmask : integer
            A bit array representing the presence or absence of taxa in the
            subtree descending from the child node of the edge of which this
            bipartition is associated. The least-significant bit corresponds to
            the first taxon, the next least-signficant bit corresponds to the
            second taxon, and so on, with the last taxon corresponding to the
            most-significant bit.
        tree_leafset_bitmask : integer
            The ``leafset_bitmask`` of the root edge of the tree with which this
            bipartition is associated. In, general, this will be $0b1111...n$,
            where $n$ is the number of taxa, *except* in cases of trees with
            incomplete leaf-sets, where the positions corresponding to the
            missing taxa will have the bits unset.
        is_rooted : bool
            Specifies whether or not the tree with which this bipartition is
            associated is rooted.
        """
        self._split_bitmask = kwargs.get("bitmask", 0)
        self._leafset_bitmask = kwargs.get("leafset_bitmask", self._split_bitmask)
        self._tree_leafset_bitmask = kwargs.get("tree_leafset_bitmask", None)
        self._lowest_relevant_bit = None
        self._is_rooted = kwargs.get("is_rooted", None)
        # self.edge = kwargs.get("edge", None)
        is_mutable = kwargs.get("is_mutable", None)
        if kwargs.get("compile_bipartition", True):
            self.is_mutable = True
            self.compile_split_bitmask(
                    leafset_bitmask=self._leafset_bitmask,
                    tree_leafset_bitmask=self._tree_leafset_bitmask)
            if is_mutable is None:
                self.is_mutable = True
            else:
                self.is_mutable = is_mutable
        elif is_mutable is not None:
            self.is_mutable = is_mutable

    ##############################################################################
    ## Identity

    def __hash__(self):
        assert not self.is_mutable, "Bipartition is mutable: hash is unstable"
        return self._split_bitmask or 0

    def __eq__(self, other):
        # return self._split_bitmask == other._split_bitmask
        return (self._split_bitmask is not None and self._split_bitmask == other._split_bitmask) or (self._split_bitmask is other._split_bitmask)

    ##############################################################################
    ## All properties are publically read-only if not mutable

    def _get_split_bitmask(self):
        return self._split_bitmask
    def _set_split_bitmask(self, value):
        assert self.is_mutable, "Bipartition instance is not mutable"
        self._split_bitmask = value
    split_bitmask = property(_get_split_bitmask, _set_split_bitmask)

    def _get_leafset_bitmask(self):
        return self._leafset_bitmask
    def _set_leafset_bitmask(self, value):
        assert self.is_mutable, "Bipartition instance is not mutable"
        self._leafset_bitmask = value
    leafset_bitmask = property(_get_leafset_bitmask, _set_leafset_bitmask)

    def _get_tree_leafset_bitmask(self):
        return self._tree_leafset_bitmask
    def _set_tree_leafset_bitmask(self, value):
        assert self.is_mutable, "Bipartition instance is not mutable"
        self.compile_tree_leafset_bitmask(value)
    tree_leafset_bitmask = property(_get_tree_leafset_bitmask, _set_tree_leafset_bitmask)

    def _get_is_rooted(self):
        return self._is_rooted
    def _set_is_rooted(self, value):
        assert self.is_mutable, "Bipartition instance is not mutable"
        self._is_rooted = value
    is_rooted = property(_get_is_rooted, _set_is_rooted)

    ##############################################################################
    ## Representation

    def __str__(self):
        return bin(self._split_bitmask)[2:].rjust(bitprocessing.bit_length(self._tree_leafset_bitmask), '0')

    def __int__(self):
        return self._split_bitmask

    def split_as_int(self):
        return self._split_bitmask

    def leafset_as_int(self):
        return self._leafset_bitmask

    def split_as_bitstring(self, symbol0="0", symbol1="1", reverse=False):
        """
        Composes and returns and representation of the bipartition as a
        bitstring.

        Parameters
        ----------
        symbol1 : str
            The symbol to represent group '0' in the bitmask.
        symbol1 : str
            The symbol to represent group '1' in the bitmask.
        reverse : bool
            If |True|, then the first taxon will correspond to the
            most-significant bit, instead of the least-significant bit, as is
            the default.

        Returns
        -------
        str
            The bitstring representing the bipartition.

        Example
        -------
        To represent a bipartition in the same scheme used by, e.g. PAUP* or
        Mr. Bayes::

            print(bipartition.split_as_bitstring('.', '*', reverse=True))
        """
        return self.bitmask_as_bitstring(
                mask=self._split_bitmask,
                symbol0=symbol0,
                symbol1=symbol1,
                reverse=reverse)

    def leafset_as_bitstring(self, symbol0="0", symbol1="1", reverse=False):
        """
        Composes and returns and representation of the bipartition leafset as a
        bitstring.

        Parameters
        ----------
        symbol1 : str
            The symbol to represent group '0' in the bitmask.
        symbol1 : str
            The symbol to represent group '1' in the bitmask.
        reverse : bool
            If |True|, then the first taxon will correspond to the
            most-significant bit, instead of the least-significant bit, as is
            the default.

        Returns
        -------
        str
            The bitstring representing the bipartition.

        Example
        -------
        To represent a bipartition in the same scheme used by, e.g. PAUP* or
        Mr. Bayes::

            print(bipartition.leafset_as_bitstring('.', '*', reverse=True))
        """
        return self.bitmask_as_bitstring(
                mask=self._leafset_bitmask,
                symbol0=symbol0,
                symbol1=symbol1,
                reverse=reverse)

    def bitmask_as_bitstring(self, mask, symbol0=None, symbol1=None, reverse=False):
        return bitprocessing.int_as_bitstring(mask,
                length=bitprocessing.bit_length(self._tree_leafset_bitmask),
                symbol0=symbol0,
                symbol1=symbol1,
                reverse=reverse)

    ##############################################################################
    ## Calculation

    def compile_tree_leafset_bitmask(self,
            tree_leafset_bitmask,
            lowest_relevant_bit=None):
        """
        Avoids recalculation of ``lowest_relevant_bit`` if specified.
        """
        assert self.is_mutable, "Bipartition instance is not mutable"
        self._tree_leafset_bitmask = tree_leafset_bitmask
        if lowest_relevant_bit is not None:
            self._lowest_relevant_bit = lowest_relevant_bit
        elif self._tree_leafset_bitmask:
            self._lowest_relevant_bit = bitprocessing.least_significant_set_bit(self._tree_leafset_bitmask)
        else:
            self._lowest_relevant_bit = None
        return self._tree_leafset_bitmask

    def compile_leafset_bitmask(self,
           leafset_bitmask=None,
           tree_leafset_bitmask=None):
        assert self.is_mutable, "Bipartition instance is not mutable"
        if tree_leafset_bitmask is not None:
            self.compile_tree_leafset_bitmask(tree_leafset_bitmask)
        if leafset_bitmask is None:
            leafset_bitmask = self._leafset_bitmask
        if self._tree_leafset_bitmask:
            self._leafset_bitmask = leafset_bitmask & self._tree_leafset_bitmask
        else:
            self._leafset_bitmask = leafset_bitmask
        return self._leafset_bitmask

    def compile_split_bitmask(self,
           leafset_bitmask=None,
           tree_leafset_bitmask=None,
           is_rooted=None,
           is_mutable=True):
        """
        Updates the values of the various masks specified and calculates the
        normalized bipartition bitmask.

        If a rooted bipartition, then this is set to the value of the leafset
        bitmask.
        If an unrooted bipartition, then the leafset bitmask is normalized such that
        the lowest-significant bit (i.e., the group to which the first taxon
        belongs) is set to '0'.

        Also makes this bipartition immutable (unless ``is_mutable`` is |False|),
        which facilitates it being used in dictionaries and sets.

        Parameters
        ----------
        leafset_bitmask : integer
            A bit array representing the presence or absence of taxa in the
            subtree descending from the child node of the edge of which this
            bipartition is associated. The least-significant bit corresponds to
            the first taxon, the next least-signficant bit corresponds to the
            second taxon, and so on, with the last taxon corresponding to the
            most-significant bit. If not specified or |None|, the current value
            of ``self.leafset_bitmask`` is used.
        tree_leafset_bitmask : integer
            The ``leafset_bitmask`` of the root edge of the tree with which this
            bipartition is associated. In, general, this will be $0b1111...n$,
            where $n$ is the number of taxa, *except* in cases of trees with
            incomplete leaf-sets, where the positions corresponding to the
            missing taxa will have the bits unset. If not specified or |None|,
            the current value of ``self.tree_leafset_bitmask`` is used.
        is_rooted : bool
            Specifies whether or not the tree with which this bipartition is
            associated is rooted. If not specified or |None|, the current value
            of ``self.is_rooted`` is used.

        Returns
        -------
        integer
            The bipartition bitmask.
        """
        assert self.is_mutable, "Bipartition instance is not mutable"
        if is_rooted is not None:
            self._is_rooted = is_rooted
        if tree_leafset_bitmask:
            self.compile_tree_leafset_bitmask(tree_leafset_bitmask=tree_leafset_bitmask)
        if leafset_bitmask:
            self.compile_leafset_bitmask(leafset_bitmask=leafset_bitmask)
        if self._leafset_bitmask is None:
            return
        if self._tree_leafset_bitmask is None:
            return
        if self._is_rooted:
            self._split_bitmask = self._leafset_bitmask
        else:
            self._split_bitmask = Bipartition.normalize_bitmask(
                    bitmask=self._leafset_bitmask,
                    fill_bitmask=self._tree_leafset_bitmask,
                    lowest_relevant_bit=self._lowest_relevant_bit)
        if is_mutable is not None:
            self.is_mutable = is_mutable
        return self._split_bitmask

    def compile_bipartition(self, is_mutable=None):
        """
        Updates the values of the various masks specified and calculates the
        normalized bipartition bitmask.

        If a rooted bipartition, then this is set to the value of the leafset
        bitmask.
        If an unrooted bipartition, then the leafset bitmask is normalized such that
        the lowest-significant bit (i.e., the group to which the first taxon
        belongs) is set to '0'.

        Also makes this bipartition immutable (unless ``is_mutable`` is |False|),
        which facilitates it being used in dictionaries and sets.

        Note that this requires full population of the following fields:
            - self._leafset_bitmask
            - self._tree_leafset_bitmask
        """
        self.compile_split_bitmask(self,
            leafset_bitmask=self._leafset_bitmask,
            tree_leafset_bitmask=self._tree_leafset_bitmask,
            is_rooted=self._is_rooted,
            is_mutable=is_mutable)

    ##############################################################################
    ## Operations

    def normalize(self, bitmask, convention="lsb0"):
        """
        Return ``bitmask`` ensuring that the bit corresponding to the first
        taxon is 1.
        """
        if convention == "lsb0":
            if self._lowest_relevant_bit & bitmask:
                return (~bitmask) & self._tree_leafset_bitmask
            else:
                return bitmask & self._tree_leafset_bitmask
        elif convention == "lsb1":
            if self._lowest_relevant_bit & bitmask:
                return bitmask & self._tree_leafset_bitmask
            else:
                return (~bitmask) & self._tree_leafset_bitmask
        else:
            raise ValueError("Unrecognized convention: {}".format(convention))

    def is_compatible_with(self, other):
        """
        Returns |True| if ``other`` is compatible with self.

        Parameters
        ----------
        other : |Bipartition|
            The bipartition to check for compatibility.

        Returns
        -------
        bool
            |True| if ``other`` is compatible with ``self``; |False| otherwise.
        """
        m1 = self._split_bitmask
        if isinstance(other, int):
            m2 = other
        else:
            m2 = other._split_bitmask
        return Bipartition.is_compatible_bitmasks(m1, m2, self._tree_leafset_bitmask)

    def is_incompatible_with(self, other):
        """
        Returns |True| if ``other`` conflicts with self.

        Parameters
        ----------
        other : |Bipartition|
            The bipartition to check for conflicts.

        Returns
        -------
        bool
            |True| if ``other`` conflicts with ``self``; |False| otherwise.
        """
        return not self.is_compatible_with(other)

    def is_nested_within(self, other, is_other_masked_for_tree_leafset=False):
        """
        Returns |True| if the current bipartition is contained
        within other.

        Parameters
        ----------
        other : |Bipartition|
            The bipartition to check.

        Returns
        -------
        bool
            |True| if the the bipartition is "contained" within ``other``
        """
        if self._is_rooted:
            m1 = self._leafset_bitmask
            m2 = other._leafset_bitmask
        else:
            m1 = self._split_bitmask
            m2 = other._split_bitmask
        if not is_other_masked_for_tree_leafset:
            m2 = self._tree_leafset_bitmask & m2
        return ( (m1 & m2) == m1 )

    def is_leafset_nested_within(self, other):
        """
        Returns |True| if the leafset of ``self`` is a subset of the leafset of
        ``other``.

        Parameters
        ----------
        other : |Bipartition|
            The bipartition to check for compatibility.

        Returns
        -------
        bool
            |True| if the leafset of ``self`` is contained in ``other``.
        """
        if isinstance(other, int):
            m2 = other
        else:
            m2 = other._leafset_bitmask
        m2 = self._tree_leafset_bitmask & m2
        return ( (m2 & self._leafset_bitmask) ==  self._leafset_bitmask )

    def is_trivial(self):
        """
        Returns
        -------
        bool
            |True| if this bipartition divides a leaf and the rest of the
            tree.
        """
        return Bipartition.is_trivial_bitmask(self._split_bitmask,
                self._tree_leafset_bitmask)

    def split_as_newick_string(self,
            taxon_namespace,
            preserve_spaces=False,
            quote_underscores=True):
        """
        Represents this bipartition split as a newick string.

        Parameters
        ----------
        taxon_namespace : |TaxonNamespace| instance
            The operational taxonomic unit concept namespace to reference.
        preserve_spaces : boolean, optional
            If |False| (default), then spaces in taxon labels will be replaced
            by underscores. If |True|, then taxon labels with spaces will be
            wrapped in quotes.
        quote_underscores : boolean, optional
            If |True| (default), then taxon labels with underscores will be
            wrapped in quotes. If |False|, then the labels will not be wrapped
            in quotes.

        Returns
        -------
        string
            NEWICK representation of split specified by ``bitmask``.
        """
        return taxon_namespace.bitmask_as_newick_string(
                bitmask=self._split_bitmask,
                preserve_spaces=preserve_spaces,
                quote_underscores=quote_underscores)

    def leafset_as_newick_string(self,
            taxon_namespace,
            preserve_spaces=False,
            quote_underscores=True):
        """
        Represents this bipartition leafset as a newick string.

        Parameters
        ----------
        taxon_namespace : |TaxonNamespace| instance
            The operational taxonomic unit concept namespace to reference.
        preserve_spaces : boolean, optional
            If |False| (default), then spaces in taxon labels will be replaced
            by underscores. If |True|, then taxon labels with spaces will be
            wrapped in quotes.
        quote_underscores : boolean, optional
            If |True| (default), then taxon labels with underscores will be
            wrapped in quotes. If |False|, then the labels will not be wrapped
            in quotes.

        Returns
        -------
        string
            NEWICK representation of split specified by ``bitmask``.
        """
        return taxon_namespace.bitmask_as_newick_string(
                bitmask=self._leafset_bitmask,
                preserve_spaces=preserve_spaces,
                quote_underscores=quote_underscores)

    def leafset_taxa(self, taxon_namespace, index=0):
        """
        Returns list of |Taxon| objects in the leafset of this
        bipartition.

        Parameters
        ----------
        taxon_namespace : |TaxonNamespace| instance
            The operational taxonomic unit concept namespace to reference.
        index : integer, optional
            Start from this |Taxon| object instead of the first
            |Taxon| object in the collection.

        Returns
        -------
        :py:class:`list` [|Taxon|]
            List of |Taxon| objects specified or spanned by
            ``bitmask``.
        """
        return taxon_namespace.bitmask_taxa_list(
                bitmask=self._leafset_bitmask,
                index=index)

    # def as_newick_string
    # def is_trivial
    # def is_non_singleton
    # def leafset_hash
    # def leafset_as_bitstring
    # def is_compatible

##############################################################################
### Edge

class Edge(
        basemodel.DataObject,
        basemodel.Annotable):
    """
    An :term:``edge`` on a :term:``tree``.
    """

    ###########################################################################
    ### Life-cycle and Identity

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------
        head_node : |Node|, optional
            Node from to which this edge links, i.e., the child node of this
            node ``tail_node``.
        length : numerical, optional
            A value representing the weight of the edge.
        rootedge : boolean, optional
            Is the child node of this edge the root or seed node of the tree?
        label : string, optional
            Label for this edge.

        """
        basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
        self._head_node = kwargs.pop("head_node", None)
        if "tail_node" in kwargs:
            raise TypeError("Setting the tail node directly is no longer supported: instead, set the parent node of the head node")
        self.rootedge = kwargs.pop("rootedge", None)
        self.length = kwargs.pop("length", None)
        if kwargs:
            raise TypeError("Unsupported keyword arguments: {}".format(kwargs))

        self._bipartition = None
        self.comments = []

    def __copy__(self, memo=None):
        raise TypeError("Cannot directly copy Edge")

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot directly copy Edge")

    def __deepcopy__(self, memo=None):
        # call Annotable.__deepcopy__()
        return basemodel.Annotable.__deepcopy__(self, memo=memo)
        # return super(Edge, self).__deepcopy__(memo=memo)

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    ###########################################################################
    ### Basic Structure

    def _get_tail_node(self):
        if self._head_node is None:
            return None
        return self._head_node._parent_node
    def _set_tail_node(self, node):
        if self._head_node is None:
            raise ValueError("'_head_node' is 'None': cannot assign 'tail_node'")
        # Go through managed property instead of
        # setting attribute to ensure book-keeping
        self._head_node.parent_node = node
    tail_node = property(_get_tail_node, _set_tail_node)

    def _get_head_node(self):
        return self._head_node
    def _set_head_node(self, node):
        # Go through managed property instead of setting attribute to ensure
        # book-keeping; following should also set ``_head_node`` of ``self``
        node.edge = self
    head_node = property(_get_head_node, _set_head_node)

    def is_leaf(self):
        "Returns True if the head node has no children"
        return self.head_node and self.head_node.is_leaf()

    def is_terminal(self):
        return self.is_leaf()

    def is_internal(self):
        "Returns True if the head node has children"
        return self.head_node and not self.head_node.is_leaf()

    def get_adjacent_edges(self):
        """
        Returns a list of all edges that "share" a node with ``self``.
        """
        he = [i for i in self.head_node.incident_edges() if i is not self]
        te = [i for i in self.tail_node.incident_edges() if i is not self]
        he.extend(te)
        return he
    adjacent_edges = property(get_adjacent_edges)

    ###########################################################################
    ### Structural Manipulation

    def collapse(self, adjust_collapsed_head_children_edge_lengths=False):
        """
        Inserts all children of the head_node of self as children of the
        tail_node of self in the same place in the child_node list that
        head_node had occupied. The edge length and head_node will no longer be
        part of the tree unless ``adjust_collapsed_head_children_edge_lengths``.
        is True.
        """
        to_del = self.head_node
        parent = self.tail_node
        if not parent:
            return
        children = to_del.child_nodes()
        if not children:
            raise ValueError('collapse_self called with a terminal.')
        pos = parent.child_nodes().index(to_del)
        parent.remove_child(to_del)
        for child in children:
            parent.insert_child(pos, child)
            pos += 1
            if adjust_collapsed_head_children_edge_lengths and self.length is not None:
                # print id(child), child.edge.length, self.length
                if child.edge.length is None:
                    child.edge.length = self.length
                else:
                    child.edge.length += self.length

    def invert(self, update_bipartitions=False):
        """
        Changes polarity of edge.
        """
        # self.head_node, self.tail_node = self.tail_node, self.head_node

        if not self.head_node:
            raise ValueError("Cannot invert edge with 'None' for head node")
        if not self.tail_node:
            raise ValueError("Cannot invert edge with 'None' for tail node")

        old_head_node = self.head_node
        new_tail_node = old_head_node
        old_tail_node = self.tail_node
        new_head_node = old_tail_node
        grandparent = old_tail_node._parent_node
        if grandparent is not None:
            for idx, ch in enumerate(grandparent._child_nodes):
                if ch is old_tail_node:
                    grandparent._child_nodes[idx] = old_head_node
                    break
            else:
                # we did not break loop: force insertion of old_head_node if
                # not already there
                if old_head_node not in grandparent._child_nodes:
                    grandparent._child_nodes.append(old_head_node)
        assert old_head_node in old_tail_node._child_nodes
        old_tail_node.remove_child(old_head_node)
        assert old_head_node not in old_tail_node._child_nodes
        old_head_node.add_child(old_tail_node)
        old_tail_node.edge.length, old_head_node.edge.length = old_head_node.edge.length, old_tail_node.edge_length

    ###########################################################################
    ### Bipartition Management

    def _get_bipartition(self):
        if self._bipartition is None:
            self._bipartition = Bipartition(
                    edge=self,
                    is_mutable=True,
                    )
        return self._bipartition
    def _set_bipartition(self, v=None):
        self._bipartition = v
    bipartition = property(_get_bipartition, _set_bipartition)

    def _get_split_bitmask(self):
        return self.bipartition._split_bitmask
    def _set_split_bitmask(self, h):
        self.bipartition._split_bitmask = h
    split_bitmask = property(_get_split_bitmask, _set_split_bitmask)

    def _get_leafset_bitmask(self):
        return self.bipartition._leafset_bitmask
    def _set_leafset_bitmask(self, h):
        self.bipartition._leafset_bitmask = h
    leafset_bitmask = property(_get_leafset_bitmask, _set_leafset_bitmask)

    def _get_tree_leafset_bitmask(self):
        return self.bipartition._tree_leafset_bitmask
    def _set_tree_leafset_bitmask(self, h):
        self.bipartition._tree_leafset_bitmask = h
    tree_leafset_bitmask = property(_get_tree_leafset_bitmask, _set_tree_leafset_bitmask)

    def split_as_bitstring(self):
        return self.bipartition.split_as_bitstring()

    def leafset_as_bitstring(self):
        return self.bipartition.leafset_as_bitstring()

    ###########################################################################
    ### Representation

    def description(self,
            depth=1,
            indent=0,
            itemize="",
            output=None,
            taxon_namespace=None):
        """
        Returns description of object, up to level ``depth``.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s, Length=%s)" % (id(self), str(self.length))
        else:
            label = " (%s: '%s', Length=%s)" % (id(self), self.label, str(self.length))
        output_strio.write('%s%sEdge object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            leader1 = ' ' * (indent + 4)
            leader2 = ' ' * (indent + 8)
            output_strio.write('\n%s[Length]' % leader1)
            if self.length is not None:
                length = self.length
            else:
                length = "None"
            output_strio.write('\n%s%s' % (leader2, length))
            output_strio.write('\n%s[Tail Node]' % leader1)
            if self.tail_node is not None:
                tn = self.tail_node.description(0)
            else:
                tn = "None"
            output_strio.write('\n%s%s' % (leader2, tn))
            output_strio.write('\n%s[Head Node]' % leader1)
            if self.head_node is not None:
                hn = self.head_node.description(0)
            else:
                hn = "None"
            output_strio.write('\n%s%s' % (leader2, hn))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

##############################################################################
### Node

class Node(
        basemodel.DataObject,
        basemodel.Annotable):
    """
    A :term:|Node| on a :term:|Tree|.
    """

    ###########################################################################
    ### Life-cycle

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------
        taxon : |Taxon|, optional
            The |Taxon| instance representing the operational taxonomic
            unit concept associated with this Node.
        label : string, optional
            A label for this node.
        edge_length : numeric, optional
            Length or weight of the edge subtending this node.

        """
        basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
        self.taxon = kwargs.pop("taxon", None)
        self.age = None
        self._edge = None
        self._child_nodes = []
        self._parent_node = None
        self.edge = Edge(head_node=self,
                length=kwargs.pop("edge_length", None))
        if kwargs:
            raise TypeError("Unsupported keyword arguments: {}".format(kwargs))
        self.comments = []

    def __copy__(self, memo=None):
        raise TypeError("Cannot directly copy Edge")

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot directly copy Node")

    def __deepcopy__(self, memo=None):
        return basemodel.Annotable.__deepcopy__(self, memo=memo)
        # if memo is None:
        #     memo = {}
        # other = basemodel.Annotable.__deepcopy__(self, memo=memo)
        # memo[id(self._child_nodes)] = other._child_nodes
        # for ch in self._child_nodes:
        #     try:
        #         och = memo[id(ch)]
        #         if och not in other._child_nodes:
        #             other._child_nodes.append(och)
        #     except KeyError:
        #         och = copy.deepcopy(ch, memo)
        #         memo[id(chd)] = och
        #         if och not in other._child_nodes:
        #             other._child_nodes.append(och)
        # return other
        # return super(Node, self).__deepcopy__(memo=memo)

    ###########################################################################
    ### Identity

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        # IMPORTANT LESSON LEARNED: if you define __hash__, you *must* define __eq__
        return self is other

    def __repr__(self):
        return "<{} object at {}: '{}' ({})>".format(self.__class__.__name__, hex(id(self)), self._label, repr(self.taxon))

    ###########################################################################
    ### Iterators

    def preorder_iter(self, filter_fn=None):
        """
        Pre-order iterator over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited before its
        children. Nodes can optionally be filtered by ``filter_fn``: only nodes
        for which ``filter_fn`` returns |True| when called with the node as an
        argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes of the subtree rooted at this node in
            pre-order sequence.
        """
        stack = [self]
        while stack:
            node = stack.pop()
            if filter_fn is None or filter_fn(node):
                yield node
            stack.extend(n for n in reversed(node._child_nodes))

    def preorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes of subtree rooted at this node.

        Visits self and all internal descendant nodes, with each node visited
        before its children. In DendroPy, "internal nodes" are nodes that have
        at least one child node, and thus the root or seed node is typically included
        unless ``exclude_seed_node`` is |True|. Nodes can optionally be filtered
        by ``filter_fn``: only nodes for which ``filter_fn`` returns |True| when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If |False| (default), then the seed node or root is visited. If
            |True|, then the seed node is skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the internal nodes of the subtree rooted at
            this node in pre-order sequence.
        """
        if exclude_seed_node:
            froot = lambda x: x._parent_node is not None
        else:
            froot = lambda x: True
        if filter_fn:
            f = lambda x: (froot(x) and x._child_nodes and filter_fn(x)) or None
        else:
            f = lambda x: (x and froot(x) and x._child_nodes) or None
        return self.preorder_iter(filter_fn=f)

    def postorder_iter(self, filter_fn=None):
        """
        Post-order iterator over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited after its
        children. Nodes can optionally be filtered by ``filter_fn``: only nodes
        for which ``filter_fn`` returns |True| when called with the node as an
        argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the nodes of the subtree rooted at
            this node in post-order sequence.
        """
        # if self._child_nodes:
        #     for nd in self._child_nodes:
        #         for ch in nd.postorder_iter(filter_fn=filter_fn):
        #             yield ch
        # if filter_fn is None or filter_fn(self):
        #     yield self
        # return

        # stack = [(self, False)]
        # while stack:
        #     node, state = stack.pop(0)
        #     if state:
        #         if filter_fn is None or filter_fn(node):
        #             yield node
        #     else:
        #         stack.insert(0, (node, True))
        #         child_nodes = [(n, False) for n in node._child_nodes]
        #         child_nodes.extend(stack)
        #         stack = child_nodes

        ## Prefer `pop()` to `pop(0)`.
        ## Thanks to Mark T. Holder
        ## From peyotl commits: d1ffef2 + 19fdea1
        stack = [(self, False)]
        while stack:
            node, state = stack.pop()
            if state:
                if filter_fn is None or filter_fn(node):
                    yield node
            else:
                stack.append((node, True))
                stack.extend([(n, False) for n in reversed(node._child_nodes)])

    def postorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes of subtree rooted at this node.

        Visits self and all internal descendant nodes, with each node visited
        after its children. In DendroPy, "internal nodes" are nodes that have
        at least one child node, and thus the root or seed node is typically
        included unless ``exclude_seed_node`` is |True|. Nodes can optionally be
        filtered by ``filter_fn``: only nodes for which ``filter_fn`` returns
        |True| when passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If |False| (default), then the seed node or root is visited. If
            |True|, then the seed node is skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the internal nodes of the subtree rooted at
            this node in post-order sequence.
        """
        if exclude_seed_node:
            froot = lambda x: x._parent_node is not None
        else:
            froot = lambda x: True
        if filter_fn:
            f = lambda x: (froot(x) and x._child_nodes and filter_fn(x)) or None
        else:
            f = lambda x: (x and froot(x) and x._child_nodes) or None
        return self.postorder_iter(filter_fn=f)

    def levelorder_iter(self, filter_fn=None):
        """
        Level-order iteration over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node and other nodes at
        the same level (distance from root) visited before their children.
        Nodes can optionally be filtered by ``filter_fn``: only nodes for which
        ``filter_fn`` returns |True| when called with the node as an argument are
        visited.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes of the subtree rooted at this node in
            level-order sequence.
        """
        if filter_fn is None or filter_fn(self):
            yield self
        remaining = self.child_nodes()
        while len(remaining) > 0:
            node = remaining.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = node.child_nodes()
            remaining.extend(child_nodes)

    def level_order_iter(self, filter_fn=None):
        """
        DEPRECATED: Use :meth:`Node.levelorder_iter()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'level_order_iter()' will no longer be supported in future releases; use 'levelorder_iter()' instead",
                stacklevel=3)
        return self.levelorder_iter(filter_fn=filter_fn)

    def inorder_iter(self, filter_fn=None):
        """
        In-order iteration over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited in-between
        its children. Only valid for strictly-bifurcating trees. Nodes can
        optionally be filtered by ``filter_fn``: only nodes for which ``filter_fn``
        returns |True| when called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes of the subtree rooted at this node in
            infix or in-order sequence.
        """
        if len(self._child_nodes) == 0:
            if filter_fn is None or filter_fn(self):
                yield self
        elif len(self._child_nodes) == 2:
            for nd in self._child_nodes[0].inorder_iter(filter_fn=filter_fn):
                yield nd
            if filter_fn is None or filter_fn(self):
                yield self
            for nd in self._child_nodes[1].inorder_iter(filter_fn=filter_fn):
                yield nd
        else:
            raise TypeError("In-order traversal only supported for binary trees")

    def leaf_iter(self, filter_fn=None):
        """
        Iterate over all tips or leaves that ultimately descend from this node.

        Visits all leaf or tip nodes descended from this node. Nodes can
        optionally be filtered by ``filter_fn``: only nodes for which ``filter_fn``
        returns |True| when called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding leaf nodes of the subtree rooted at this node.
        """
        if filter_fn:
            ff = lambda x: x.is_leaf() and filter_fn(x) or None
        else:
            ff = lambda x: x.is_leaf() and x or None
        for node in self.postorder_iter(ff):
            yield node

    def child_node_iter(self, filter_fn=None):
        """
        Iterator over all nodes that are the (immediate) children of this node.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes that have this node as a parent.
        """
        for node in self._child_nodes:
            if filter_fn is None or filter_fn(node):
                yield node

    def child_edge_iter(self, filter_fn=None):
        """
        Iterator over all edges that are the (immediate) children of this edge.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all edges visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Edge|]
            An iterator yielding edges that have this edge as a parent.
        """
        for node in self._child_nodes:
            if filter_fn is None or filter_fn(node.edge):
                yield node.edge

    def ancestor_iter(self, filter_fn=None, inclusive=False):
        """
        Iterator over all ancestors of this node.

        Visits all nodes that are the ancestors of this node.  If ``inclusive``
        is |True|, ``self`` is returned as the first item of the sequence;
        otherwise ``self`` is skipped. Nodes can optionally be filtered by
        ``filter_fn``: only nodes for which ``filter_fn`` returns |True| when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        inclusive : boolean, optional
            If |True|, includes this node in the sequence. If |False|, this is
            skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            Iterator over all predecessor/ancestor nodes of this node.
        """
        if inclusive and (filter_fn is None or filter_fn(self)):
            yield self
        node = self
        while node is not None:
            node = node._parent_node
            if node is not None \
                   and (filter_fn is None or filter_fn(node)):
                yield node

    def ageorder_iter(self, filter_fn=None, include_leaves=True, descending=False):
        """
        Iterator over nodes of subtree rooted at this node in order of the age
        of the node (i.e., the time since the present).

        Iterates over nodes in order of age ('age' is as given by the ``age``
        attribute, which is usually the sum of edge lengths from tips
        to node, i.e., time since present).
        If ``include_leaves`` is |True| (default), leaves are included in the
        iteration; if ``include_leaves`` is |False|, leaves will be skipped.
        If ``descending`` is |False| (default), younger nodes will be returned
        before older ones; if |True|, older nodes will be returned before
        younger ones.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (defau
        include_leaves : boolean, optional
            If |True| (default), then leaf nodes are included in the iteration.
            If |False|, then leaf nodes are skipped.lt), then all nodes visited will be yielded.
        descending : boolean, optional
            If |False| (default), then younger nodes are visited before older
            ones. If |True|, then older nodes are visited before younger ones.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            Iterator over age-ordered sequence of nodes in subtree rooted at
            this node.
        """
        # if not descending:
        #     leaves = [nd for nd in self.leaf_iter()]
        #     queued_pairs = []
        #     in_queue = set()
        #     for leaf in leaves:
        #         age_nd_tuple = (leaf.age, leaf)
        #         queued_pairs.insert(bisect.bisect(queued_pairs, age_nd_tuple), age_nd_tuple)
        #         in_queue.add(leaf)
        #     while queued_pairs:
        #         next_el = queued_pairs.pop(0)
        #         age, nd = next_el
        #         in_queue.remove(nd)
        #         p = nd._parent_node
        #         if p and p not in in_queue:
        #             age_nd_tuple = (p.age, p)
        #             queued_pairs.insert(bisect.bisect(queued_pairs, age_nd_tuple), age_nd_tuple)
        #             in_queue.add(p)
        #         if include_leaves or nd.is_internal():
        #             yield nd
        # else:
        #     nds = [(nd.age, nd) for nd in self.preorder_iter()]
        #     nds.sort(reverse=True)
        #     for nd in nds:
        #         if include_leaves or nd[1].is_internal():
        #             yield nd[1]
        nds = [nd for nd in self.preorder_iter()]
        if descending:
            reverse = True
        else:
            reverse = False
        nds.sort(key=lambda x: x.age, reverse=reverse)
        for nd in nds:
            if (include_leaves or nd._child_nodes) and (filter_fn is None or filter_fn(nd)):
                yield nd

    def age_order_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Deprecated: use :meth:`Node.ageorder_iter()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'age_order_iter()' will no longer be supported in future releases; use 'ageorder_iter()' instead",
                stacklevel=3)
        return self.ageorder_iter(include_leaves=include_leaves,
                filter_fn=filter_fn,
                descending=descending)

    ###########################################################################
    ### Node Processesor

    def apply(self, before_fn=None, after_fn=None, leaf_fn=None):
        """
        Applies function ``before_fn`` and ``after_fn`` to all internal nodes and
        ``leaf_fn`` to all terminal nodes in subtree starting with ``self``, with
        nodes visited in pre-order.

        Given a tree with preorder sequence of nodes of
        [a,b,i,e,j,k,c,g,l,m,f,n,h,o,p,]::

                           a
                          / \
                         /   \
                        /     \
                       /       \
                      /         \
                     /           \
                    /             c
                   b             / \
                  / \           /   \
                 /   e         /     f
                /   / \       /     / \
               /   /   \     g     /   h
              /   /     \   / \   /   / \
             i   j       k l   m n   o   p


        the following order of function calls results:

            before_fn(a)
            before_fn(b)
            leaf_fn(i)
            before_fn(e)
            leaf_fn(j)
            leaf_fn(k)
            after_fn(e)
            after_fn(b)
            before_fn(c)
            before_fn(g)
            leaf_fn(l)
            leaf_fn(m)
            after_fn(g)
            before_fn(f)
            leaf_fn(n)
            before_fn(h)
            leaf_fn(o)
            leaf_fn(p)
            after_fn(h)
            after_fn(f)
            after_fn(c)
            after_fn(a)

        Parameters
        ----------
        before_fn : function object or |None|
            A function object that takes a |Node| as its argument.
        after_fn : function object or |None|
            A function object that takes a |Node| as its argument.
        leaf_fn : function object or |None|
            A function object that takes a |Node| as its argument.

        Notes
        -----
        Adapted from work by Mark T. Holder (the ``peyotl`` module of the Open
        Tree of Life Project):

            https://github.com/OpenTreeOfLife/peyotl.git

        """
        stack = [self]
        while stack:
            node = stack.pop()
            if not node._child_nodes:
                if leaf_fn:
                    leaf_fn(node)
                # (while node is the last child of parent ...)
                while (
                        (node._parent_node is None)
                        or (node._parent_node._child_nodes[-1] is node)
                      ):
                    node = node._parent_node
                    if node is not None:
                        if after_fn is not None:
                            after_fn(node)
                    else:
                        break
            else:
                if before_fn is not None:
                    before_fn(node)
                stack.extend([i for i in reversed(node._child_nodes)])
        return

    ###########################################################################
    ### Child Node Access and Manipulation

    def set_child_nodes(self, child_nodes):
        """
        Assigns the set of child nodes for this node.

        Results in the ``parent_node`` attribute of each |Node| in ``nodes``
        as well as the ``tail_node`` attribute of corresponding |Edge|
        objects being assigned to ``self``.

        Parameters
        ----------
        child_nodes : collections.Iterable[|Node|]
            The (iterable) collection of child nodes to be assigned this node
            as a parent.
        """
        self.clear_child_nodes()
        # Go through add to ensure book-keeping
        # (e.g. avoiding multiple adds) takes
        # place.
        for nd in child_nodes:
            self.add_child(nd)

    def set_children(self, child_nodes):
        """Deprecated: use :meth:`Node.set_child_nodes()` instead."""
        return self.set_child_nodes(child_nodes)

    def add_child(self, node):
        """
        Adds a child node to this node if it is not already a child.

        Results in the ``parent_node`` attribute of ``node`` as well as the
        ``tail_node`` attribute of ``node.edge`` being assigned to ``self``.

        Parameters
        ----------
        node : |Node|
            The node to be added as a child of this node.

        Returns
        -------
        |Node|
            The node that was added.
        """
        assert node is not self, "Cannot add node as child of itself"
        assert self._parent_node is not node, "Cannot add a node's parent as its child: remove the node from its parent's child set first"
        node._parent_node = self
        if node not in self._child_nodes:
            self._child_nodes.append(node)
        return node

    def insert_child(self, index, node):
        """
        Adds a child node to this node.

        If the node is already a child of this node, then it is moved
        to the specified position.
        Results in the ``parent_node`` attribute of ``node`` as well as the
        ``tail_node`` attribute of ``node.edge`` being assigned to ``self``.

        Parameters
        ----------
        index : integer
            The index before which to insert the new node.
        node : |Node|
            The node to be added as a child of this node.

        Returns
        -------
        |Node|
            The node that was added.
        """
        node._parent_node = self
        try:
            cur_index = self._child_nodes.index(node)
        except ValueError:
            pass
        else:
            if cur_index == index:
                return
            self._child_nodes.remove(node)
        self._child_nodes.insert(index, node)
        return node

    def new_child(self, **kwargs):
        """
        Create and add a new child to this node.

        Parameters
        ----------
        \*\*kwargs : keyword arguments
            Keyword arguments will be passed directly to the |Node|
            constructor (:meth:`Node.__init()__`).

        Returns
        -------
        |Node|
            The new child node that was created and added.
        """
        node = self.__class__(**kwargs)
        return self.add_child(node=node)

    def insert_new_child(self, index, **kwargs):
        """
        Create and add a new child to this node at a particular position.

        Results in the ``parent_node`` attribute of ``node`` as well as the
        ``tail_node`` attribute of ``node.edge`` being assigned to ``self``.

        Parameters
        ----------
        index : integer
            The index before which to insert the new node.
        \*\*kwargs : keyword arguments, optional
            Keyword arguments will be passed directly to the |Node|
            constructor (:meth:`Node.__init()__`).

        Returns
        -------
        |Node|
            The new child node that was created and added.
        """
        node = self.__class__(**kwargs)
        return self.insert_child(index=index, node=node)

    def remove_child(self, node, suppress_unifurcations=False):
        """
        Removes a node from the child set of this node.

        Results in the parent of the node being removed set to |None|.  If
        ``suppress_unifurcations`` is |True|, if this node ends up having only one
        child after removal of the specified node, then this node will be
        removed from the tree, with its single child added to the child node
        set of its parent and the edge length adjusted accordingly.
        ``suppress_unifurcations`` should only be |True| for unrooted trees.

        Parameters
        ----------
        node : |Node|
            The node to be removed.
        suppress_unifurcations : boolean, optional
            If |False| (default), no action is taken. If |True|, then if the
            node removal results in a node with degree of two (i.e., a single
            parent and a single child), then it will be removed from
            the tree and its (sole) child will be added as a child of its
            parent (with edge lengths adjusted accordingly).

        Returns
        -------
        |Node|
            The node removed.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        children = self._child_nodes
        if node in children:
            node._parent_node = None
            node.edge.tail_node = None
            index = children.index(node)
            children.remove(node)
            if suppress_unifurcations:
                if self._parent_node:
                    if len(children) == 1:
                        child = children[0]
                        pos = self._parent_node._child_nodes.index(self)
                        self._parent_node.insert_child(pos, child)
                        self._parent_node.remove_child(self, suppress_unifurcations=False)
                        try:
                            child.edge.length += self.edge.length
                        except:
                            pass
                        self._child_nodes = []
                else:
                    to_remove = None
                    if len(children) == 2:
                        if children[0].is_internal():
                            to_remove = children[0]
                            other = children[1]
                        elif children[1].is_internal():
                            to_remove = children[1]
                            other = children[0]
                    if to_remove is not None:
                        try:
                            other.edge.length += to_remove.edge.length
                        except:
                            pass
                        pos = self._child_nodes.index(to_remove)
                        self.remove_child(to_remove, suppress_unifurcations=False)
                        tr_children = to_remove._child_nodes
                        tr_children.reverse()
                        for c in tr_children:
                            self.insert_child(pos, c)
                        to_remove._child_nodes = []
        else:
            raise ValueError("Tried to remove a node that is not listed as a child")
        return node

    def clear_child_nodes(self):
        """
        Removes all child nodes.
        """
        del self._child_nodes[:] # list.clear() is not in Python 2.7

    def reversible_remove_child(self, node, suppress_unifurcations=False):
        """
        This function is a (less-efficient) version of remove_child that also
        returns the data needed by reinsert_nodes to "undo" the removal.

        Returns a list of tuples.  The first element of each tuple is the
        node removed, the other elements are the information needed by
        ``reinsert_nodes`` in order to restore the tree to the same topology as
        it was before the call to ``remove_child.`` If ``suppress_unifurcations`` is False
        then the returned list will contain only one item.

        ``suppress_unifurcations`` should only be called on unrooted trees.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        children = self._child_nodes
        try:
            pos = children.index(node)
        except:
            raise ValueError("Tried to remove a node that is not listed as a child")
        removed = [(node, self, pos, [], None)]
        node._parent_node = None
        node.edge.tail_node = None
        children.remove(node)
        if suppress_unifurcations:
            p = self._parent_node
            if p:
                if len(children) == 1:
                    child = children[0]
                    pos = p._child_nodes.index(self)
                    p.insert_child(pos, child)
                    self._child_nodes = []
                    p.remove_child(self, suppress_unifurcations=False)
                    e = child.edge
                    try:
                        e.length += self.edge.length
                    except:
                        e = None
                    t = (self, p, pos, [child], e)
                    removed.append(t)
            else:
                to_remove = None
                if len(children) == 2:
                    if children[0].is_internal():
                        to_remove = children[0]
                        other = children[1]
                    elif children[1].is_internal():
                        to_remove = children[1]
                        other = children[0]
                if to_remove is not None:
                    e = other.edge
                    try:
                        e.length += to_remove.edge.length
                    except:
                        e = None
                    pos = self._child_nodes.index(to_remove)
                    self.remove_child(to_remove, suppress_unifurcations=False)
                    tr_children = to_remove._child_nodes
                    to_remove._child_nodes = []
                    for n, c in enumerate(tr_children):
                        new_pos = pos + n
                        self.insert_child(pos, c)
                    t = (to_remove, self, pos, tr_children, e)
                    removed.append(t)

        return removed

    def reinsert_nodes(self, nd_connection_list):
        """
        This function should be used to "undo" the effects of
        Node.reversible_remove_child NOTE: the behavior is only
        guaranteed if the tree has not been modified between the
        remove_child and reinsert_nodes calls! (or the tree has been
        restored such that the node/edge identities are identical to the
        state before the remove_child call.

        The order of info in each tuple is:

            0 - node removed
            1 - parent of node removed
            2 - pos in parent matrix
            3 - children of node removed that were "stolen"
            4 - edge that was lengthened by "stealing" length from node's edge
        """
        # we unroll the stack of operations
        for blob in nd_connection_list[-1::-1]:
            #_LOG.debug(blob)
            n, p, pos, children, e = blob
            for c in children:
                cp = c._parent_node
                if cp:
                    cp.remove_child(c)
                n.add_child(c)
            p.insert_child(pos, n)
            if e is not None:
                e.length -= n.edge.length

    def collapse_neighborhood(self, dist):
        if dist < 1:
            return
        children = self.child_nodes()
        for ch in children:
            if not ch.is_leaf():
                ch.edge.collapse()
        if self._parent_node:
            p = self._parent_node
            self.edge.collapse()
            p.collapse_neighborhood(dist -1)
        else:
            self.collapse_neighborhood(dist - 1)

    def collapse_clade(self):
        """Collapses all internal edges that are descendants of self."""
        if self.is_leaf():
            return
        leaves = [i for i in self.leaf_iter()]
        self.set_child_nodes(leaves)

    def collapse_conflicting(self, bipartition):
        """
        Collapses every edge in the subtree that conflicts with the given
        bipartition. This can include the edge subtending subtree_root.
        """
        to_collapse_head_nodes = []
        for nd in self.postorder_iter():
            if nd._child_nodes and nd.edge.bipartition.is_incompatible_with(bipartition):
                to_collapse_head_nodes.append(nd)
        for nd in to_collapse_head_nodes:
            e = nd.edge
            e.collapse()

    ###########################################################################
    ### Edge Access and Manipulation

    def _get_edge(self):
        """
        Returns the edge subtending this node.
        """
        return self._edge
    def _set_edge(self, new_edge):
        """
        Sets the edge subtending this node, and sets head_node of
        ``edge`` to point to self.
        """
        # if edge is None:
        #     raise ValueError("A Node cannot have 'None' for an edge")
        if new_edge is self._edge:
            return
        if self._parent_node is not None:
            try:
                self._parent_node._child_nodes.remove(self)
            except ValueError:
                pass

        ## Minimal management
        self._edge = new_edge
        if self._edge:
            self._edge._head_node = self

    edge = property(_get_edge, _set_edge)

    def _get_edge_length(self):
        """
        Returns the length of the edge subtending this node.
        """
        return self._edge.length
    def _set_edge_length(self, v=None):
        """
        Sets the edge subtending this node, and sets head_node of
        ``edge`` to point to self.
        """
        self._edge.length = v
    edge_length = property(_get_edge_length, _set_edge_length)

    def _get_bipartition(self):
        """
        Returns the bipartition for the edge subtending this node.
        """
        return self._edge.bipartition
    def _set_bipartition(self, v=None):
        """
        Sets the bipartition for the edge subtending this node.
        """
        self._edge.bipartition = v
    bipartition = property(_get_bipartition, _set_bipartition)

    def _get_split_bitmask(self):
        return self._edge.bipartition._split_bitmask
    def _set_split_bitmask(self, h):
        self._edge.bipartition._split_bitmask = h
    split_bitmask = property(_get_split_bitmask, _set_split_bitmask)

    def _get_leafset_bitmask(self):
        return self._edge.bipartition._leafset_bitmask
    def _set_leafset_bitmask(self, h):
        self._edge.bipartition._leafset_bitmask = h
    leafset_bitmask = property(_get_leafset_bitmask, _set_leafset_bitmask)

    def _get_tree_leafset_bitmask(self):
        return self._edge.bipartition._tree_leafset_bitmask
    def _set_tree_leafset_bitmask(self, h):
        self._edge.bipartition._tree_leafset_bitmask = h
    tree_leafset_bitmask = property(_get_tree_leafset_bitmask, _set_tree_leafset_bitmask)

    def split_as_bitstring(self):
        return self._edge.bipartition.split_as_bitstring()

    def leafset_as_bitstring(self):
        return self._edge.bipartition.leafset_as_bitstring()

    ###########################################################################
    ### Parent Access and Manipulation

    def _get_parent_node(self):
        """Returns the parent node of this node."""
        return self._parent_node
    def _set_parent_node(self, parent):
        """Sets the parent node of this node."""
        if self._parent_node is not None:
            try:
                self._parent_node._child_nodes.remove(self)
            except ValueError:
                pass
        self._parent_node = parent
        if self._parent_node is not None:
            if self not in self._parent_node._child_nodes:
                self._parent_node._child_nodes.append(self)
    parent_node = property(_get_parent_node, _set_parent_node)

    ###########################################################################
    ### General Structural Access and Information

    def is_leaf(self):
        """
        Returns |True| if the node is a tip or a leaf node, i.e. has no child
        nodes.

        Returns
        -------
        boolean
            |True| if the node is a leaf, i.e., has no child nodes. |False|
            otherwise.
        """
        return bool(not self._child_nodes)

    def is_internal(self):
        """
        Returns |True| if the node is *not* a tip or a leaf node.

        Returns
        -------
        boolean
            |True| if the node is not a leaf. |False| otherwise.
        """
        return bool(self._child_nodes)

    def leaf_nodes(self):
        """
        Returns list of all leaf_nodes descended from this node (or just
        list with ``self`` as the only member if ``self`` is a leaf).

        Note
        ----
        Usage of  `leaf_iter()` is preferable for efficiency reasons unless
        actual list is required.

        Returns
        -------
        :py:class:`list` [|Node|]
           A ``list`` of |Node| objects descended from this node
           (inclusive of ``self``) that are the leaves.
        """
        return [node for node in \
                self.postorder_iter(lambda x: bool(len(x.child_nodes())==0))]

    def num_child_nodes(self):
        """
        Returns number of child nodes.

        Returns
        -------
        int
            Number of children in ``self``.
        """
        return len(self._child_nodes)

    def child_nodes(self):
        """
        Returns a shallow-copy list of all child nodes of this node.

        Note
        ----
        Unless an actual ``list`` is needed, iterating over the child nodes using
        :meth:`Node.child_node_iter()` is preferable to avoid the overhead of
        list construction.

        Returns
        -------
        :py:class:`list` [|Node|]
           A ``list`` of |Node| objects that have ``self`` as a parent.
        """
        return list(self._child_nodes)

    def child_edges(self):
        """
        Returns a shallow-copy list of all child edges of this node.

        Note
        ----
        Unless an actual ``list`` is needed, iterating over the child edges using
        :meth:`Node.child_edge_iter()` is preferable to avoid the overhead of
        list construction.

        Returns
        -------
        :py:class:`list` [|Edge|]
           A ``list`` of |Edge| objects that have ``self`` as a tail node.
        """
        return list(ch.edge for ch in self._child_nodes)

    def incident_edges(self):
        """
        Return parent and child edges.

        Returns
        -------
        :py:class:`list` [|Edge|]
            A list of edges linking to this node, with outgoing edges (edges
            connecting to child nodes) followed by the edge connecting
            this node to its parent.
        """
        e = [c.edge for c in self._child_nodes]
        e.append(self.edge)
        return e

    def get_incident_edges(self):
        """Legacy synonym for :meth:`Node.incident_edges()`."""
        return self.incident_edges()

    def adjacent_nodes(self):
        """
        Return parent and child nodes.

        Returns
        -------
        :py:class:`list` [|Node|]
            A list with all child nodes and parent node of this node.
        """
        n = [c for c in self._child_nodes]
        if self._parent_node:
            n.append(self._parent_node)
        return n

    def get_adjacent_nodes(self):
        """Legacy synonym for :meth:`Node.adjacent_edges()`"""
        return self.adjacent_nodes()

    def sibling_nodes(self):
        """
        Return all other children of parent, excluding self.

        Returns
        -------
        :py:class:`list` [|Node|]
            A list of all nodes descended from the same parent as ``self``,
            excluding ``self``.
        """
        p = self._parent_node
        if not p:
            return []
        sisters = [nd for nd in p.child_nodes() if nd is not self]
        return sisters

    def sister_nodes(self):
        """Legacy synonym for :meth:`Node.sister_nodes()`"""
        return self.sibling_nodes()

    def extract_subtree(self,
            extraction_source_reference_attr_name="extraction_source",
            node_filter_fn=None,
            suppress_unifurcations=True,
            is_apply_filter_to_leaf_nodes=True,
            is_apply_filter_to_internal_nodes=False,
            node_factory=None,
            ):
        """
        Returns a clone of the structure descending from this node.

        Parameters
        ----------
        extraction_source_reference_attr_name : str
            Name of attribute to set on cloned nodes that references
            corresponding original node. If ``None``, then attribute (and
            reference) will not be created.
        node_filter_fn : None or function object
            If ``None``, then entire tree structure is cloned.
            If not ``None``, must be a function object that returns ``True``
            if a particular |Node| instance on the original tree should
            be included in the cloned tree, or ``False`` otherwise.
        is_apply_filter_to_leaf_nodes : bool
            If ``True`` then the above filter will be applied to leaf nodes. If
            ``False`` then it will not (and all leaf nodes will be
            automatically included, unless excluded by an ancestral node being
            filtered out).
        is_apply_filter_to_internal_nodes : bool
            If ``True`` then the above filter will be applied to internal nodes. If
            ``False`` then it will not (internal nodes without children will
            still be filtered out).
        node_factory : function
            If not ``None``, must be a function that takes no arguments and
            returns a new |Node| (or equivalent) instance.

        Returns
        -------
        nd : |Node|
            A node with descending subtree mirroring this one.

        """
        memo = {}
        is_excluded_nodes = False
        start_node = None
        start_node_to_match = self
        if node_factory is None:
            node_factory = self.__class__
        for nd0 in self.postorder_iter():
            if node_filter_fn is not None:
                if nd0._child_nodes:
                    if is_apply_filter_to_internal_nodes:
                        is_apply_filter = True
                    else:
                        is_apply_filter = False
                else:
                    if is_apply_filter_to_leaf_nodes:
                        is_apply_filter = True
                    else:
                        is_apply_filter = False
                if is_apply_filter and not node_filter_fn(nd0):
                    is_excluded_nodes = True
                    continue
            original_node_has_children = False
            children_to_add = []
            for ch_nd0 in nd0.child_node_iter():
                original_node_has_children = True
                ch_nd1 = memo.get(ch_nd0, None)
                if ch_nd1 is not None:
                    children_to_add.append(ch_nd1)
            if not children_to_add and original_node_has_children:
                # filter removes all descendents of internal node,
                # so this internal node is not added
                if nd0.parent_node is None:
                    raise error.SeedNodeDeletionException("Attempting to remove seed node or node without parent")
                if nd0 is self:
                    start_node_to_match = nd0.parent_node
                continue
            elif len(children_to_add) == 1 and suppress_unifurcations:
                if nd0.edge.length is not None:
                    if children_to_add[0].edge.length is None:
                        children_to_add[0].edge.length = nd0.edge.length
                    else:
                        children_to_add[0].edge.length += nd0.edge.length
                else:
                    nd1.edge.length = children_to_add[0].edge.length
                if nd0.parent_node is None:
                    start_node = children_to_add[0]
                    break
                if nd0 is self:
                    start_node_to_match = nd0.parent_node
                memo[nd0] = children_to_add[0]
            else:
                nd1 = node_factory()
                nd1.label = nd0.label
                nd1.taxon = nd0.taxon
                nd1.edge.length = nd0.edge.length
                nd1.edge.label = nd0.edge.label
                for ch_nd1 in children_to_add:
                    nd1.add_child(ch_nd1)
                if nd0 is start_node_to_match:
                    start_node = nd1
                memo[nd0] = nd1
                if extraction_source_reference_attr_name:
                    setattr(nd1, extraction_source_reference_attr_name, nd0)
        if start_node is not None:
            return start_node
        else:
            ## TODO: find a replacement node
            raise ValueError

    ###########################################################################
    ### Metrics

    def level(self):
        """
        Returns the number of nodes between ``self`` and the seed node of the tree.

        Returns
        -------
        integer
            The number of nodes between ``self`` and the seed node of the tree,
            or 0 if ``self`` has no parent.
        """
        if self._parent_node:
            return self._parent_node.level() + 1
        else:
            return 0

    def distance_from_root(self):
        """
        Weighted path length of ``self`` from root.

        Returns
        -------
        numeric
            Total weight of all edges connecting ``self`` with the root of the
            tree.
        """
        if self._parent_node and self.edge.length != None:
            if self._parent_node.distance_from_root == None:
                return float(self.edge.length)
            else:
                distance_from_root = float(self.edge.length)
                parent_node = self._parent_node
                # The root is identified when a node with no
                # parent is encountered. If we want to use some
                # other criteria (e.g., where a is_root property
                # is True), we modify it here.
                while parent_node:
                    if parent_node.edge.length != None:
                        distance_from_root = distance_from_root + float(parent_node.edge.length)
                    parent_node = parent_node._parent_node
                return distance_from_root
        elif not self._parent_node and self.edge.length != None:
            return float(self.edge.length)
        elif self._parent_node and self.edge.length == None:
            # what do we do here: parent node exists, but my
            # length does not?
            return float(self._parent_node.edge.length)
        elif not self._parent_node and self.edge.length == None:
            # no parent node, and no edge length
            return 0.0
        else:
            # WTF????
            return 0.0

    def distance_from_tip(self):
        """
        Maximum weighted length of path of ``self`` to tip.

        If tree is not ultrametric (i.e., descendent edges have different
        lengths), then count the maximum of edge lengths. Note that
        :meth:`Tree.calc_node_ages()` is a more efficient way of doing this
        over the whole tree if this value is need for many or all the nodes on
        the tree.

        Returns
        -------
        numeric
            Maximum weight of edges connecting ``self`` to tip.
        """
        if not self._child_nodes:
            return 0.0
        else:
            distance_from_tips = []
            for ch in self._child_nodes:
                if ch.edge.length is not None:
                    curr_edge_length = ch.edge_length
                else:
                    curr_edge_length = 0.0
                if not hasattr(ch, "_distance_from_tip"):
                    ch._distance_from_tip = ch.distance_from_tip()
                distance_from_tips.append(ch._distance_from_tip + curr_edge_length)
            self._distance_from_tip = float(max(distance_from_tips))
            return self._distance_from_tip

    ###########################################################################
    ### Representation

    def description(self, depth=1, indent=0, itemize="", output=None, taxon_namespace=None):
        """
        Returns description of object, up to level ``depth``.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        label = str(self)
        output_strio.write('%s%sNode object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            leader1 = ' ' * (indent + 4)
            leader2 = ' ' * (indent + 8)
            output_strio.write('\n%s[Edge]' % leader1)
            if self.edge is not None:
                edge_desc = self.edge.description(0)
            else:
                edge_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, edge_desc))

            output_strio.write('\n%s[Taxon]' % leader1)
            if self.taxon is not None:
                taxon_desc = self.taxon.description(0)
            else:
                taxon_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, taxon_desc))

            output_strio.write('\n%s[Parent]' % leader1)
            if self._parent_node is not None:
                parent_node_desc = self._parent_node.description(0)
            else:
                parent_node_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, parent_node_desc))
            output_strio.write('\n%s[Children]' % leader1)
            if len(self._child_nodes) == 0:
                output_strio.write('\n%sNone' % leader2)
            else:
                for i, cnd in enumerate(self._child_nodes):
                    output_strio.write('\n%s[%d] %s' % (leader2, i, cnd.description(0)))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    ###########################################################################
    ### Native NEWICK printer
    ## For debugging we build-in a full-fledged NEWICK composition independent
    ## of the nexus/newick family of modules. Client code should prefer to
    ## use Newick/Nexus readers/writers, or Tree.write(), TreeList.write(),
    ## DataSet.write() etc.

    def _as_newick_string(self, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.
        For production purposes, use the the full-fledged 'as_string()'
        method of the object.
        """
        out = StringIO()
        self._write_newick(out, **kwargs)
        return out.getvalue()

    def _write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.  For
        production purposes, use the the full-fledged 'write_to_stream()'
        method of the object.
        """
        edge_lengths = not kwargs.get('suppress_edge_lengths', False)
        edge_lengths = kwargs.get('edge_lengths', edge_lengths)
        child_nodes = self.child_nodes()
        if child_nodes:
            out.write('(')
            f_child = child_nodes[0]
            for child in child_nodes:
                if child is not f_child:
                    out.write(',')
                child._write_newick(out, **kwargs)
            out.write(')')

        out.write(self._get_node_token(**kwargs))
        if edge_lengths:
            e = self.edge
            if e:
                sel = e.length
                if sel is not None:
                    fmt = kwargs.get('edge_length_formatter', None)
                    if fmt:
                        out.write(":%s" % fmt(sel))
                    else:
                        s = ""
                        try:
                            s = float(sel)
                            s = str(s)
                        except ValueError:
                            s = str(sel)
                        if s:
                            out.write(":%s" % s)

    def _get_node_token(self, **kwargs):
        """returns a string that is an identifier for the node.  This is called
        by the newick-writing functions, so the kwargs that affect how node
        labels show up in a newick string are the same ones used here:
        ``suppress_internal_labels`` is a Boolean, and defaults to False.
        """
        is_leaf = (len(self._child_nodes) == 0)
        if not is_leaf:
            if kwargs.get("suppress_internal_labels", False) \
                    or not kwargs.get("include_internal_labels", True):
                return ""
        if self.taxon is not None:
            if self.taxon.label:
                label = self.taxon.label
            else:
                # return "_" # taxon, but no label: anonymous
                label = "" # "_" is not anonoymous/unnamed, but a name of <blank>; so we return nothing instead
        else:
            if self.label:
                label = self.label
            else:
                label = ""
        if not label or kwargs.get("raw_labels", False):
            return label
        elif " " in label and "_" in label:
            if "'" in label:
                label.replace("'", "''")
            return "'{}'".format(label)
        elif " " in label and not kwargs.get("preserve_spaces"):
            return label.replace(" ", "_")
        else:
            return label

    ###########################################################################
    ### alternate representation of tree structure for debugging

    def _get_indented_form(self, **kwargs):
        out = StringIO()
        self._write_indented_form(out, **kwargs)
        return out.getvalue()

    def _write_indented_form(self, out, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        level = kwargs.get("level", 0)
        ancestors = []
        siblings = []
        n = self
        while n is not None:
            n._write_indented_form_line(out, level, **kwargs)
            n, lev = _preorder_list_manip(n, siblings, ancestors)
            level += lev

    def _get_indented_form_line(self, level, **kwargs):
        out = StringIO()
        self._write_indented_form_line(out, level, **kwargs)
        return out.getvalue()

    def _write_indented_form_line(self, out, level, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        label = _format_node(self, **kwargs)
        if kwargs.get("bipartitions"):
            cm = "%s " % _format_bipartition(self.edge.bipartition, **kwargs)
        else:
            cm = ""
        out.write("%s%s%s\n" % ( cm, indentation*level, label))

##############################################################################
### Tree

class Tree(
        taxonmodel.TaxonNamespaceAssociated,
        basemodel.Annotable,
        basemodel.Deserializable,
        basemodel.NonMultiReadable,
        basemodel.Serializable,
        basemodel.DataObject):
    """
    An arborescence, i.e. a fully-connected directed acyclic graph with all
    edges directing away from the root and toward the tips. The "root" of the
    tree is represented by the :attr:`Tree.seed_node` attribute.  In unrooted
    trees, this node is an algorithmic artifact. In rooted trees this node is
    semantically equivalent to the root.
    """

    def _parse_and_create_from_stream(cls,
            stream,
            schema,
            collection_offset=None,
            tree_offset=None,
            **kwargs):
        """
        Constructs a new |Tree| object and populates it with data from
        file-like object ``stream``.

        If the source defines multiple tree collections (e.g. multiple NEXUS
        "Trees" blocks), then the ``collection_offset`` argument can be
        used to specify the 0-based index of the tree collection, and
        ``tree_offset`` argument can be used to specify the 0-based index of
        the tree within the collection, as the source. If ``collection_offset``
        is not specified or |None|, then the first collection (offset=0) is
        assumed. If ``tree_offset`` is not specified, then the first tree
        (offset=0) is returned.

        Keyword arguments `**kwargs` are passed directly to
        :meth:|TreeList|.read()`, which wraps the actual parsing.

        If no tree is found in the source according to the specified criteria,
        then |None| is returned.

        Notes
        -----
        *All* operational taxonomic unit concepts in the data source will be included
        in the |TaxonNamespace| object associated with the new
        |TreeList| object and its contained |Tree| objects, even those
        not associated with tree being retrieved.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in ``stream``

        collection_offset : integer
            0-based index of tree block or collection in source to be parsed.

        tree_offset : integer
            0-based index of tree in source to be parsed.  This is the 0-based
            index of the tree within the collection specified by
            ``collection_offset`` to be retrieved.

        \*\*kwargs : keyword arguments
            Arguments to customize parsing and instantiation this |Tree|
            from the data source, including schema- or format-specific
            handling. The following optional keyword arguments are recognized
            and handled by this constructor:

                ``label``
                    The label or description of the new |Tree| object.
                ``taxon_namespace``
                   Specifies the |TaxonNamespace| object to be attached
                   to the new |TreeList| object. Note that *all*
                   operational taxonomic unit concepts in the data source will
                   be accessioned into the specified |TaxonNamespace|
                   instance. This includes the operation taxonomic unit
                   definitions associated with all tree collections and
                   character matrices in the data source.

            Other keyword arguments may be available, depending on the
            implementation of the reader specialized to handle ``schema``
            formats. See documentation for details on keyword arguments
            supported by readers of various schemas.

        Returns
        -------
        |Tree| or |None|
            The |Tree| object corresponding to the tree in the data
            source, or |None| if no valid tree description was found.

        """
        from dendropy.datamodel.treecollectionmodel import TreeList

        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
        if taxon_namespace is None:
            taxon_namespace = taxonmodel.TaxonNamespace()

        def tns_factory(label):
            if label is not None and taxon_namespace.label is None:
                taxon_namespace.label = label
            return taxon_namespace

        tree_list_factory = lambda label, taxon_namespace: TreeList(label=label, taxon_namespace=taxon_namespace, tree_type=cls)
        label = kwargs.pop("label", None)
        reader = dataio.get_reader(schema, **kwargs)
        # if collection_offset is None and tree_offset is not None:
        #     raise TypeError("Cannot specify ``tree_offset`` without specifying ``collection_offset``")
        if collection_offset is None:
            collection_offset = 0
        if tree_offset is None:
            tree_offset = 0
        tree_lists = reader.read_tree_lists(
                    stream=stream,
                    taxon_namespace_factory=tns_factory,
                    tree_list_factory=tree_list_factory,
                    global_annotations_target=None)
        if not tree_lists:
            raise ValueError("No trees in data source")
        tree_list = tree_lists[collection_offset]
        if not tree_list:
            raise ValueError("No trees available at requested location in data source")
        tree = tree_list[tree_offset]
        tree.label = label
        return tree
    _parse_and_create_from_stream = classmethod(_parse_and_create_from_stream)

    @classmethod
    def get(cls, **kwargs):
        """
        Instantiate and return a *new* |Tree| object from a data source.

        **Mandatory Source-Specification Keyword Argument (Exactly One of the Following Required):**

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

            - **label** (*str*) -- Name or identifier to be assigned to the new
              object; if not given, will be assigned the one specified in the
              data source, or |None| otherwise.
            - **taxon_namespace** (|TaxonNamespace|) -- The |TaxonNamespace|
              instance to use to :doc:`manage the taxon names </primer/taxa>`.
              If not specified, a new one will be created.
            - **collection_offset** (*int*) -- 0-based index of tree block or
              collection in source to be parsed. If not specified then the
              first collection (offset = 0) is assumed.
            - **tree_offset** (*int*) -- 0-based index of tree within the
              collection specified by ``collection_offset`` to be parsed. If
              not specified, then the first tree (offset = 0) is assumed.
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

            # From a URL
            t1 = dendropy.Tree.get(
                    url="http://api.opentreeoflife.org/v2/study/pg_1144/tree/tree2324.nex",
                    schema="nexus")

            # From a file-like object
            t2 = Tree.get(file=open('treefile.tre', 'r'),
                            schema="newick",
                            tree_offset=0)

            # From a path
            t3 = Tree.get(path='sometrees.nexus',
                    schema="nexus",
                    collection_offset=2,
                    tree_offset=1)

            # From a string
            s = "((A,B),(C,D));((A,C),(B,D));"
            # tree will be '((A,B),(C,D))'
            t4 = Tree.get(data=s,
                    schema="newick")
            # tree will be '((A,C),(B,D))'
            t5 = Tree.get(data=s,
                    schema="newick",
                    tree_offset=1)
            # passing keywords to underlying tree parser
            t7 = dendropy.Tree.get(
                    data="((A,B),(C,D));",
                    schema="newick",
                    taxon_namespace=t3.taxon_namespace,
                    suppress_internal_node_taxa=False,
                    preserve_underscores=True)

        """
        return cls._get_from(**kwargs)

    def yield_from_files(cls,
            files,
            schema,
            taxon_namespace=None,
            **kwargs):
        """
        Iterates over trees from files, returning them one-by-one instead of
        instantiating all of them in memory at once.

        For operations where it is sufficient to process each tree individually
        (e.g., performing a calculation or set of calculations on a tree and
        storing the result, after the which the entire tree itself is
        not needed), this approach is *far* more performant that reading in the
        trees using a |TreeList|. This is because a full tree structure
        requires significant memory overhead, and as memory gets used up and
        the OS starts page faulting, performance starts taking some serious
        hits.

        Parameters
        ----------
        files : iterable of file paths or file-like objects.
            Iterable of sources, which can either be strings specifying file
            paths or file-like objects open for reading. If a source element is
            a string (``isinstance(i,str) == True``), then it is assumed to be
            a path to a file. Otherwise, the source is assumed to be a file-like
            object.
        schema : string
            The name of the data format (e.g., "newick" or "nexus").
        taxon_namespace : |TaxonNamespace| instance
            The operational taxonomic unit concept namespace to use to manage
            taxon definitions.
        \*\*kwargs : keyword arguments
            These will be passed directly to the schema-parser implementation.

        Yields
        ------
        t : |Tree|
            Trees as read from the file.

        Examples
        --------

        ::

            taxon_namespace = dendropy.TaxonNamespace()
            f1 = open("path/to/trees1.nex", "r")
            f2 = open("path/to/trees2.nex", "r")
            tree_yielder = dendropy.Tree.yield_from_files(
                    files=[f1, f2, "path/to/trees3.nex", "path/to/trees4.nex"],
                    schema="nexus",
                    taxon_namespace=taxon_namespace,
                    store_tree_weights=True,
                    preserve_underscores=True,
                    rooting="default-unrooted",
                    ignore_unrecognized_keyword_arguments=True,
                    )
            lengths = []
            root_ages = []
            for tree in tree_yielder:
                length = 0.0
                for edge in tree:
                    length += edge.length
                lengths.append(length)
                tree.calc_node_ages()
                root_ages.append(tree.seed_node.age)

        """
        if taxon_namespace is None:
            taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None)
            if taxon_namespace is None:
                taxon_namespace = taxonmodel.TaxonNamespace()
        else:
            assert "taxon_set" not in kwargs
        if "tree_offset" in kwargs:
            raise TypeError("'tree_offset' is not supported: trees should be skipped/discarded on the client code side")
        tree_yielder = dataio.get_tree_yielder(
                files,
                schema,
                taxon_namespace=taxon_namespace,
                tree_type=cls,
                **kwargs)
        return tree_yielder
    yield_from_files = classmethod(yield_from_files)

    def from_bipartition_encoding(
            cls,
            bipartition_encoding,
            taxon_namespace,
            is_rooted=False,
            edge_lengths=None,
            ):
        """
        Reconstructs a tree from a bipartition encoding.

        Parameters
        ----------
        bipartition_encoding : iterable[|Bipartition|]
            An iterable of |Bipartition| instances representing a tree.
            Bipartitions will be added to the tree in the order given by
            iterating over this. Bipartitions that are incompatible with the
            tree will be skipped. So if not all bipartitions are compatible
            with each other, then the sequence of bipartitions given in
            ``bipartition_encoding`` should be in order of their support values
            or some other preference criteria.
        taxon_namespace : |TaxonNamespace| instance
            The operational taxonomic unit concept namespace to use to manage
            taxon definitions.
        is_rooted : bool
            Specifies whether or not the tree is rooted.
        edge_lengths : iterable or |None|
            An iterable of edge lengths. This should be in the same order
            as the bipartitions in the bipartition encoding.

        Returns
        -------
        |Tree|
            The tree reconstructed from the given bipartition encoding.
        """
        split_bitmasks = [b.split_bitmask for b in bipartition_encoding]
        if edge_lengths:
            split_edge_lengths = dict(zip(split_bitmasks, edge_lengths))
        else:
            split_edge_lengths = None
        # elif edge_lengths is not False:
        #     split_edge_lengths = dict(zip(split_bitmasks,
        #         [b.edge.length for b in bipartition_encoding]))
        return cls.from_split_bitmasks(
                split_bitmasks=split_bitmasks,
                taxon_namespace=taxon_namespace,
                split_edge_lengths=split_edge_lengths,
                is_rooted=is_rooted)
    from_bipartition_encoding = classmethod(from_bipartition_encoding)

    def from_split_bitmasks(
            cls,
            split_bitmasks,
            taxon_namespace,
            is_rooted=False,
            split_edge_lengths=None,
            ):
        """
        Reconstructs a tree from a collection of splits represented as bitmasks.

        Parameters
        ----------
        split_bitmasks : iterable[int]
            An iterable of split bitmasks representing a tree.
            Splits will be added to the tree in the order given by
            iterating over this. Splits that are incompatible with the
            tree will be skipped. So if not all splits are compatible
            with each other, then the sequence of splits given in
            ``bipartition_encoding`` should be in order of their support values
            or some other preference criteria.
        taxon_namespace : |TaxonNamespace| instance
            The operational taxonomic unit concept namespace to use to manage
            taxon definitions.
        is_rooted : bool
            Specifies whether or not the tree is rooted.
        edge_lengths : dict or |False| or |None|
            If |False| or |None|, then no edge lengths will be added.
            Otherwise, this should be a dictionary mapping splits to edge
            lengths.
        Returns
        -------
        |Tree|
            The tree reconstructed from the given bipartition encoding.
        """
        leaf_to_root_search = True
        reconstructed_tree = cls(taxon_namespace=taxon_namespace)
        # reconstructed_tree.is_rooted = True
        reconstructed_tree.is_rooted = is_rooted
        for taxon in taxon_namespace:
            reconstructed_tree.seed_node.new_child(taxon=taxon)
        all_taxa_bitmask = taxon_namespace.all_taxa_bitmask()
        reconstructed_tree.encode_bipartitions()
        reconstructed_tree.bipartition_encoding = []
        leaves = reconstructed_tree.leaf_nodes()
        if leaf_to_root_search:
            to_leaf_dict = {}
            for leaf in leaves:
                to_leaf_dict[leaf.edge.bipartition.leafset_bitmask] = leaf
        root = reconstructed_tree.seed_node
        root_edge = root.edge

        split_bitmasks_to_add = []
        for s in split_bitmasks:
            m = s & all_taxa_bitmask
            if (m != all_taxa_bitmask) and ((m-1) & m): # if not root (i.e., all "1's") and not singleton (i.e., one "1")
                if is_rooted:
                    split_bitmasks_to_add.append(m)
                else:
                    if 1 & m:
                        split_bitmasks_to_add.append( (~m) & all_taxa_bitmask )
                    else:
                        # "denormalize" split_bitmasks
                        split_bitmasks_to_add.append(m)

        # Now when we add split_bitmasks in order, we will do a greedy, extended majority-rule consensus tree
        #for freq, split_to_add, split_in_dict in to_try_to_add:
        _get_mask = lambda x: getattr(x.bipartition, "_leafset_bitmask")
        for split_to_add in split_bitmasks_to_add:
            if (split_to_add & root_edge.bipartition.leafset_bitmask) != split_to_add:
                # incompatible
                continue
            elif leaf_to_root_search:
                lb = bitprocessing.least_significant_set_bit(split_to_add)
                one_leaf = to_leaf_dict[lb]
                parent_node = one_leaf
                while (split_to_add & parent_node.edge.bipartition.leafset_bitmask) != split_to_add:
                    parent_node = parent_node.parent_node
            else:
                parent_node = reconstructed_tree.mrca(split_bitmask=split_to_add)
            if parent_node is None or parent_node.edge.bipartition.leafset_bitmask == split_to_add:
                continue # split is not in tree, or already in tree.
            new_node = cls.node_factory()
            #self.map_split_support_to_node(node=new_node, split_support=freq)
            new_node_children = []
            new_edge = new_node.edge
            new_mask = 0
            for child in parent_node.child_nodes():
                # might need to modify the following if rooted split_bitmasks
                # are used
                cecm = child.edge.bipartition.leafset_bitmask
                if (cecm & split_to_add):
                    assert cecm != split_to_add
                    new_mask |= cecm
                    new_node_children.append(child)
                    new_edge.bipartition = Bipartition(
                            leafset_bitmask=new_mask,
                            tree_leafset_bitmask=all_taxa_bitmask,
                            is_mutable=False,
                            compile_bipartition=True)
                    reconstructed_tree.bipartition_encoding.append(new_edge.bipartition)
            # Check to see if we have accumulated all of the bits that we
            #   needed, but none that we don't need.
            if new_edge.bipartition.leafset_bitmask == split_to_add:
                if split_edge_lengths:
                    new_edge.length = split_edge_lengths[split_to_add]
                    #old_split = new_old_split_map[split_to_add]
                    #new_edge.length = split_edge_lengths[old_split]
                for child in new_node_children:
                    parent_node.remove_child(child)
                    new_node.add_child(child)
                parent_node.add_child(new_node)
                # reconstructed_tree.split_edge_map[split_to_add] = new_edge
        return reconstructed_tree
    from_split_bitmasks = classmethod(from_split_bitmasks)

    def node_factory(cls, **kwargs):
        """
        Creates and returns a |Node| object.

        Derived classes can override this method to provide support for
        specialized or different types of nodes on the tree.

        Parameters
        ----------

        \*\*kwargs : keyword arguments
            Passed directly to constructor of |Node|.

        Returns
        -------
        |Node|
            A new |Node| object.

        """
        return Node(**kwargs)
    node_factory = classmethod(node_factory)

    ###########################################################################
    ### Special/Lifecycle methods

    def __init__(self, *args, **kwargs):
        """
        The constructor can optionally construct a |Tree| object by
        cloning another |Tree| object passed as the first positional
        argument, or out of a data source if ``stream`` and ``schema`` keyword
        arguments are passed with a file-like object and a schema-specification
        string object values respectively.

        Parameters
        ----------

        \*args : positional argument, optional
            If given, should be exactly one |Tree| object. The new
            |Tree| will then be a structural clone of this argument.

        \*\*kwargs : keyword arguments, optional
            The following optional keyword arguments are recognized
            and handled by this constructor:

                ``label``
                    The label or description of the new |Tree| object.
                ``taxon_namespace``
                    Specifies the |TaxonNamespace| object to be
                    that the new |Tree| object will reference.
        Examples
        --------

        Tree objects can be instantiated in the following ways::

            # /usr/bin/env python

            try:
                from StringIO import StringIO
            except ImportError:
                from io import StringIO
            from dendropy import Tree, TaxonNamespace

            # empty tree
            t1 = Tree()

            # Tree objects can be instantiated from an external data source
            # using the 'get()' factory class method

            # From a file-like object
            t2 = Tree.get(file=open('treefile.tre', 'r'),
                            schema="newick",
                            tree_offset=0)

            # From a path
            t3 = Tree.get(path='sometrees.nexus',
                    schema="nexus",
                    collection_offset=2,
                    tree_offset=1)

            # From a string
            s = "((A,B),(C,D));((A,C),(B,D));"
            # tree will be '((A,B),(C,D))'
            t4 = Tree.get(data=s,
                    schema="newick")
            # tree will be '((A,C),(B,D))'
            t5 = Tree.get(data=s,
                    schema="newick",
                    tree_offset=1)
            # passing keywords to underlying tree parser
            t7 = dendropy.Tree.get(
                    data="((A,B),(C,D));",
                    schema="newick",
                    taxon_namespace=t3.taxon_namespace,
                    suppress_internal_node_taxa=False,
                    preserve_underscores=True)

            # Tree objects can be written out using the 'write()' method.
            t1.write(file=open('treefile.tre', 'r'),
                    schema="newick")
            t1.write(path='treefile.nex',
                    schema="nexus")

            # Or returned as a string using the 'as_string()' method.
            s = t1.as_string("nexml")

            # tree structure deep-copied from another tree
            t8 = dendropy.Tree(t7)
            assert t8 is not t7                             # Trees are distinct
            assert t8.symmetric_difference(t7) == 0         # and structure is identical
            assert t8.taxon_namespace is t7.taxon_namespace             # BUT taxa are not cloned.
            nds3 = [nd for nd in t7.postorder_node_iter()]  # Nodes in the two trees
            nds4 = [nd for nd in t8.postorder_node_iter()]  # are distinct objects,
            for i, n in enumerate(nds3):                    # and can be manipulated
                assert nds3[i] is not nds4[i]               # independentally.
            egs3 = [eg for eg in t7.postorder_edge_iter()]  # Edges in the two trees
            egs4 = [eg for eg in t8.postorder_edge_iter()]  # are also distinct objects,
            for i, e in enumerate(egs3):                    # and can also be manipulated
                assert egs3[i] is not egs4[i]               # independentally.
            lves7 = t7.leaf_nodes()                         # Leaf nodes in the two trees
            lves8 = t8.leaf_nodes()                         # are also distinct objects,
            for i, lf in enumerate(lves3):                  # but order is the same,
                assert lves7[i] is not lves8[i]             # and associated Taxon objects
                assert lves7[i].taxon is lves8[i].taxon     # are the same.

            # To create deep copy of a tree with a different taxon namespace,
            # Use 'copy.deepcopy()'
            t9 = copy.deepcopy(t7)

            # Or explicitly pass in a new TaxonNamespace instance
            taxa = TaxonNamespace()
            t9 = dendropy.Tree(t7, taxon_namespace=taxa)
            assert t9 is not t7                             # As above, the trees are distinct
            assert t9.symmetric_difference(t7) == 0         # and the structures are identical,
            assert t9.taxon_namespace is not t7.taxon_namespace         # but this time, the taxa *are* different
            assert t9.taxon_namespace is taxa                     # as the given TaxonNamespace is used instead.
            lves3 = t7.leaf_nodes()                         # Leaf nodes (and, for that matter other nodes
            lves5 = t9.leaf_nodes()                         # as well as edges) are also distinct objects
            for i, lf in enumerate(lves3):                  # and the order is the same, as above,
                assert lves7[i] is not lves9[i]             # but this time the associated Taxon
                assert lves7[i].taxon is not lves9[i].taxon # objects are distinct though the taxon
                assert lves7[i].taxon.label == lves9[i].taxon.label # labels are the same.

            # to 'switch out' the TaxonNamespace of a tree, replace the reference and
            # reindex the taxa:
            t11 = Tree.get(data='((A,B),(C,D));', 'newick')
            taxa = TaxonNamespace()
            t11.taxon_namespace = taxa
            t11.reindex_subcomponent_taxa()

            # You can also explicitly pass in a seed node:
            seed = Node(label="root")
            t12 = Tree(seed_node=seed)
            assert t12.seed_node is seed

        """
        if len(args) > 1:
            # only allow 1 positional argument
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        elif len(args) == 1:
            if "seed_node" in kwargs:
                raise TypeError("Cannot specify 'seed_node' if passing in a Tree object to clone")
            if "stream" in kwargs or "schema" in kwargs:
                raise TypeError("Constructing from an external stream is no longer supported: use the factory method 'Tree.get(file=...)'")
            if isinstance(args[0], Node):
                raise TypeError("Constructing a tree around a Node passed as a position argument is no longer supported; a keyword argument is now required for this approach: use Tree(seed_node=node)")
            if isinstance(args[0], Tree):
                self._clone_from(args[0], kwargs)
            else:
                raise error.InvalidArgumentValueError(func_name=self.__class__.__name__, arg=args[0])
        else:
            basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
            taxonmodel.TaxonNamespaceAssociated.__init__(self,
                    taxon_namespace=taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs, None))
            self.comments = []
            self._is_rooted = None
            self.weight = None
            self.length_type = None
            self._seed_node = None
            self.seed_node = None
            self.bipartition_encoding = None
            self._split_bitmask_edge_map = None
            self._bipartition_edge_map = None
            seed_node = kwargs.pop("seed_node", None)
            if seed_node is None:
                self.seed_node = self.node_factory()
            else:
                self.seed_node = seed_node
                self.update_taxon_namespace()
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    ##############################################################################
    ## Bipartitions

    def _get_split_edges(self):
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'Tree.split_edges' will no longer be supported in future releases; use 'Tree.bipartition_encoding' for a list of bipartitions on the tree, or dereference the edge through the 'Tree.bipartition_edge_map' attribute.",
                stacklevel=3)
        return self.bipartition_encoding
    def _set_split_edges(self, m):
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'Tree.split_edges' will no longer be supported in future releases; use 'Tree.bipartition_encoding' for a list of bipartitions on the tree, or dereference the edge through the 'Tree.bipartition_edge_map' attribute.",
                stacklevel=3)
        self.bipartition_encoding = m
    split_edges = property(_get_split_edges, _set_split_edges)

    ##############################################################################
    ## Identity

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    ##############################################################################
    ## Copying/cloning

    def _clone_from(self, tree, kwargs_dict):
        # super(Tree, self).__init__()
        memo = {}
        # memo[id(tree)] = self
        taxon_namespace = taxonmodel.process_kwargs_dict_for_taxon_namespace(kwargs_dict, tree.taxon_namespace)
        memo[id(tree.taxon_namespace)] = taxon_namespace
        if taxon_namespace is not tree.taxon_namespace:
            for t1 in tree.taxon_namespace:
                t2 = taxon_namespace.require_taxon(label=t1.label)
                memo[id(t1)] = t2
        else:
            for t1 in tree.taxon_namespace:
                memo[id(t1)] = t1
        t = copy.deepcopy(tree, memo)
        self.__dict__ = t.__dict__
        self.label = kwargs_dict.pop("label", tree.label)
        return self
        # for k in tree.__dict__:
        #     if k == "_annotations":
        #         continue
        #     if k in self.__dict__:
        #         # do not copy if already populated, perhaps by a derived class
        #         continue
        #     self.__dict__[k] = copy.deepcopy(tree.__dict__[k], memo)
        #     memo[id(tree.__dict__[k])] = self.__dict__[k]
        # self.deep_copy_annotations_from(tree)

    def __copy__(self):
        return self.taxon_namespace_scoped_copy()

    def taxon_namespace_scoped_copy(self, memo=None):
        if memo is None:
            memo = {}
        # this populates ``memo`` with references to the
        # the TaxonNamespace and Taxon objects
        self.taxon_namespace.populate_memo_for_taxon_namespace_scoped_copy(memo)
        return self.__deepcopy__(memo=memo)

    def __deepcopy__(self, memo=None):
        # ensure clone map
        return basemodel.Annotable.__deepcopy__(self, memo=memo)
        # if memo is None:
        #     memo = {}
        # # get or create clone of self
        # try:
        #     other = memo[id(self)]
        # except KeyError:
        #     # create object without initialization
        #     other = self.__class__.__new__(self.__class__)
        #     # store
        #     memo[id(self)] = other
        # # copy other attributes first, skipping annotations
        # for k in self.__dict__:
        #     if k == "_annotations":
        #         continue
        #     if k in other.__dict__:
        #         # do not copy if already populated, perhaps by a derived class
        #         continue
        #     other.__dict__[k] = copy.deepcopy(self.__dict__[k], memo)
        #     memo[id(self.__dict__[k])] = other.__dict__[k]
        #     # assert id(self.__dict__[k]) in memo
        # # create annotations
        # other.deep_copy_annotations_from(self, memo)
        # # return
        # return other

    ###########################################################################
    ### Extracting Trees and Subtrees

    def extract_tree(self,
            extraction_source_reference_attr_name="extraction_source",
            node_filter_fn=None,
            suppress_unifurcations=True,
            is_apply_filter_to_leaf_nodes=True,
            is_apply_filter_to_internal_nodes=False,
            tree_factory=None,
            node_factory=None,
            ):
        """
        Returns a copy of this tree that only includes the basic structure
        (nodes, edges), and minimal attributes (edge lengths, node labels, and
        taxon associations). Annotations, comments, and other attributes are
        not copied.

        Parameters
        ----------
        extraction_source_reference_attr_name : str
            Name of attribute to set on cloned nodes that references
            corresponding original node. If ``None``, then attribute (and
            reference) will not be created.
        node_filter_fn : None or function object
            If ``None``, then entire tree structure is cloned.
            If not ``None``, must be a function object that returns ``True``
            if a particular |Node| instance on the original tree should
            be included in the cloned tree, or ``False`` otherwise.
        suppress_unifurcations : bool
            If |True|, nodes of outdegree 1 will be deleted. Only will
            be done if some nodes are excluded from the cloned tree.
        is_apply_filter_to_leaf_nodes : bool
            If ``True`` then the above filter will be applied to leaf nodes. If
            ``False`` then it will not (and all leaf nodes will be
            automatically included, unless excluded by an ancestral node being
            filtered out).
        is_apply_filter_to_internal_nodes : bool
            If ``True`` then the above filter will be applied to internal nodes. If
            ``False`` then it will not (internal nodes without children will
            still be filtered out).
        tree_factory : function
            If not ``None``, must be a function that optionally takes a
            |TaxonNamespace| as an argument and returns a new |Tree| (or
            equivalent) instance.
        node_factory : function
            If not ``None``, must be a function that takes no arguments and
            returns a new |Node| (or equivalent) instance.

        Examples
        --------

        A simple clone::

            tree0 = dendropy.Tree.get(
                        path="mammals.tre",
                        schema="newick")
            tree1 = tree0.extract_tree()

        A clone that only extracts a subtree with taxa in the genus
        "Rhacophorus"::

            tree0 = dendropy.Tree.get(
                        path="old_world_frogs.tre",
                        schema="newick")
            # Include taxa only if label starts with "Rhacophorus"
            node_filter_fn = lambda nd: nd.is_internal() or \
                        nd.taxon.label.startswith("Rhacophorus")
            tree1 = tree0.extract_tree(node_filter_fn=node_filter_fn)

            # Above is equivalent to, but more efficient than:
            #   inclusion_set = [nd.taxon for nd in tree0.leaf_node_iter()
            #           if nd.taxon.label.startswith("Rhacophorus)]
            #   tree1 = dendropy.Tree(tree0)
            #   tree1.retain_taxa(inclusion_set)

        A clone that only extracts a subtree with nodes with taxa associated
        with the habitat "mountain" or "forest"::

            tree0 = dendropy.Tree.get(
                        path="birds.tre",
                        schema="newick")
            include_habitats = set(["mountain", "forest"])
            node_filter_fn = lambda nd: nd.taxon is None or \
                        nd.taxon.annotations["habitat"] in include_habitats
            tree1 = tree0.extract_tree(node_filter_fn=node_filter_fn)

        Returns
        -------
        t : |Tree|
            A new tree based on this one, with nodes filtered out if specified.

        """
        if tree_factory is None:
            other = self.__class__(taxon_namespace=self.taxon_namespace)
        else:
            other = tree_factory(taxon_namespace=self.taxon_namespace)
            if node_factory is None:
                try:
                    node_factory = other.node_factory
                except AttributeError:
                    pass
        other._is_rooted = self._is_rooted
        other.weight = self.weight
        other.length_type = self.length_type
        other.label = self.label
        other.seed_node = self.seed_node.extract_subtree(
                extraction_source_reference_attr_name=extraction_source_reference_attr_name,
                node_filter_fn=node_filter_fn,
                suppress_unifurcations=suppress_unifurcations,
                is_apply_filter_to_leaf_nodes=is_apply_filter_to_leaf_nodes,
                is_apply_filter_to_internal_nodes=is_apply_filter_to_internal_nodes,
                node_factory=node_factory,
                )
        return other

    def extract_tree_with_taxa(self,
            taxa,
            extraction_source_reference_attr_name="extraction_source",
            suppress_unifurcations=True,
            ):
        """
        Returns a copy of this tree that only includes leaf nodes if they
        are associated with the taxon objects listed in ``taxa``. Note that
        this copy will be a "thin" copy, including just the basic structure
        (nodes, edges) and minimal attributes (edge lengths, node labels, and
        taxon associations). Annotations, comments, and other attributes are
        not copied.

        Parameters
        ----------
        taxa : iterable of |Taxon| instances
            List or some other iterable of |Taxon| objects to include.
        suppress_unifurcations : bool
            If |True|, nodes of outdegree 1 will be deleted. Only will
            be done if some nodes are excluded from the cloned tree.
        is_apply_filter_to_leaf_nodes : bool
            If ``True`` then the above filter will be applied to leaf nodes. If
            ``False`` then it will not (and all leaf nodes will be
            automatically included, unless excluded by an ancestral node being
            filtered out).
        is_apply_filter_to_internal_nodes : bool
            If ``True`` then the above filter will be applied to internal nodes. If
            ``False`` then it will not (internal nodes without children will
            still be filtered out).

        Examples
        --------

        A clone that only extracts a subtree with taxa in the genus
        "Rhacophorus"::

            tree0 = dendropy.Tree.get(
                        path="old_world_frogs.tre",
                        schema="newick")
            # Include taxa only if label starts with "Rhacophorus"
            taxa_to_retain = set([taxon for taxon in tree0.taxon_namespace
                    if taxon.label.startswith("Rhacophorus")])
            tree1 = tree0.extract_tree_with_taxa(taxa=taxa_to_retain)

            # Above is equivalent to, but more efficient than:
            #   inclusion_set = [nd.taxon for nd in tree0.leaf_node_iter()
            #           if nd.taxon.label.startswith("Rhacophorus)]
            #   tree1 = dendropy.Tree(tree0)
            #   tree1.retain_taxa(inclusion_set)

        Returns
        -------
        t : |Tree|
            A new tree based on this one, with nodes filtered out if specified.

        """
        node_filter_fn = lambda nd: nd.taxon is None or nd.taxon in set(taxa)
        return self.extract_tree(
                node_filter_fn=node_filter_fn,
                extraction_source_reference_attr_name=extraction_source_reference_attr_name,
                is_apply_filter_to_leaf_nodes=True,
                is_apply_filter_to_internal_nodes=False,
                )

    def extract_tree_with_taxa_labels(self,
            labels,
            extraction_source_reference_attr_name="extraction_source",
            suppress_unifurcations=True,
            ):
        """
        Returns a copy of this tree that only includes leaf nodes if they are
        associated with taxon objects with labels matching those listed in
        ``labels``. Note that this copy will be a "thin" copy, including just
        the basic structure (nodes, edges) and minimal attributes (edge
        lengths, node labels, and taxon associations). Annotations,
        comments, and other attributes are not copied.

        Parameters
        ----------
        labels : iterable of str instances
            List or some other iterable of strings to match.
        suppress_unifurcations : bool
            If |True|, nodes of outdegree 1 will be deleted. Only will
            be done if some nodes are excluded from the cloned tree.
        is_apply_filter_to_leaf_nodes : bool
            If ``True`` then the above filter will be applied to leaf nodes. If
            ``False`` then it will not (and all leaf nodes will be
            automatically included, unless excluded by an ancestral node being
            filtered out).
        is_apply_filter_to_internal_nodes : bool
            If ``True`` then the above filter will be applied to internal nodes. If
            ``False`` then it will not (internal nodes without children will
            still be filtered out).

        Examples
        --------

        A clone that only extracts a subtree with taxa in the genus
        "Rhacophorus"::

            tree0 = dendropy.Tree.get(
                        path="old_world_frogs.tre",
                        schema="newick")
            # Include taxa only if label starts with "Rhacophorus"
            labels = set([taxon.label for taxon in tree0.taxon_namespace
                    if taxon.label.startswith("Rhacophorus")])
            tree1 = tree0.extract_tree_with_taxa_labels(labels=labels)

            # Above is equivalent to, but more efficient than:
            #   inclusion_set = [nd.taxon for nd in tree0.leaf_node_iter()
            #           if nd.taxon.label.startswith("Rhacophorus)]
            #   tree1 = dendropy.Tree(tree0)
            #   tree1.retain_taxa(inclusion_set)

        Returns
        -------
        t : |Tree|
            A new tree based on this one, with nodes filtered out if specified.

        """
        node_filter_fn = lambda nd: nd.taxon is None or nd.taxon.label in set(labels)
        return self.extract_tree(
                node_filter_fn=node_filter_fn,
                extraction_source_reference_attr_name=extraction_source_reference_attr_name,
                is_apply_filter_to_leaf_nodes=True,
                is_apply_filter_to_internal_nodes=False,
                )

    def extract_tree_without_taxa(self,
            taxa,
            extraction_source_reference_attr_name="extraction_source",
            suppress_unifurcations=True,
            ):
        """
        Returns a copy of this tree that only includes leaf nodes if they
        are NOT associated with the taxon objects listed in ``taxa``. Note that
        this copy will be a "thin" copy, including just the basic structure
        (nodes, edges) and minimal attributes (edge lengths, node labels, and
        taxon associations). Annotations, comments, and other attributes are
        not copied.

        Parameters
        ----------
        taxa : iterable of |Taxon| instances
            List or some other iterable of |Taxon| objects to exclude.
        suppress_unifurcations : bool
            If |True|, nodes of outdegree 1 will be deleted. Only will
            be done if some nodes are excluded from the cloned tree.
        is_apply_filter_to_leaf_nodes : bool
            If ``True`` then the above filter will be applied to leaf nodes. If
            ``False`` then it will not (and all leaf nodes will be
            automatically included, unless excluded by an ancestral node being
            filtered out).
        is_apply_filter_to_internal_nodes : bool
            If ``True`` then the above filter will be applied to internal nodes. If
            ``False`` then it will not (internal nodes without children will
            still be filtered out).

        Examples
        --------

        A clone that only extracts a subtree with taxa NOT in the genus
        "Rhacophorus"::

            tree0 = dendropy.Tree.get(
                        path="old_world_frogs.tre",
                        schema="newick")
            # Exclude taxa if their name starts with "Rhacophorus"
            taxa_to_exclude = set([taxon for taxon in tree0.taxon_namespace
                    if taxon.label.startswith("Rhacophorus")])
            tree1 = tree0.extract_tree_with_taxa(taxa=taxa_to_retain)

            # Above is equivalent to, but more efficient than:
            #   inclusion_set = [nd.taxon for nd in tree0.leaf_node_iter()
            #           if nd.taxon.label.startswith("Rhacophorus)]
            #   tree1 = dendropy.Tree(tree0)
            #   tree1.retain_taxa(inclusion_set)

        Returns
        -------
        t : |Tree|
            A new tree based on this one, with nodes filtered out if specified.

        """
        node_filter_fn = lambda nd: nd.taxon is None or nd.taxon not in set(taxa)
        return self.extract_tree(
                node_filter_fn=node_filter_fn,
                extraction_source_reference_attr_name=extraction_source_reference_attr_name,
                is_apply_filter_to_leaf_nodes=True,
                is_apply_filter_to_internal_nodes=False,
                )

    def extract_tree_without_taxa_labels(self,
            labels,
            extraction_source_reference_attr_name="extraction_source",
            suppress_unifurcations=True,
            ):
        """
        Returns a copy of this tree that only includes leaf nodes if they
        are NOT associated with the taxon objects listed in ``taxa``. Note that
        this copy will be a "thin" copy, including just the basic structure
        (nodes, edges) and minimal attributes (edge lengths, node labels, and
        taxon associations). Annotations, comments, and other attributes are
        not copied.

        Parameters
        ----------
        labels : iterable of str instances
            List or some other iterable of strings to match.
        suppress_unifurcations : bool
            If |True|, nodes of outdegree 1 will be deleted. Only will
            be done if some nodes are excluded from the cloned tree.
        is_apply_filter_to_leaf_nodes : bool
            If ``True`` then the above filter will be applied to leaf nodes. If
            ``False`` then it will not (and all leaf nodes will be
            automatically included, unless excluded by an ancestral node being
            filtered out).
        is_apply_filter_to_internal_nodes : bool
            If ``True`` then the above filter will be applied to internal nodes. If
            ``False`` then it will not (internal nodes without children will
            still be filtered out).

        Examples
        --------

        A clone that only extracts a subtree with taxa NOT in the genus
        "Rhacophorus"::

            tree0 = dendropy.Tree.get(
                        path="old_world_frogs.tre",
                        schema="newick")
            # Exclude taxa if label starts with "Rhacophorus"
            labels = set([taxon.label for taxon in tree0.taxon_namespace
                    if taxon.label.startswith("Rhacophorus")])
            tree1 = tree0.extract_tree_without_taxa_labels(labels=labels)

            # Above is equivalent to, but more efficient than:
            #   inclusion_set = [nd.taxon for nd in tree0.leaf_node_iter()
            #           if nd.taxon.label.startswith("Rhacophorus)]
            #   tree1 = dendropy.Tree(tree0)
            #   tree1.prune_taxa(inclusion_set)

        Returns
        -------
        t : |Tree|
            A new tree based on this one, with nodes filtered out if specified.

        """
        node_filter_fn = lambda nd: nd.taxon is None or nd.taxon.label not in set(labels)
        return self.extract_tree(
                node_filter_fn=node_filter_fn,
                extraction_source_reference_attr_name=extraction_source_reference_attr_name,
                is_apply_filter_to_leaf_nodes=True,
                is_apply_filter_to_internal_nodes=False,
                )

    ###########################################################################
    ### I/O

    def _format_and_write_to_stream(self, stream, schema, **kwargs):
        """
        Writes out ``self`` in ``schema`` format to a destination given by
        file-like object ``stream``.

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
        from dendropy.datamodel.treecollectionmodel import TreeList
        tree_list = TreeList(taxon_namespace=self.taxon_namespace)
        tree_list.append(self, taxon_import_strategy="add")
        # Go through TreeList.write() to reduce testing targets (i.e., testing
        # Tree.write() tests TreeList.write())
        tree_list.write_to_stream(stream, schema, **kwargs)
        # writer.write_tree_list(tree_list, stream)

    ###########################################################################
    ### Node and Edge Collection Access

    def nodes(self, filter_fn=None):
        """
        Returns list of nodes on tree.

        Parameters
        ----------

        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be included in
            the list, or |False| if not. If ``filter_fn`` is |None| (default),
            then all nodes visited will be included.

        Returns
        -------
        :py:class:`list` [|Node|]
            List of |Node| objects in the tree.
        """
        nodes = [node for node in self.preorder_node_iter(filter_fn)]
        return nodes

    def leaf_nodes(self):
        """
        Returns list of leaf nodes on the tree.

        Returns
        -------
        :py:class:`list` [|Node|]
            List of leaf |Node| objects in ``self``.
        """
        return [leaf for leaf in self.leaf_node_iter()]

    def internal_nodes(self, exclude_seed_node=False):
        """
        Returns list of internal nodes in the tree.

        Root or seed node is included unless ``exclude_seed_node`` is |True|.

        Parameters
        ----------
        exclude_seed_node : boolean, optional
            If |False| (default), then the seed node or root is included. If
            |True|, then the seed node is omitted.

        Returns
        -------
        :py:class:`list` [|Node|]
            List of internal |Node| objects in ``self``.
        """
        return [nd for nd in self.preorder_internal_node_iter(exclude_seed_node=exclude_seed_node)]

    def edges(self, filter_fn=None):
        """
        Returns list of edges on tree.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be included,
            or |False| if not. If ``filter_fn`` is |None| (default), then all
            edges will be included.

        Returns
        -------
        :py:class:`list` [|Edge|]
            List of |Edge| objects in ``self``.
        """
        edges = [edge for edge in self.preorder_edge_iter(filter_fn)]
        return edges

    def leaf_edges(self):
        """
        Returns list of leaf edges on the tree.

        Returns
        -------
        :py:class:`list` [|Edge|]
            List of leaf |Edge| objects in ``self``.
        """
        return [leaf.edge for leaf in self.leaf_node_iter()]

    def internal_edges(self, exclude_seed_edge=False):
        """
        Returns list of internal edges on tree.

        Parameters
        ----------
        exclude_seed_node : boolean, optional
            If |False| (default), then the edge subtending the seed node or
            root is included. If |True|, then the seed node is omitted.

        Returns
        -------
        :py:class:`list` [|Edge|]
            List of internal |Edge| objects in ``self``.
        """
        return [nd.edge for nd in self.preorder_internal_node_iter(exclude_seed_node=exclude_seed_edge)]

    ###########################################################################
    ### Node Finders

    def find_node(self, filter_fn):
        """
        Finds the first node for which ``filter_fn(node) == True``.

        For example, if::

            filter_fn = lambda n: hasattr(n, 'genes') and n.genes is not None

        then::

            node = t.find_node(filter_fn=filter_fn)

        will return the first node which has the attribute 'genes' and this
        value is not None.

        Parameters
        ----------
        filter_fn : function object
            Takes a single |Node| object as an argument and returns
            |True| if the node should be returned.

        Returns
        -------
        |Node| or |None|
            Returns first |Node| object for which the filter function
            ``filter_fn`` returns |True|, or |None| if no such node exists on
            this tree.
        """
        for node in self.preorder_node_iter(filter_fn):
            return node
        return None

    def find_nodes(self, filter_fn):
        """
        Find all nodes for which ``filter_fn(node) == True``.

        For example, if::

            filter_fn = lambda n: hasattr(n, 'genes') and n.genes is not None

        then::

            nodes = t.find_node(filter_fn=filter_fn)

        will return all nodes which have the attribute 'genes' and this
        value is not None.

        Parameters
        ----------
        filter_fn : function object
            Takes a single |Node| object as an argument and returns
            |True| if the node should be returned.

        Returns
        -------
        nodes : list of |Node| instances
            Returns list of |Node| objects for which the filter function
            ``filter_fn`` returns |True|.
        """
        return [node for node in self.preorder_node_iter(filter_fn)]

    def find_node_with_label(self, label):
        """
        Returns first node with ``label`` attribute matching ``label`` argument.

        Parameters
        ----------
        label : string
            Value for ``label`` attribute of |Node| object in this tree.

        Returns
        -------
        |Node| or |None|
            Returns first |Node| object with ``label`` attribute having value
            given in ``label``, or |None| if no such node is found.

        """
        for node in self.preorder_node_iter():
            if node.label == label:
                return node
        return None

    def find_node_for_taxon(self, taxon):
        """
        Returns node associated with |Taxon| object ``taxon``.

        Parameters
        ----------
        taxon : |Taxon| object
            |Taxon| object that should be associated with the node to be
            returned.

        Returns
        -------
        |Node| or |None|
            Returns first |Node| object with ``taxon`` attribute referencing same
            object as ``taxon`` argument, or |None| if no such node exists.
        """
        for node in self.postorder_node_iter():
            try:
                if node.taxon is taxon:
                    return node
            except AttributeError:
                pass
        return None

    def find_node_with_taxon(self, taxon_filter_fn=None):
        """
        Returns node associated with |Taxon| object for which ``taxon_filter_fn``
        returns |True|.

        Parameters
        ----------
        taxon_filter_fn : function object
            Takes a single |Taxon| object as an argument and returns
            |True| if the node associated with that |Taxon| should be
            returned.

        Returns
        -------
        |Node| or |None|
            Returns first |Node| object with ``taxon`` attribute passing filter
            function ``taxon_filter_fn``, or |None| if no such node is found.
        """
        for node in self.preorder_node_iter():
            if hasattr(node, "taxon") and node.taxon is not None:
                if taxon_filter_fn(node.taxon):
                    return node
        return None

    def find_node_with_taxon_label(self, label):
        """
        Returns node associated with |Taxon| object with the specified label.

        Parameters
        ----------
        label : string
            Label of |Taxon| object associated with the node to be returned.

        Returns
        -------
        |Node| or |None|
            Returns first |Node| object with ``taxon`` attribute having label
            ``label``, or|None| if no such node is found.

        """
        return self.find_node_with_taxon(lambda x: x.label == label)
        # taxon = self.taxon_namespace.get_taxon(label=label)
        # if taxon is None:
        #     return None
        # return self.find_node_with_taxon(lambda x: x is taxon)

    def mrca(self, **kwargs):
        """
        Returns most-recent common ancestor node of a set of taxa on the tree.

        Returns the shallowest node in the tree (the node nearest the tips)
        that has all of the taxa that:

            * are specified by the leafset bitmask given by the keyword argument
              ``leafset_bitmask``
            * are in the list of Taxon objects given by the keyword argument
              ``taxa``
            * have the labels specified by the list of strings given by the
              keyword argument ``taxon_labels``

        Returns |None| if no appropriate node is found. Assumes that
        bipartitions have been encoded on the tree. It is possible that the
        leafset bitmask is not compatible with the subtree that is returned!
        (compatibility tests are not fully performed).  This function is used
        to find the "insertion point" for a new bipartition via a root to tip
        search.

        Parameters
        ----------
        \*\*kwargs : keyword arguments
            Exactly one of the following must be specified:

                ``leafset_bitmask`` : integer
                    Node object subtended by the first edge compatible with this
                    leafset bitmask will be returned.
                ``taxa`` : collections.Iterable [|Taxon|]
                    Shallowest node object with descendent nodes associated with
                    all the |Taxon| objects specified will be returned.
                ``taxon_labels`` : collections.Iterable [string]
                    Shallowest node object with descendent nodes associated
                    with the minimal set of Taxon objects that
                    collectively have all the labels specified in
                    ``taxon_labels`` will be returned.

            In addition, the following optional keywords are supported:

                ``start_node`` : |Node|, optional
                    If given, specifies the node at which to start searching.
                    If not, defaults to the root or ``seed_node``.

        Returns
        -------
        |Node| or |None|
            The most-recent common ancestor of the nodes specified, or |None|
            if no such node exists.
        """
        start_node = kwargs.get("start_node", self.seed_node)
        leafset_bitmask = None
        if "leafset_bitmask" in kwargs:
            leafset_bitmask = kwargs["leafset_bitmask"]
        else:
            taxa = kwargs.get("taxa", None)
            if taxa is None:
                if "taxon_labels" in kwargs:
                    taxa = self.taxon_namespace.get_taxa(labels=kwargs["taxon_labels"])
                    if len(taxa) != len(kwargs["taxon_labels"]):
                        raise KeyError("Not all labels matched to taxa")
                else:
                    raise TypeError("Must specify one of: 'leafset_bitmask', 'taxa' or 'taxon_labels'")
            if taxa is None:
                raise ValueError("No taxa matching criteria found")
            leafset_bitmask = self.taxon_namespace.taxa_bitmask(taxa=taxa)

        if leafset_bitmask is None or leafset_bitmask == 0:
            raise ValueError("Null leafset bitmask (0)")

        if start_node.edge.bipartition.leafset_bitmask == 0 or not kwargs.get("is_bipartitions_updated", True):
            self.encode_bipartitions(suppress_unifurcations=False)

        if (start_node.edge.bipartition.leafset_bitmask & leafset_bitmask) != leafset_bitmask:
            return None

        curr_node = start_node
        last_match = start_node
        nd_source = iter(start_node.child_nodes())
        try:
            while True:
                cm = curr_node.edge.bipartition.leafset_bitmask
                cms = (cm & leafset_bitmask)
                if cms:
                    # for at least one taxon cm has 1 and bipartition has 1
                    if cms == leafset_bitmask:
                        # curr_node has all of the 1's that bipartition has
                        if cm == leafset_bitmask:
                            return curr_node
                        last_match = curr_node
                        nd_source = iter(curr_node.child_nodes())
                    else:
                        # we have reached a child that has some, but not all of the
                        #   required taxa as descendants, so we return the last_match
                        return last_match
                curr_node = next(nd_source)
        except StopIteration:
            # we shouldn't reach this if all of the descendants are properly
            #   decorated with leafset_bitmask attributes, but there may be some hacky
            #   context in which we want to allow the function to be called with
            #   leaves that have not been encoded with leafset_bitmasks.
            return last_match

    ###########################################################################
    ### Node iterators

    def __iter__(self):
        """
        Iterate over nodes on tree in pre-order.

        Example
        -------

        >>> for nd in tree:
        ...    print(nd.label)
        ...

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the internal nodes of the subtree rooted at
            this node in post-order sequence.
        """
        return self.preorder_node_iter()

    def preorder_node_iter(self, filter_fn=None):
        """
        Pre-order iterator over nodes in tree.

        Visits nodes in ``self``, with each node visited before its children.
        Nodes can optionally be filtered by ``filter_fn``: only nodes for which
        ``filter_fn`` returns |True| when called with the node as an argument are
        yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes in ``self`` in pre-order sequence.
        """
        return self.seed_node.preorder_iter(filter_fn=filter_fn)

    def preorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes in tree.

        Visits internal nodes in ``self``, with each node visited before its
        children. In DendroPy, "internal nodes" are nodes that have at least
        one child node, and thus the root or seed node is typically included
        unless ``exclude_seed_node`` is |True|. Nodes can optionally be filtered
        by ``filter_fn``: only nodes for which ``filter_fn`` returns |True| when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If |False| (default), then the seed node or root is visited. If
            |True|, then the seed node is skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the internal nodes of ``self``.
        """
        return self.seed_node.preorder_internal_node_iter(filter_fn=filter_fn,
                exclude_seed_node=exclude_seed_node)

    def postorder_node_iter(self, filter_fn=None):
        """
        Post-order iterator over nodes of tree.

        Visits self and all descendant nodes, with each node visited after its
        children. Nodes can optionally be filtered by ``filter_fn``: only nodes
        for which ``filter_fn`` returns |True| when called with the node as an
        argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the nodes in ``self`` in post-order sequence.
        """
        return self.seed_node.postorder_iter(filter_fn=filter_fn)

    def postorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes tree.

        Visits internal nodes in ``self``, with each node visited after its
        children. In DendroPy, "internal nodes" are nodes that have at least
        one child node, and thus the root or seed node is typically included
        unless ``exclude_seed_node`` is |True|. Nodes can optionally be filtered
        by ``filter_fn``: only nodes for which ``filter_fn`` returns |True| when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If |False| (default), then the seed node or root is visited. If
            |True|, then the seed node is skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the internal nodes of ``self`` in post-order
            sequence.
        """
        return self.seed_node.postorder_internal_node_iter(filter_fn=filter_fn,
                exclude_seed_node=exclude_seed_node)

    def levelorder_node_iter(self, filter_fn=None):
        """
        Level-order iteration over nodes of tree.

        Visits nodes in ``self``, with each node and other nodes at the same
        level (distance from root) visited before their children.  Nodes can
        optionally be filtered by ``filter_fn``: only nodes for which ``filter_fn``
        returns |True| when called with the node as an argument are visited.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes of ``self`` in level-order sequence.
        """
        return self.seed_node.levelorder_iter(filter_fn=filter_fn)

    def level_order_node_iter(self, filter_fn=None):
        """
        Deprecated: use :meth:`Tree.levelorder_node_iter()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'level_order_node_iter()' will no longer be supported in future releases; use 'levelorder_node_iter()' instead",
                stacklevel=3)
        return self.seed_node.levelorder_iter(filter_fn=filter_fn)

    def inorder_node_iter(self, filter_fn=None):
        """
        In-order iteration over nodes of tree.

        Visits nodes in ``self``, with each node visited in-between its children.
        Only valid for strictly-bifurcating trees. Nodes can optionally be
        filtered by ``filter_fn``: only nodes for which ``filter_fn`` returns
        |True| when called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes of ``self`` in infix or in-order sequence.
        """
        return self.seed_node.inorder_iter(filter_fn=filter_fn)

    def leaf_node_iter(self, filter_fn=None):
        """
        Iterate over all tips or leaves of tree.

        Visits all leaf or tip in ``self``. Nodes can optionally be filtered by
        ``filter_fn``: only nodes for which ``filter_fn`` returns |True| when
        called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding leaf nodes in ``self``.
        """
        return self.seed_node.leaf_iter(filter_fn=filter_fn)

    def leaf_iter(self, filter_fn=None):
        """
        Deprecated: use :meth:`Tree.leaf_node_iter()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'leaf_iter()' will no longer be supported in future releases; use 'leaf_node_iter()' instead",
                stacklevel=3)
        return self.seed_node.leaf_iter(filter_fn=filter_fn)

    def ageorder_node_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Iterator over nodes of tree in order of the age of the node (i.e., the
                time since the present).

        Iterates over nodes in order of age ('age' is as given by the ``age``
        attribute, which is usually the sum of edge lengths from tips
        to node, i.e., time since present).
        If ``include_leaves`` is |True| (default), leaves are included in the
        iteration; if ``include_leaves`` is |False|, leaves will be skipped.
        If ``descending`` is |False| (default), younger nodes will be returned
        before older ones; if |True|, older nodes will be returned before
        younger ones.

        Parameters
        ----------
        include_leaves : boolean, optional
            If |True| (default), then leaf nodes are included in the iteration.
            If |False|, then leaf nodes are skipped.
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        descending : boolean, optional
            If |False| (default), then younger nodes are visited before older
            ones. If |True|, then older nodes are visited before younger ones.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            Iterator over age-ordered sequence of nodes of ``self``.
        """
        if self.seed_node.age is None:
            self.calc_node_ages()
        return self.seed_node.ageorder_iter(include_leaves=include_leaves,
                filter_fn=filter_fn,
                descending=descending)

    def age_order_node_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Deprecated: use :meth:`Tree.ageorder_node_iter()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'age_order_node_iter()' will no longer be supported in future releases; use 'ageorder_node_iter()' instead",
                stacklevel=3)
        return self.ageorder_node_iter(include_leaves=include_leaves,
                filter_fn=filter_fn,
                descending=descending)

    def apply(self, before_fn=None, after_fn=None, leaf_fn=None):
        """
        Applies function ``before_fn`` and ``after_fn`` to all internal nodes and
        ``leaf_fn`` to all terminal nodes in subtree starting with ``self``, with
        nodes visited in pre-order.

        Given a tree with preorder sequence of nodes of
        [a,b,i,e,j,k,c,g,l,m,f,n,h,o,p,]::

                           a
                          / \
                         /   \
                        /     \
                       /       \
                      /         \
                     /           \
                    /             c
                   b             / \
                  / \           /   \
                 /   e         /     f
                /   / \       /     / \
               /   /   \     g     /   h
              /   /     \   / \   /   / \
             i   j       k l   m n   o   p


        the following order of function calls results:

            before_fn(a)
            before_fn(b)
            leaf_fn(i)
            before_fn(e)
            leaf_fn(j)
            leaf_fn(k)
            after_fn(e)
            after_fn(b)
            before_fn(c)
            before_fn(g)
            leaf_fn(l)
            leaf_fn(m)
            after_fn(g)
            before_fn(f)
            leaf_fn(n)
            before_fn(h)
            leaf_fn(o)
            leaf_fn(p)
            after_fn(h)
            after_fn(f)
            after_fn(c)
            after_fn(a)

        Parameters
        ----------
        before_fn : function object or |None|
            A function object that takes a |Node| as its argument.
        after_fn : function object or |None|
            A function object that takes a |Node| as its argument.
        leaf_fn : function object or |None|
            A function object that takes a |Node| as its argument.

        Notes
        -----
        Adapted from work by Mark T. Holder (the ``peyotl`` module of the Open
        Tree of Life Project):

            https://github.com/OpenTreeOfLife/peyotl.git

        """
        self.seed_node.apply(before_fn, after_fn, leaf_fn)

    ###########################################################################
    ### Edge iterators

    def preorder_edge_iter(self, filter_fn=None):
        """
        Pre-order iterator over nodes in tree.

        Visits nodes in ``self``, with each node visited before its children.
        Nodes can optionally be filtered by ``filter_fn``: only nodes for which
        ``filter_fn`` returns |True| when called with the node as an argument are
        yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes in ``self`` in pre-order sequence.
        """
        # NOTE: from-scratch implementation here instead of wrapping
        # `preorder_node_iter()`for efficiency
        stack = [self.seed_node._edge]
        while stack:
            edge = stack.pop()
            if filter_fn is None or filter_fn(edge):
                yield edge
            stack.extend(n._edge for n in reversed(edge._head_node._child_nodes))

    def preorder_internal_edge_iter(self, filter_fn=None, exclude_seed_edge=False):
        """
        Pre-order iterator over internal edges in tree.

        Visits internal edges in ``self``, with each edge visited before its
        children. In DendroPy, "internal edges" are edges that have at least
        one child edge, and thus the root or seed edge is typically included
        unless ``exclude_seed_edge`` is |True|. Edges can optionally be filtered
        by ``filter_fn``: only edges for which ``filter_fn`` returns |True| when
        passed the edge as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all edges visited will be yielded.
        exclude_seed_edge : boolean, optional
            If |False| (default), then the edge subtending the seed node or
            root is visited. If |True|, then this edge is skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Edge|]
            An iterator yielding the internal edges of ``self``.
        """
        # NOTE: from-scratch implementation here instead of wrapping
        # `preorder_internal_node_iter()`for efficiency
        if exclude_seed_edge:
            froot = lambda e: e._head_node._parent_node is not None
        else:
            froot = lambda e: True
        if filter_fn:
            f = lambda x: (froot(x) and x._head_node._child_nodes and filter_fn(x)) or None
        else:
            f = lambda x: (x and froot(x) and x._head_node._child_nodes) or None
        return self.preorder_edge_iter(filter_fn=f)


    def postorder_edge_iter(self, filter_fn=None):
        """
        Post-order iterator over edges of tree.

        Visits edges in ``self``, with each edge visited after its children.
        Edges can optionally be filtered by ``filter_fn``: only edges for which
        ``filter_fn`` returns |True| when called with the edge as an argument
        are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all edges visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Edge|]
            An iterator yielding the edges in ``self`` in post-order sequence.

        """
        # NOTE: custom implementation here instead of wrapping
        # `postorder_node_iter()`for efficiency

        # stack = [(self.seed_node._edge, False)]
        # while stack:
        #     edge, state = stack.pop(0)
        #     if state:
        #         if filter_fn is None or filter_fn(edge):
        #             yield edge
        #     else:
        #         stack.insert(0, (edge, True))
        #         child_edges = [(n._edge, False) for n in edge._head_node._child_nodes]
        #         child_edges.extend(stack)
        #         stack = child_edges

        ## Prefer `pop()` to `pop(0)`.
        ## Thanks to Mark T. Holder
        ## From peyotl commits: d1ffef2 + 19fdea1
        stack = [(self.seed_node._edge, False)]
        while stack:
            edge, state = stack.pop()
            if state:
                if filter_fn is None or filter_fn(edge):
                    yield edge
            else:
                stack.append((edge, True))
                stack.extend([(n._edge, False) for n in reversed(edge._head_node._child_nodes)])

    def postorder_internal_edge_iter(self, filter_fn=None, exclude_seed_edge=False):
        """
        Pre-order iterator over internal edges tree.

        Visits internal edges in ``self``, with each edge visited after its
        children. In DendroPy, "internal edges" are edges that have at least
        one child edge, and thus the root or seed edge is typically included
        unless ``exclude_seed_edge`` is |True|. Edges can optionally be filtered
        by ``filter_fn``: only edges for which ``filter_fn`` returns |True| when
        passed the edge as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all edges visited will be yielded.
        exclude_seed_edge : boolean, optional
            If |False| (default), then the seed edge or root is visited. If
            |True|, then the seed edge is skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Edge|]
            An iterator yielding the internal edges of ``self`` in post-order
            sequence.
        """
        # NOTE: from-scratch implementation here instead of wrapping
        # `preorder_internal_node_iter()`for efficiency
        if exclude_seed_edge:
            froot = lambda e: e._head_node._parent_node is not None
        else:
            froot = lambda e: True
        if filter_fn:
            f = lambda x: (froot(x) and x._head_node._child_nodes and filter_fn(x)) or None
        else:
            f = lambda x: (x and froot(x) and x._head_node._child_nodes) or None
        return self.postorder_edge_iter(filter_fn=f)

    def levelorder_edge_iter(self, filter_fn=None):
        """
        Level-order iteration over edges of tree.

        Visits edges in ``self``, with each edge and other edges at the same
        level (distance from root) visited before their children.  Edges can
        optionally be filtered by ``filter_fn``: only edges for which ``filter_fn``
        returns |True| when called with the edge as an argument are visited.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all edges visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Edge|]
            An iterator yielding edges of ``self`` in level-order sequence.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.levelorder_iter(filter_fn=f):
            yield nd.edge

    def level_order_edge_iter(self, filter_fn=None):
        """
        Deprecated: use :meth:`Tree.levelorder_edge_iter()` instead.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'level_order_edge_iter()' will no longer be supported in future releases; use 'levelorder_edge_iter()' instead",
                stacklevel=3)
        return self.levelorder_edge_iter(filter_fn=filter_fn)

    def inorder_edge_iter(self, filter_fn=None):
        """
        In-order iteration over edges of tree.

        Visits edges in ``self``, with each edge visited in-between its children.
        Only valid for strictly-bifurcating trees. Edges can optionally be
        filtered by ``filter_fn``: only edges for which ``filter_fn`` returns
        |True| when called with the edge as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all edges visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Edge|]
            An iterator yielding edges of ``self`` in infix or in-order sequence.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.inorder_iter(filter_fn=f):
            yield nd.edge

    def leaf_edge_iter(self, filter_fn=None):
        """
        Iterate over all tips or leaves of tree.

        Visits all leaf or tip in ``self``. Edges can optionally be filtered by
        ``filter_fn``: only edges for which ``filter_fn`` returns |True| when
        called with the edge as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all edges visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Edge|]
            An iterator yielding leaf edges in ``self``.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.leaf_iter(filter_fn=f):
            yield nd.edge

    ###########################################################################
    ### Taxa Management

    def reconstruct_taxon_namespace(self,
            unify_taxa_by_label=True,
            taxon_mapping_memo=None):
        if taxon_mapping_memo is None:
            taxon_mapping_memo = {}
        for node in self:
            if (node.taxon is not None
                    and (unify_taxa_by_label or node.taxon not in self.taxon_namespace)):
                t = taxon_mapping_memo.get(node.taxon, None)
                if t is None:
                    # taxon to use not given and
                    # we have not yet created a counterpart
                    if unify_taxa_by_label:
                        # this will force usage of any taxon with
                        # a label that matches the current taxon
                        t = self.taxon_namespace.require_taxon(label=node.taxon.label)
                    else:
                        # this will unconditionally create a new taxon
                        t = self.taxon_namespace.new_taxon(label=node.taxon.label)
                    taxon_mapping_memo[node.taxon] = t
                else:
                    # taxon to use is given by mapping
                    self.taxon_namespace.add_taxon(t)
                node.taxon = t

    def update_taxon_namespace(self):
        """
        All |Taxon| objects in ``self`` that are not in
        ``self.taxon_namespace`` will be added.
        """
        for nd in self:
            if nd.taxon is not None:
                self.taxon_namespace.add_taxon(nd.taxon)
        return self.taxon_namespace

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
        set[|Taxon|]
            Set of taxa associated with ``self``.
        """
        if taxa is None:
            taxa = set()
        for nd in self:
            if nd.taxon is not None:
                taxa.add(nd.taxon)
        return taxa

    def infer_taxa(self):
        """
        Creates (and returns) a new TaxonNamespace object for ``self`` populated
        with taxa from this tree.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'infer_taxa()' will no longer be supported in future releases; use 'update_taxon_namespace()' instead",
                stacklevel=3)
        taxon_namespace = taxonmodel.TaxonNamespace()
        for node in self.postorder_node_iter():
            if node.taxon is not None:
                taxon_namespace.add_taxon(node.taxon)
        self.taxon_namespace = taxon_namespace
        return taxon_namespace

    def reindex_subcomponent_taxa(self):
        """
        Remaps node taxon objects
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'reindex_subcomponent_taxa()' will no longer be supported in future releases; use 'reconstruct_taxon_namespace()' instead",
                stacklevel=3)
        for node in self.postorder_node_iter():
            t = node.taxon
            if t:
                node.taxon = self.taxon_namespace.require_taxon(label=t.label)

    def unassign_taxa(self, exclude_leaves=False, exclude_internal=False):
        """
        Strips taxon assignments from tree. If ``exclude_leaves`` is True,
        then taxa on leaves will be retained. If ``exclude_internal`` is True,
        then taxa on internal nodes will be retained. The ``taxon_namespace`` is not
        affected by this operation.
        """
        for nd in self.postorder_node_iter():
            if (len(nd._child_nodes) == 0) and not exclude_leaves:
                nd.taxon = None
            elif (len(nd._child_nodes) > 0) and not exclude_internal:
                nd.taxon = None

    def randomly_assign_taxa(self, create_required_taxa=True, rng=None):
        """
        Randomly assigns taxa to leaf nodes. If the number of taxa defined in
        the taxon set of the tree is more than the number of tips, then a random
        subset of taxa in ``taxon_namespace`` will be assigned to the tips of tree.
        If the number of tips is more than the number of taxa in the ``taxon_namespace``,
        and ``add_extra_taxa`` is not True [default], then new Taxon
        objects will be created and added to the ``taxon_namespace``; if ``create_required_taxa``
        is False, then an exception is raised.

        In addition, a Random() object or equivalent can be passed using ``rng``;
        otherwise GLOBAL_RNG is used.
        """
        if rng is None:
            rng = GLOBAL_RNG
        if len(self.taxon_namespace) == 0:
            for i, nd in enumerate(self.leaf_nodes()):
                nd.taxon = self.taxon_namespace.require_taxon(label=("T%d" % (i+1)))
        else:
            taxa = [t for t in self.taxon_namespace]
            for i, nd in enumerate(self.leaf_nodes()):
                if len(taxa) > 0:
                    nd.taxon = taxa.pop(rng.randint(0, len(taxa)-1))
                else:
                    if not create_required_taxa:
                        raise ValueError("TaxonNamespace has %d taxa, but tree has %d tips" % (len(self.taxon_namespace), len(self.leaf_nodes())))
                    label = "T%d" % (i+1)
                    k = 0
                    while self.taxon_namespace.has_taxon(label=label):
                        label = "T%d" % (i+1+k)
                        k += 1
                    nd.taxon = self.taxon_namespace.require_taxon(label=label)

    ###########################################################################
    ### Structure

    def _get_is_rootedness_undefined(self):
        return self._is_rooted is None
    is_rootedness_undefined = property(_get_is_rootedness_undefined)
    # legacy:
    rooting_state_is_undefined = property(_get_is_rootedness_undefined)

    def _get_is_rooted(self):
        return None if self._is_rooted is None else self._is_rooted
    def _set_is_rooted(self, val):
        self._is_rooted = val
    is_rooted = property(_get_is_rooted, _set_is_rooted)

    def _get_is_unrooted(self):
        return None if self._is_rooted is None else (not self._is_rooted)
    def _set_is_unrooted(self, val):
        self._is_rooted = not val
    is_unrooted = property(_get_is_unrooted, _set_is_unrooted)

    def collapse_basal_bifurcation(self, set_as_unrooted_tree=True):
        "Converts a degree-2 node at the root to a degree-3 node."
        seed_node = self.seed_node
        if not seed_node:
            return
        child_nodes = seed_node.child_nodes()
        if len(child_nodes) != 2:
            return

        if len(child_nodes[1].child_nodes()) >= 2:
            to_keep, to_del = child_nodes
        elif len(child_nodes[0].child_nodes()) >= 2:
            to_del, to_keep = child_nodes
        else:
            return
        to_del_edge = to_del.edge
        try:
            to_keep.edge.length += to_del_edge.length
        except:
            pass
        # print to_keep.edge.length, to_del_edge.length, [id(c) for c in to_del_edge.head_node.child_nodes()]
        to_del_edge.collapse(adjust_collapsed_head_children_edge_lengths=False)
        if set_as_unrooted_tree:
            self.is_rooted = False
        return self.seed_node

    def _get_seed_node(self):
        return self._seed_node
    def _set_seed_node(self, node):
        self._seed_node = node
        if self._seed_node is not None:
            self._seed_node.parent_node = None
    seed_node = property(_get_seed_node, _set_seed_node)

    def deroot(self):
        self.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    def reseed_at(self,
            new_seed_node,
            update_bipartitions=False,
            collapse_unrooted_basal_bifurcation=True,
            suppress_unifurcations=True):
        """
        Reseeds the tree at a different (existing) node.

        Takes an internal node, ``new_seed_node`` that must already be in the
        tree and rotates the tree such that ``new_seed_node`` is the ``seed_node``
        of the tree. This is a 'soft' rerooting -- i.e., changes the tree
        representation so tree traversal behaves as if the tree is rooted at
        'new_seed_node', but it does not actually change the tree's rooting
        state.  If ``update_bipartitions`` is True, then the edges'
        ``bipartition_bitmask`` and the tree's ``bipartition_edge_map`` attributes
        will be updated. If the *old* root of the tree had an outdegree of 2,
        then after this operation, it will have an outdegree of one. In this
        case, unless ``suppress_unifurcations`` is False, then it will be removed
        from the tree.
        """

        # def _dump_node(nd, name):
        #     print("- {}: {}".format(name, nd.label))
        #     if nd._parent_node:
        #         print("    Node Parent: {}".format(nd._parent_node.label))
        #     else:
        #         print("    Node Parent: None")
        #     if nd._edge.tail_node:
        #         print("    Edge Parent: {}".format(nd._edge.tail_node.label))
        #     else:
        #         print("    Edge Parent: None")
        #     debug_children = []
        #     for ch in nd._child_nodes:
        #         parts = []
        #         if ch._parent_node:
        #             parts.append(ch._parent_node.label)
        #         else:
        #             parts.append("None")
        #         if ch.edge.tail_node:
        #             parts.append(ch.edge.tail_node.label)
        #         else:
        #             parts.append("None")
        #         debug_children.append("{} ({})".format(ch.label, "/".join(parts)))
        #     debug_children = ", ".join(debug_children)
        #     print("    Children (Node Parent, Edge Tail Node Parent): {}".format(debug_children))

        if self.seed_node is new_seed_node:
            # do not just return: allow for updating of bipartitions,
            # collapsing of unifurcations, collapsing of unrooted basal
            # bifurcations
            pass
        else:
            old_seed_node = self.seed_node
            old_parent_node = new_seed_node._parent_node
            if old_parent_node is None:
                return

            if new_seed_node._child_nodes:
                new_seed_node_is_leaf = False
            else:
                new_seed_node_is_leaf = True

            edges_to_invert = []
            current_node = new_seed_node
            while current_node:
                if current_node._parent_node is not None:
                    edges_to_invert.append(current_node.edge)
                current_node = current_node._parent_node
            while edges_to_invert:
                edge = edges_to_invert.pop()
                edge.invert(update_bipartitions=update_bipartitions)

            if new_seed_node_is_leaf and suppress_unifurcations:
                ## Cannot just suppress_unifurcations, because wrong node will be deleted
                ## need to remove child (i.e. new seed node's old parent, which is now its child, needs to be deleted)
                # self.suppress_unifurcations(update_bipartitions=update_bipartitions)
                if len(new_seed_node._child_nodes) == 1:
                    nsn_ch = new_seed_node._child_nodes[0]
                    new_seed_node.remove_child(nsn_ch)
                    for ch in nsn_ch._child_nodes:
                        new_seed_node.add_child(ch)
            self.seed_node = new_seed_node

        if update_bipartitions:
            self.encode_bipartitions(
                    suppress_unifurcations=suppress_unifurcations,
                    collapse_unrooted_basal_bifurcation=collapse_unrooted_basal_bifurcation)
        else:
            if (collapse_unrooted_basal_bifurcation
                    and not self._is_rooted
                    and len(self.seed_node._child_nodes) == 2):
                self.collapse_basal_bifurcation()
            if suppress_unifurcations:
                self.suppress_unifurcations()

        return self.seed_node

    def to_outgroup_position(self, outgroup_node, update_bipartitions=False, suppress_unifurcations=True):
        """Reroots the tree at the parent of ``outgroup_node`` and makes ``outgroup_node`` the first child
        of the new root.  This is just a convenience function to make it easy
        to place a clade as the first child under the root.
        Assumes that ``outgroup_node`` and ``outgroup_node._parent_node`` and are in the tree/
        If ``update_bipartitions`` is True, then the edges' ``bipartition`` and the tree's
        ``bipartition_encoding`` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        ``suppress_unifurcations`` is False, then it will be
        removed from the tree.
        """
        p = outgroup_node._parent_node
        assert p is not None
        self.reseed_at(p, update_bipartitions=update_bipartitions, suppress_unifurcations=suppress_unifurcations)
        p.remove_child(outgroup_node)
        _ognlen = outgroup_node.edge.length
        p.insert_child(0, outgroup_node)
        assert outgroup_node.edge.length == _ognlen
        return self.seed_node

    def reroot_at_node(self, new_root_node, update_bipartitions=False, suppress_unifurcations=True):
        """
        Takes an internal node, ``new_seed_node`` that must already be in the tree and
        roots the tree at that node.
        This is a 'hard' rerooting -- i.e., changes the tree
        representation so tree traversal behaves as if the tree is rooted at
        'new_seed_node', *and* changes the tree's rooting state.
        If ``update_bipartitions`` is True, then the edges' ``bipartition`` and the tree's
        ``bipartition_encoding`` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        ``suppress_unifurcations`` is False, then it will be
        removed from the tree.
        """
        self.reseed_at(new_seed_node=new_root_node,
                update_bipartitions=False,
                suppress_unifurcations=suppress_unifurcations)
        self.is_rooted = True
        if update_bipartitions:
            self.update_bipartitions(suppress_unifurcations=suppress_unifurcations)
        return self.seed_node

    def reroot_at_edge(self,
            edge,
            length1=None,
            length2=None,
            update_bipartitions=False,
            suppress_unifurcations=True):
        """
        Takes an internal edge, ``edge``, adds a new node to it, and then roots
        the tree on the new node.
        ``length1`` and ``length2`` will be assigned to the new (sub-)edge leading
        to the old parent of the original edge, while ``length2`` will be
        assigned to the old child of the original edge.
        If ``update_bipartitions`` is True, then the edges' ``bipartition`` and the tree's
        ``bipartition_encoding`` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        ``suppress_unifurcations`` is False, then it will be
        removed from the tree.
        """
        old_tail = edge.tail_node
        old_head = edge.head_node
        new_seed_node = old_tail.new_child(edge_length=length1)
        old_tail.remove_child(old_head)
        # new_seed_node.add_child(old_head, edge_length=length2)
        new_seed_node.add_child(old_head)
        old_head.edge.length = length2
        self.reroot_at_node(new_seed_node,
                update_bipartitions=update_bipartitions,
                suppress_unifurcations=suppress_unifurcations)
        return self.seed_node

    def reroot_at_midpoint(self, update_bipartitions=False, suppress_unifurcations=True):
        """
        Reroots the tree at the the mid-point of the longest distance between
        two taxa in a tree.
        Sets the rooted flag on the tree to True.
        If ``update_bipartitions`` is True, then the edges' ``bipartition`` and the tree's
        ``bipartition_encoding`` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        ``suppress_unifurcations`` is False, then it will be
        removed from the tree.
        """
        from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix
        pdm = PhylogeneticDistanceMatrix.from_tree(self)

        ## ugly, ugly, ugly code to find two nodes that span the midpoint
        maxtax1, maxtax2 = pdm.max_pairwise_distance_taxa()
        spanning_nodes = [None, None]
        found = 0
        for nd in self.leaf_node_iter():
            for tax in (maxtax1, maxtax2):
                if nd.taxon is tax:
                    spanning_nodes[found] = nd
                    found +=1
                    break
            if found == 2:
                break
        if spanning_nodes[0].distance_from_root() < spanning_nodes[1].distance_from_root():
            n1 = spanning_nodes[1]
            n2 = spanning_nodes[0]
        else:
            n1 = spanning_nodes[0]
            n2 = spanning_nodes[1]

        plen = float(pdm.patristic_distance(maxtax1, maxtax2)) / 2
        mrca_node = pdm.mrca(n1.taxon, n2.taxon)
        #assert mrca_node is self.mrca(taxa=[n1.taxon, n2.taxon])
        #mrca_node = self.mrca(taxa=[n1.taxon, n2.taxon])
        cur_node = n1

        break_on_node = None # populated *iff* midpoint is exactly at an existing node
        target_edge = None
        head_node_edge_len = None

        # going up ...
        while cur_node is not mrca_node:
            if cur_node.edge.length > plen:
                target_edge = cur_node.edge
                head_node_edge_len = plen #cur_node.edge.length - plen
                plen = 0
                break
            elif cur_node.edge.length < plen:
                plen -= cur_node.edge.length
                cur_node = cur_node._parent_node
            else:
                break_on_node = cur_node
                break

        assert break_on_node is not None or target_edge is not None

        if break_on_node:
            self.reseed_at(break_on_node, update_bipartitions=False, suppress_unifurcations=suppress_unifurcations)
            new_seed_node = break_on_node
        else:
            tail_node_edge_len = target_edge.length - head_node_edge_len
            old_head_node = target_edge.head_node
            old_tail_node = target_edge.tail_node
            old_tail_node.remove_child(old_head_node)
            new_seed_node = Node()
            # new_seed_node.add_child(old_head_node, edge_length=head_node_edge_len)
            new_seed_node.add_child(old_head_node)
            old_head_node.edge.length = head_node_edge_len
            # old_tail_node.add_child(new_seed_node, edge_length=tail_node_edge_len)
            old_tail_node.add_child(new_seed_node)
            new_seed_node.edge.length = tail_node_edge_len
            self.reseed_at(new_seed_node, update_bipartitions=False, suppress_unifurcations=suppress_unifurcations)
        self.is_rooted = True
        if update_bipartitions:
            self.update_bipartitions(suppress_unifurcations=False)
        return self.seed_node

    def suppress_unifurcations(self, update_bipartitions=False):
        """
        Delete all nodes of outdegree-one from this tree.

        Parameters
        ----------
        update_bipartitions : bool
            If |True| then the bipartitions encoding will be calculated.

        """
        if update_bipartitions and self.bipartition_encoding:
            bipartitions_to_delete = set()
        else:
            bipartitions_to_delete = None
        remapped_nodes = []
        for nd in self.postorder_node_iter():
            children = nd._child_nodes
            if len(children) == 1:
                remapped_nodes.append((nd, children[0]))
                if nd.edge.length is not None:
                    if children[0].edge.length is None:
                        children[0].edge.length = nd.edge.length
                    else:
                        children[0].edge.length += nd.edge.length
                if bipartitions_to_delete is not None:
                    bipartitions_to_delete.add(id(nd.edge.bipartition))
                if nd._parent_node is not None:
                    parent = nd._parent_node
                    pos = parent._child_nodes.index(nd)
                    parent.remove_child(nd)
                    parent.insert_child(index=pos, node=children[0])
                    # assert children[0]._parent_node is parent
                    # assert children[0] in parent._child_nodes
                    # assert children[0].edge.tail_node is parent
                    # assert children[0].edge.head_node is children[0]
                    nd._parent_node = None
                else:
                    # assert nd is self.seed_node
                    self.seed_node = children[0]
                    self.seed_node._parent_node = None
        if bipartitions_to_delete:
            old_encoding = self.bipartition_encoding
            self.bipartition_encoding = [b for b in old_encoding if id(b) not in bipartitions_to_delete]
        return remapped_nodes

    def delete_outdegree_one_nodes(self):
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'delete_outdegree_one_nodes()' has been replaced by 'suppress_unifurcations()'",
                stacklevel=3)
        return self.suppress_unifurcations()

    def collapse_unweighted_edges(self,
            threshold=0.0000001,
            update_bipartitions=False):
        """
        Collapse all *internal* edges with edge lengths less than or equal to
        ``threshold`` (or with |None| for edge length).
        """
        for e in self.postorder_edge_iter():
            if e.length is None or (e.length <= threshold) and e.is_internal():
               e.collapse()
        if update_bipartitions:
            self.update_bipartitions()

    def resolve_polytomies(self,
            limit=2,
            update_bipartitions=False,
            rng=None):
        """
        Arbitrarily resolve polytomies using 0-length edges.

        Parameters
        ----------
        limit : int
            The maximum number of children a node can have before being
            resolved.
        update_bipartitions : bool
            If |True|, then bipartitions will be calculated.
        rng : ``random.Random`` object or |None|
            If ``rng`` is an object with a ``sample()`` method then the polytomy
            will be resolved by sequentially adding, generating all tree
            topologies equiprobably. ``rng.sample()`` should behave like
            ``random.sample()``
            If ``rng`` is |None|, then polytomy is broken deterministically by
            repeatedly joining pairs of children.
        """
        polytomies = []
        for node in self.postorder_node_iter():
            if len(node._child_nodes) > limit:
                polytomies.append(node)
        for node in polytomies:
            if rng:
                to_attach = rng.sample(node._child_nodes, len(node._child_nodes)-limit)
                for child in to_attach:
                    node.remove_child(child)
                attachment_points = list(node._child_nodes)
                while len(to_attach) > 0:
                    next_child = to_attach.pop()
                    next_sib = rng.choice(attachment_points)
                    next_attachment = Node()
                    p = next_sib._parent_node
                    p.add_child(next_attachment)
                    next_attachment.edge.length = 0.0
                    p.remove_child(next_sib)
                    next_attachment.add_child(next_sib)
                    next_attachment.add_child(next_child)
                    attachment_points.append(next_attachment)
                    attachment_points.append(next_child)
            else:
                while len(node._child_nodes) > limit:
                    nn1 = Node()
                    nn1.edge.length = 0.0
                    c1 = node._child_nodes[0]
                    c2 = node._child_nodes[1]
                    node.remove_child(c1)
                    node.remove_child(c2)
                    nn1.add_child(c1)
                    nn1.add_child(c2)
                    node.add_child(nn1)
        if update_bipartitions:
            self.update_bipartitions()

    def prune_subtree(self,
            node,
            update_bipartitions=False,
            suppress_unifurcations=True):
        """
        Removes subtree starting at ``node`` from tree.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        if node._parent_node is None:
            raise TypeError('Node has no parent and is implicit root: cannot be pruned')
        node._parent_node.remove_child(node)
        if suppress_unifurcations:
            self.suppress_unifurcations()
        if update_bipartitions:
            self.update_bipartitions()

    def filter_leaf_nodes(
            self,
            filter_fn,
            recursive=True,
            update_bipartitions=False,
            suppress_unifurcations=True):
        """
        Removes all leaves for which ``filter_fn`` returns |False|. If recursive
        is |True|, then process is repeated until all leaf nodes in the tree will
        evaluate to |True| when passed to ``filter_fn``.

        Parameters
        ----------
        ``filter_fn`` : function object
            A function that takes a |Node| object and returns |True| if
            the object is to be allowed as a leaf node, and |False| if otherwise.
        recursive : bool
            If |True|, then filter is repeatedly applied until all leaf nodes
            evaluate to |True| under ``filter_fn``. If |False|, then only a
            single pass is made on the current leaf set. This may result in new
            leaves for which the ``filter_fn`` is |False| (e.g., the parent node
            of a cherry in which both children evaluated to |False|
            under ``filter_fn`` now is a leaf node which may be |False|
            under ``filter_fn``).
        suppress_unifurcations : bool
            If |True|, nodes of outdegree 1 will be deleted as they are
            encountered.
        update_bipartitions : bool
            If |True|, then bipartitions will be calculated.

        Returns
        -------
        nds : list[|Node|]
            List of nodes removed.
        """
        nodes_removed = []
        while True:
            is_nodes_deleted = False
            nodes_to_remove = [nd for nd in self.leaf_node_iter() if not filter_fn(nd)]
            for nd in nodes_to_remove:
                if nd.edge.tail_node is None:
                    raise error.SeedNodeDeletionException("Attempting to remove seed node or node without parent")
                nd.edge.tail_node.remove_child(nd)
            if nodes_to_remove:
                nodes_removed += nodes_to_remove
                is_nodes_deleted = True
            if not is_nodes_deleted or not recursive:
                break
        if suppress_unifurcations:
            self.suppress_unifurcations()
        if update_bipartitions:
            self.update_bipartitions()
        return nodes_removed

    def prune_leaves_without_taxa(self,
            recursive=True,
            update_bipartitions=False,
            suppress_unifurcations=True):
        """
        Removes all terminal nodes that have their ``taxon`` attribute set to
        |None|.
        """
        nodes_removed = []
        while True:
            nodes_to_remove = []
            for nd in self.leaf_node_iter():
                if nd.taxon is None:
                    nodes_to_remove.append(nd)
            for nd in nodes_to_remove:
                nd.edge.tail_node.remove_child(nd)
            nodes_removed += nodes_to_remove
            if not nodes_to_remove or not recursive:
                break
        if suppress_unifurcations:
            self.suppress_unifurcations()
        if update_bipartitions:
            self.update_bipartitions()
        return nodes_removed

    def prune_nodes(self, nodes, prune_leaves_without_taxa=False, update_bipartitions=False, suppress_unifurcations=True):
        for nd in nodes:
            if nd.edge.tail_node is None:
                raise Exception("Attempting to remove root node or node without parent")
            nd.edge.tail_node.remove_child(nd)
        if prune_leaves_without_taxa:
            self.prune_leaves_without_taxa(update_bipartitions=update_bipartitions,
                    suppress_unifurcations=suppress_unifurcations)

    def prune_taxa(self,
            taxa,
            update_bipartitions=False,
            suppress_unifurcations=True,
            is_apply_filter_to_leaf_nodes=True,
            is_apply_filter_to_internal_nodes=False):
        """
        Removes terminal nodes associated with Taxon objects given by the container
        ``taxa`` (which can be any iterable, including a TaxonNamespace object) from ``self``.
        """
        taxa = set(taxa)
        nodes_to_remove = []
        for nd in self.postorder_node_iter():
            if (
                ((is_apply_filter_to_internal_nodes and nd._child_nodes)
                or (is_apply_filter_to_leaf_nodes and not nd._child_nodes))
                and (nd.taxon and nd.taxon in taxa)
                ):
                    nd.edge.tail_node.remove_child(nd)
        self.prune_leaves_without_taxa(update_bipartitions=update_bipartitions,
                suppress_unifurcations=suppress_unifurcations)

    def prune_taxa_with_labels(self,
            labels,
            update_bipartitions=False,
            suppress_unifurcations=True,
            is_apply_filter_to_leaf_nodes=True,
            is_apply_filter_to_internal_nodes=False):
        """
        Removes terminal nodes that are associated with Taxon objects with
        labels given by ``labels``.
        """
        taxa = self.taxon_namespace.get_taxa(labels=labels)
        self.prune_taxa(taxa=taxa,
                update_bipartitions=update_bipartitions,
                suppress_unifurcations=suppress_unifurcations,
                is_apply_filter_to_leaf_nodes=is_apply_filter_to_leaf_nodes,
                is_apply_filter_to_internal_nodes=is_apply_filter_to_internal_nodes)

    def retain_taxa(self,
            taxa,
            update_bipartitions=False,
            suppress_unifurcations=True):
        """
        Removes terminal nodes that are not associated with any
        of the Taxon objects given by ``taxa`` (which can be any iterable, including a
        TaxonNamespace object) from the ``self``.
        """
        to_prune = [t for t in self.taxon_namespace if t not in taxa]
        self.prune_taxa(to_prune,
                update_bipartitions=update_bipartitions,
                suppress_unifurcations=suppress_unifurcations)

    def retain_taxa_with_labels(self,
            labels,
            update_bipartitions=False,
            suppress_unifurcations=True):
        """
        Removes terminal nodes that are not associated with Taxon objects with
        labels given by ``labels``.
        """
        taxa = self.taxon_namespace.get_taxa(labels=labels)
        self.retain_taxa(taxa=taxa,
                update_bipartitions=update_bipartitions,
                suppress_unifurcations=suppress_unifurcations)

    def randomly_reorient(self, rng=None, update_bipartitions=False):
        """
        Randomly picks a new rooting position and rotates the branches around all
        internal nodes in the ``self``. If ``update_bipartitions`` is True, the the ``bipartition_bitmask``
        and ``bipartition_edge_map`` attributes kept valid.
        """
        if rng is None:
            rng = GLOBAL_RNG # use the global rng by default
        nd = rng.sample(self.nodes(), 1)[0]
        if nd.is_leaf():
            self.to_outgroup_position(nd, update_bipartitions=update_bipartitions)
        else:
            self.reseed_at(nd, update_bipartitions=update_bipartitions)
        self.randomly_rotate(rng=rng)

    def randomly_rotate(self, rng=None):
        "Randomly rotates the branches around all internal nodes in ``self``"
        if rng is None:
            rng = GLOBAL_RNG # use the global rng by default
        internal_nodes = self.internal_nodes()
        for nd in internal_nodes:
            c = nd.child_nodes()
            rng.shuffle(c)
            nd.set_child_nodes(c)

    def shuffle_taxa(self, include_internal_nodes=False, rng=None):
        """
        Randomly re-assigns taxa associated with nodes. Note that in the case
        of not all nodes being associated with taxa, this will NOT assign taxa
        to nodes that currently do not have them, nor will nodes currently
        associated with taxa end up not being associated with taxa.
        Returns a dictionary mapping the old taxa to their new counterparts.
        """
        if rng is None:
            rng = GLOBAL_RNG # use the global rng by default
        if include_internal_nodes:
            nd_iterator = self.preorder_node_iter
        else:
            nd_iterator = self.leaf_node_iter
        current_node_taxon_map = {}
        node_taxa = set()
        for nd in nd_iterator():
            if nd.taxon is not None:
                current_node_taxon_map[nd] = nd.taxon
                assert nd.taxon not in node_taxa
                node_taxa.add(nd.taxon)
        assert len(current_node_taxon_map) == len(node_taxa)
        current_to_shuffled_taxon_map = {}
        for nd in current_node_taxon_map:
            new_taxon = rng.sample(node_taxa, 1)[0]
            current_to_shuffled_taxon_map[nd.taxon] = new_taxon
            nd.taxon = new_taxon
            node_taxa.remove(new_taxon)
        assert len(node_taxa) == 0, node_taxa
        assert len(current_to_shuffled_taxon_map) == len(current_node_taxon_map)
        return current_to_shuffled_taxon_map

    def ladderize(self, ascending=True):
        """
        Sorts child nodes in ascending (if ``ascending`` is |False|) or
        descending (if ``ascending`` is |False|) order in terms of the number of
        children each child node has.
        """
        node_desc_counts = {}
        for nd in self.postorder_node_iter():
            if len(nd._child_nodes) == 0:
                node_desc_counts[nd] = 0
            else:
                total = 0
                for child in nd._child_nodes:
                    total += node_desc_counts[child]
                total += len(nd._child_nodes)
                node_desc_counts[nd] = total
                nd._child_nodes.sort(key=lambda n: node_desc_counts[n], reverse=not ascending)

    def truncate_from_root(self, distance_from_root):
        self.calc_node_root_distances()
        new_terminals = []
        for nd in self.preorder_node_iter():
            if not nd._parent_node:
                # root node
                # TODO: strictly speaking, this might be a terminal if distance_from_root == 0
                pass
            else:
                if nd.root_distance == distance_from_root:
                    new_terminals.append(nd)
                elif nd.root_distance > distance_from_root and nd._parent_node.root_distance < distance_from_root:
                    # cut above current node
                    nd.edge.length = distance_from_root - nd._parent_node.root_distance
                    nd.root_distance = distance_from_root
                    new_terminals.append(nd)
        for nd in new_terminals:
            for ch in nd.child_nodes():
                nd.remove_child(ch)

    ###########################################################################
    ### Ages, depths, branch lengths etc. (mutation)

    def scale_edges(self, edge_len_multiplier):
        """Multiplies every edge length in ``self`` by ``edge_len_multiplier``"""
        for e in self.postorder_edge_iter():
            if e.length is not None:
                e.length *= edge_len_multiplier

    def set_edge_lengths_from_node_ages(self,
            minimum_edge_length=0.0,
            error_on_negative_edge_lengths=False,
            ):
        """
        Sets the edge lengths of the tree so that the path lengths from the
        tips equal the value of the ``age`` attribute of the nodes.

        Parameters
        ----------
        minimum_edge_length : numeric
            All edge lengths calculated to have a value less than this will be
            set to this.
        error_on_negative_edge_lengths : bool
            If |True|, an inferred edge length that is less than 0 will result
            in a ValueError.
        """
        for nd in self.preorder_node_iter():
            if nd._parent_node is not None:
                #if nd._parent_node.age < nd.age:
                #    nd.edge.length = 0.0
                #else:
                #    nd.edge.length = nd._parent_node.age - nd.age
                edge_length = nd._parent_node.age - nd.age
                if minimum_edge_length is not None and edge_length < minimum_edge_length:
                    edge_length = minimum_edge_length
                if error_on_negative_edge_lengths and edge_length < 0.0:
                    raise ValueError("Negative edge length: {}".foramt(edge_length))
                nd.edge.length = edge_length

    ###########################################################################
    ### Ages, depths, branch lengths etc. (calculation)

    def phylogenetic_distance_matrix(self):
        """
        Returns a |PhylogeneticDistanceMatrix| instance based
        on the tree (in its current state).

        Returns
        -------
        pdc : a |PhylogeneticDistanceMatrix| instance
            A |PhylogeneticDistanceMatrix| instance corresponding to the
            tree in its current state.
        """
        from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix
        return PhylogeneticDistanceMatrix.from_tree(tree=self)

    def node_distance_matrix(self):
        from dendropy.calculate.phylogeneticdistance import NodeDistanceMatrix
        return NodeDistanceMatrix.from_tree(tree=self)

    def calc_node_ages(self,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            is_force_max_age=False,
            is_force_min_age=False,
            set_node_age_fn=None,
            is_return_internal_node_ages_only=False):
        """
        Adds an attribute called "age" to  each node, with the value equal to
        the sum of edge lengths from the node to the tips.

        Parameters
        ----------
        ultrametricity_precision : numeric or bool or None
            If the lengths of different paths to the node differ by more than
            ``ultrametricity_precision``, then a ValueError exception will be
            raised indicating deviation from ultrametricity. If
            ``ultrametricity_precision`` is negative or False, then this check
            will be skipped.
        is_force_max_age: bool
            If ``is_force_max_age`` is |True|, then each node will be set to the
            maximum possible age, by being set to the oldest age given its
            child set and the subtending edge lengths. This option only makes a
            difference if the tree is not ultrametric, and so the
            ultrametricity precision check is ignore if this option is set to
            True.
        is_force_min_age: bool
            If ``is_force_min_age`` is |True| then each node will be set to the
            minimum possible age, by being set to the youngest age given its
            child set and the subtending edge lengths. This option only makes a
            difference if the tree is not ultrametric, and so the
            ultrametricity precision check is ignore if this option is set to
            True.
        set_node_age_fn: function object
            If not |None|, then this should be a function that takes a node as
            an argument and returns |None| or a non-|None| value. If
            |None|, then this indicates that the node's age should be
            calculated by this function. If not |None|, then this is the
            value that this node's age should be set to. This can be used to
            set non-contemporary tip ages by passing something like:

                f = lambda nd: None if not nd.is_leaf else nd.annotations["height"]

            which returns |None| if the node is an internal node, but
            otherwise returns the value in the ``height`` annotation.

        Returns
        -------
        a : iterable[numeric]
            Returns collection of node ages.

        """
        ages = []
        if is_force_max_age and is_force_min_age:
            raise ValueError("Cannot specify both 'is_force_max_age' and 'is_force_min_age'")
        for node in self.postorder_node_iter():
            child_nodes = node.child_nodes()
            if set_node_age_fn is not None:
                node.age = set_node_age_fn(node)
                # print("Setting node age: {} = {}".format(node.taxon, node.age))
                if node.age is not None:
                    continue
            if len(child_nodes) == 0:
               node.age = 0.0
               if not is_return_internal_node_ages_only:
                   ages.append(node.age)
            else:
                if is_force_max_age:
                    age_to_set = max([ (child.age + child.edge.length) for child in child_nodes ])
                elif is_force_min_age:
                    age_to_set = min([ (child.age + child.edge.length) for child in child_nodes ])
                else:
                    first_child = child_nodes[0]
                    if first_child.edge.length is not None and first_child.age is not None:
                        age_to_set = first_child.age + first_child.edge.length
                    elif first_child.edge.length is None:
                        first_child.edge.length = 0.0
                        age_to_set = first_child.age
                    elif first_child.age is None:
                        first_child.age = 0.0
                        age_to_set = first_child.edge.length
                    else:
                        age_to_set = 0.0
                node.age = age_to_set
                if not (is_force_max_age or is_force_min_age or ultrametricity_precision is None or ultrametricity_precision is False or ultrametricity_precision < 0):
                    for nnd in child_nodes[1:]:
                        try:
                            ocnd = nnd.age + nnd.edge.length
                        except TypeError:
                            nnd.edge.length = 0.0
                            ocnd = nnd.age
                        d = abs(node.age - ocnd)
                        if  d > ultrametricity_precision:
                            # try:
                            #     self.encode_bipartitions()
                            #     node_id = nnd.bipartition.split_as_newick_string(taxon_namespace=self.taxon_namespace)
                            # except OSError:
                            #     node_id = str(nnd)
                            node_id = str(node)
                            subtree = node._as_newick_string()
                            desc = []
                            for desc_nd in child_nodes:
                                desc.append("-   {}: has age of {} and edge length of {}, resulting in parent node age of {}".format(
                                    desc_nd,
                                    desc_nd.age,
                                    desc_nd.edge.length,
                                    desc_nd.edge.length + desc_nd.age))
                            desc = "\n".join(desc)
                            raise error.UltrametricityError(
                                    ("Tree is not ultrametric within threshold of {threshold}: {deviance}.\n"
                                     "Encountered in subtree of node {node} (edge length of {length}):\n"
                                     "\n    {subtree}\n\n"
                                     "Age of children:\n"
                                     "{desc}"
                                     ).format(
                                threshold=ultrametricity_precision,
                                deviance=d,
                                node=node_id,
                                length=node.edge.length,
                                desc=desc,
                                subtree=subtree,
                                ))
                ages.append(node.age)
        return ages

    def calc_node_root_distances(self, return_leaf_distances_only=True):
        """
        Adds attribute "root_distance" to each node, with value set to the
        sum of edge lengths from the node to the root. Returns list of
        distances. If ``return_leaf_distances_only`` is True, then only
        leaf distances will be true.
        """
        dists = []
        for node in self.preorder_node_iter():
            if node._parent_node is None:
                node.root_distance = 0.0
            else:
                node.root_distance = node.edge.length + node._parent_node.root_distance
            if (not return_leaf_distances_only or node.is_leaf()):
                dists.append(node.root_distance)
        return dists

    def internal_node_ages(self,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            is_force_max_age=False,
            is_force_min_age=False,
            set_node_age_fn=None,
            ):
        """
        Returns list of ages of speciation events / coalescence times on tree.
        """
        ages = self.calc_node_ages(
                ultrametricity_precision=ultrametricity_precision, is_return_internal_node_ages_only=True,
                is_force_max_age=is_force_max_age,
                is_force_min_age=is_force_min_age,
                set_node_age_fn=set_node_age_fn,
                )
        ages.sort()
        return ages

    def node_ages(self,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            is_force_max_age=False,
            is_force_min_age=False,
            set_node_age_fn=None,
            internal_only=False):
        """
        Returns list of ages of all nodes on tree.
        NOTE: Changed from DendroPy3: this function now returns the ages of
        *ALL* nodes. To get only internal node ages, use
        `Tree.internal_node_ages`.
        """
        ages = self.calc_node_ages(
                ultrametricity_precision=ultrametricity_precision,
                is_force_max_age=is_force_max_age,
                is_force_min_age=is_force_min_age,
                set_node_age_fn=set_node_age_fn,
                is_return_internal_node_ages_only=internal_only)
        ages.sort()
        return ages

    def length(self):
        """
        Returns sum of edge lengths of self. Edges with no lengths defined
        (None) will be considered to have a length of 0.
        Note that we do not overrride ``__len__`` as this requires an integer
        return value.
        """
        total = 0
        for edge in self.postorder_edge_iter():
            if edge.length is not None:
                total += edge.length
        return total

    def max_distance_from_root(self):
        """
        Returns distance of node furthest from root.
        """
        dists = self.calc_node_root_distances()
        return max(dists)

    def minmax_leaf_distance_from_root(self):
        """
        Returns pair of values, representing the distance of the leaf closest
        to a furthest from the root.
        """
        dists = self.calc_node_root_distances(return_leaf_distances_only=True)
        return min(dists), max(dists)

    def coalescence_intervals(self):
        """
        Returns list of coalescence intervals of self., i.e., the waiting
        times between successive coalescence events.
        """
        ages = self.node_ages()
        intervals = []
        intervals.append(ages[0])
        for i, d in enumerate(ages[1:]):
            intervals.append(d - ages[i])
        return intervals

    def num_lineages_at(self, distance_from_root):
        """
        Returns the number of lineages on the tree at a particular distance
        from the root.
        """
        self.calc_node_root_distances()
        num_lineages = 0
        for nd in self.preorder_node_iter():
            if not nd._parent_node:
                # root node
                pass
            else:
                if nd.root_distance == distance_from_root:
                    num_lineages += 1
                elif nd.root_distance >= distance_from_root and nd._parent_node.root_distance < distance_from_root:
                    num_lineages += 1
        return num_lineages

    ###########################################################################
    ### Bipartition Management

    def _compile_mutable_bipartition_for_edge(self, edge):
        edge.bipartition.compile_split_bitmask(
                tree_leafset_bitmask=self.seed_node.edge.bipartition._leafset_bitmask,
                is_mutable=True)
        return edge.bipartition

    def _compile_immutable_bipartition_for_edge(self, edge):
        edge.bipartition.compile_split_bitmask(
                tree_leafset_bitmask=self.seed_node.edge.bipartition._leafset_bitmask,
                is_mutable=False)
        return edge.bipartition

    def encode_bipartitions(self,
            suppress_unifurcations=True,
            collapse_unrooted_basal_bifurcation=True,
            suppress_storage=False,
            is_bipartitions_mutable=False):
        """
        Calculates the bipartitions of this tree.

        Parameters
        ----------
        suppress_unifurcations : bool
            If |True|, nodes of outdegree 1 will be deleted as they are
            encountered.
        collapse_unrooted_basal_bifurcation: bool
            If |True|, then a basal bifurcation on an unrooted tree will be
            collapsed to a trifurcation. This mean that an unrooted tree like
            '(A,(B,C))' will be changed to '(A,B,C)' after this.
        suppress_storage : bool
            By default, the bipartition encoding is stored as a list (assigned
            to ``self.bipartition_encoding``) and returned. If ``suppress_storage``
            is |True|, then the list is not created.
        is_bipartitions_mutable : bool
            By default, the |Bipartition| instances coded will be locked
            or frozen, allowing their use in hashing containers such as
            dictionary (keys) and sets. To allow modification of values, the
            ``is_mutable`` attribute must be set to |True|.

        Returns
        -------
        list[|Bipartition|] or |None|
            A list of |Bipartition| objects of this |Tree|
            representing the structure of this tree, or, if ``suppress_storage``
            is |True|, then |None|.

        """
        self._bipartition_edge_map = None
        taxon_namespace = self._taxon_namespace
        seed_node = self.seed_node
        if not seed_node:
            return
        if (collapse_unrooted_basal_bifurcation
                and not self._is_rooted
                and len(seed_node._child_nodes) == 2):
            # We do this because an unrooted tree
            # has no *true* degree-3 internal nodes:
            #
            #      \  | |  /
            #       +-+-+-+
            #      /       \
            #
            # (whereas, with a rooted tree, the basal bipartition is a true
            # degree-3 node: the edge subtending it does not really
            # exist in the graph -- it is not a true link connecting
            # two nodes).
            self.collapse_basal_bifurcation()
        tree_edges = []
        for edge in self.postorder_edge_iter():
            leafset_bitmask = 0
            head_node = edge._head_node
            child_nodes = head_node._child_nodes
            num_children = len(child_nodes)
            if num_children == 1 and suppress_unifurcations:
                # collapsing node: remove, and do not process/add edge
                if head_node.edge.length is not None:
                    if child_nodes[0].edge.length is None:
                        child_nodes[0].edge.length = head_node.edge.length
                    else:
                        child_nodes[0].edge.length += head_node.edge.length
                if head_node._parent_node is not None:
                    parent = head_node._parent_node
                    pos = parent._child_nodes.index(head_node)
                    parent.remove_child(head_node)
                    parent.insert_child(index=pos, node=child_nodes[0])
                    head_node._parent_node = None
                else:
                    self.seed_node = child_nodes[0]
                    self.seed_node._parent_node = None
            else:
                if num_children == 0:
                    tree_edges.append(edge)
                    taxon = edge._head_node.taxon
                    if taxon:
                        leafset_bitmask = taxon_namespace.taxon_bitmask(taxon)
                else:
                    tree_edges.append(edge)
                    for child in child_nodes:
                        leafset_bitmask |= child.edge.bipartition._leafset_bitmask
                edge.bipartition = Bipartition(compile_bipartition=False, is_mutable=True)
                edge.bipartition._leafset_bitmask = leafset_bitmask
                edge.bipartition._is_rooted = self._is_rooted
        # Create normalized bitmasks, where the full (self) bipartition mask is *not*
        # all the taxa, but only those found on the self; this is to handle
        # cases where we are dealing with selfs with incomplete leaf-sets.
        tree_leafset_bitmask = self.seed_node.edge.bipartition._leafset_bitmask
        if is_bipartitions_mutable:
            _compile_bipartition = self._compile_mutable_bipartition_for_edge
        else:
            _compile_bipartition = self._compile_immutable_bipartition_for_edge
        if suppress_storage:
            self.bipartition_encoding = None
            for x in map(_compile_bipartition, tree_edges):
                pass
        else:
            # self.bipartition_encoding = dict(zip(map(self._compile_bipartition_for_edge, tree_edges), tree_edges))
            self.bipartition_encoding = list(map(_compile_bipartition, tree_edges))
        return self.bipartition_encoding

    def update_bipartitions(self, *args, **kwargs):
        """
        Recalculates bipartition hashes for tree.
        """
        self.encode_bipartitions(*args, **kwargs)

    def encode_splits(self, *args, **kwargs):
        """
        Recalculates bipartition hashes for tree.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'Tree.encode_splits()' will no longer be supported in future releases; use 'Tree.encode_bipartitions()' instead",
                stacklevel=3)
        return self.encode_bipartitions(*args, **kwargs)

    def update_splits(self, *args, **kwargs):
        """
        Recalculates bipartition hashes for tree.
        """
        deprecate.dendropy_deprecation_warning(
                message="Deprecated since DendroPy 4: 'Tree.encode_splits()' will no longer be supported in future releases; use 'Tree.update_bipartitions()' instead",
                stacklevel=3)
        return self.encode_bipartitions(*args, **kwargs)

    def _get_bipartition_edge_map(self):
        if not self._bipartition_edge_map:
            if not self.bipartition_encoding:
                self.encode_bipartitions()
            self._bipartition_edge_map = {}
            for edge in self.postorder_edge_iter():
                self._bipartition_edge_map[edge.bipartition] = edge
        return self._bipartition_edge_map
    bipartition_edge_map = property(_get_bipartition_edge_map)

    ###########################################################################
    ### Metrics -- Unary

    def B1(self):
        """DEPRECATED: Use :func:`dendropy.calculate.treemeasure.B1()`."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Unary statistics on trees are now implemented in the 'dendropy.calculate.treemeasure' module.",
                old_construct="tree.B1()",
                new_construct="from dendropy.calculate import treemeasure\ntreemeasure.B1(tree)")
        from dendropy.calculate import treemeasure
        return treemeasure.B1(self)

    def colless_tree_imbalance(self, normalize="max"):
        """DEPRECATED: Use 'dendropy.calculate.treemeasure.colless_tree_imbalance()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Unary statistics on trees are now implemented in the 'dendropy.calculate.treemeasure' module.",
                old_construct="tree.colless_tree_imbalance()",
                new_construct="from dendropy.calculate import treemeasure\ntreemeasure.colless_tree_imbalance(tree)")
        from dendropy.calculate import treemeasure
        return treemeasure.colless_tree_imbalance(self, normalize)

    def pybus_harvey_gamma(self, prec=0.00001):
        """DEPRECATED: Use 'dendropy.calculate.treemeasure.pybus_harvey_gamma()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Unary statistics on trees are now implemented in the 'dendropy.calculate.treemeasure' module.",
                old_construct="tree.pybus_harvey_gamma()",
                new_construct="from dendropy.calculate import treemeasure\ntreemeasure.pybus_harvey_gamma(tree)")
        from dendropy.calculate import treemeasure
        return treemeasure.pybus_harvey_gamma(self, prec)

    def N_bar(self):
        """DEPRECATED: Use 'dendropy.calculate.treemeasure.N_bar()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Unary statistics on trees are now implemented in the 'dendropy.calculate.treemeasure' module.",
                old_construct="tree.N_bar()",
                new_construct="from dendropy.calculate import treemeasure\ntreemeasure.N_bar(tree)")
        from dendropy.calculate import treemeasure
        return treemeasure.N_bar(self)

    def sackin_index(self, normalize=True):
        """DEPRECATED: Use 'dendropy.calculate.treemeasure.sackin_index()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Unary statistics on trees are now implemented in the 'dendropy.calculate.treemeasure' module.",
                old_construct="tree.sackin_index()",
                new_construct="from dendropy.calculate import treemeasure\ntreemeasure.sackin_index(tree)")
        from dendropy.calculate import treemeasure
        return treemeasure.sackin_index(self, normalize)

    def treeness(self):
        """DEPRECATED: Use 'dendropy.calculate.treemeasure.treeness()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Unary statistics on trees are now implemented in the 'dendropy.calculate.treemeasure' module.",
                old_construct="tree.treeness()",
                new_construct="from dendropy.calculate import treemeasure\ntreemeasure.treeness(tree)")
        from dendropy.calculate import treemeasure
        return treemeasure.treeness(self)

    ###########################################################################
    ### Comparisons with Another Tree

    def is_compatible_with_bipartition(self, bipartition, is_bipartitions_updated=False):
        """
        Returns true if the |Bipartition| ``bipartition`` is compatible
        with all bipartitions of this tree.
        """
        if not is_bipartitions_updated or not self.bipartition_encoding:
            self.encode_bipartitions()
        if bipartition in self.bipartition_encoding:
            return True
        else:
            for b in self.bipartition_encoding:
                if not b.is_compatible_with(bipartition):
                    return False
            return True

    def is_compatible_with_tree(self, other):
        raise NotImplementedError

    def find_missing_splits(self, other_tree):
        """DEPRECATED: Use 'dendropy.treecompare.find_missing_bipartitions()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Statistics comparing two trees are now implemented in the 'dendropy.calculate.treecompare' module.",
                old_construct="tree1.find_missing_splits(tree2)",
                new_construct="from dendropy.calculate import treecompare\ntreecompare.find_missing_bipartitions(tree1, tree2)")
        from dendropy.calculate import treecompare
        return treecompare.find_missing_splits(self, other_tree)

    def symmetric_difference(self, other_tree):
        """DEPRECATED: Use 'dendropy.treecompare.symmetric_difference()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Statistics comparing two trees are now implemented in the 'dendropy.calculate.treecompare' module.",
                old_construct="tree1.symmetric_difference(tree2)",
                new_construct="from dendropy.calculate import treecompare\ntreecompare.symmetric_difference(tree1, tree2)")
        from dendropy.calculate import treecompare
        return treecompare.symmetric_difference(self, other_tree)

    def false_positives_and_negatives(self, other_tree):
        """DEPRECATED: Use 'dendropy.treecompare.false_positives_and_negatives()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Statistics comparing two trees are now implemented in the 'dendropy.calculate.treecompare' module.",
                old_construct="tree1.false_positives_and_negatives(tree2)",
                new_construct="from dendropy.calculate import treecompare\ntreecompare.false_positives_and_negatives(tree1, tree2)")
        from dendropy.calculate import treecompare
        return treecompare.false_positives_and_negatives(self, other_tree)

    def robinson_foulds_distance(self, other_tree):
        """DEPRECATED: Use 'dendropy.treecompare.weighted_robinson_foulds_distance()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Statistics comparing two trees are now implemented in the 'dendropy.calculate.treecompare' module, and this method's functionality is available through the 'weighted_robinson_foulds_distance()' function. For the *unweighted* RF distance, see 'dendropy.calculate.treecompare.symmetric_difference()'.",
                old_construct="tree1.robinson_foulds_distance(tree2)",
                new_construct="from dendropy.calculate import treecompare\ntreecompare.weighted_robinson_foulds_distance(tree1, tree2)")
        from dendropy.calculate import treecompare
        return treecompare.weighted_robinson_foulds_distance(self, other_tree)

    def euclidean_distance(self, other_tree):
        """DEPRECATED: Use 'dendropy.treecompare.euclidean_distance()'."""
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: Statistics comparing two trees are now implemented in the 'dendropy.calculate.treecompare' module.",
                old_construct="tree1.euclidean_distance(tree2)",
                new_construct="from dendropy.calculate import treecompare\ntreecompare.euclidean_distance(tree1, tree2)")
        from dendropy.calculate import treecompare
        return treecompare.euclidean_distance(self, other_tree)

    ###########################################################################
    ### Metadata

    def strip_comments(self):
        """
        Remove comments from tree/nodes.
        """
        self.comments = []
        for nd in self.postorder_node_iter():
            nd.comments = []
            nd.edge.comments = []

    ###########################################################################
    ### Representation

    def __str__(self):
        "Dump Newick string."
        return "%s" % self._as_newick_string()

    def __repr__(self):
        return "<{} object at {}>".format(self.__class__.__name__, hex(id(self)))

    def description(self, depth=1, indent=0, itemize="", output=None):
        """
        Returns description of object, up to level ``depth``.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s)" % id(self)
        else:
            label = " (%s: '%s')" % (id(self), self.label)
        output_strio.write('%s%sTree object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            newick_str = self._as_newick_string()
            if not newick_str:
                newick_str = "()"
            if depth == 1:
                output_strio.write(': %s' % newick_str)
            elif depth >= 2:
                num_nodes = len([nd for nd in self.preorder_node_iter()])
                num_edges = len([ed for ed in self.preorder_edge_iter()])
                output_strio.write(': %d Nodes, %d Edges' % (num_nodes, num_edges))
                if self.taxon_namespace is not None:
                    output_strio.write("\n%s[Taxon Set]\n" % (" " * (indent+4)))
                    self.taxon_namespace.description(depth=depth-1, indent=indent+8, itemize="", output=output_strio)
                output_strio.write('\n%s[Tree]' % (" " * (indent+4)))
                output_strio.write('\n%s%s' % (" " * (indent+8), newick_str))
                if depth >= 3:
                    output_strio.write("\n%s[Nodes]" % (" " * (indent+4)))
                    for i, nd in enumerate(self.preorder_node_iter()):
                        output_strio.write('\n')
                        nd.description(depth=depth-3, indent=indent+8, itemize="[%d] " % i, output=output_strio, taxon_namespace=self.taxon_namespace)
                    output_strio.write("\n%s[Edges]" % (" " * (indent+4)))
                    for i, ed in enumerate(self.preorder_edge_iter()):
                        output_strio.write('\n')
                        ed.description(depth=depth-3, indent=indent+8, itemize="[%d] " % i, output=output_strio, taxon_namespace=self.taxon_namespace)

        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    def as_python_source(self, tree_obj_name=None, tree_args=None, oids=False):
        """
        Returns string that will rebuild this tree in Python.
        """
        p = []

        if tree_obj_name is None:
            tree_obj_name = "tree_%s" % id(self)
        if self.label is not None:
            label = "'" + self.label + "'"
        else:
            label = "None"
        if tree_args is None:
            tree_args = ""
        else:
            tree_args = ", " + tree_args
        p.append("%s = dendropy.Tree(label=%s%s%s)" \
            % (tree_obj_name,
               label,
               oid_str,
               tree_args))

        taxon_obj_namer = lambda x: "tax_%s" % id(x)
        for taxon in self.taxon_namespace:
            tobj_name = taxon_obj_namer(taxon)
            if taxon.label is not None:
                label = "'" + taxon.label + "'"
            else:
                label = "None"
            p.append("%s = %s.taxon_namespace.require_taxon(label=%s%s)" \
                % (tobj_name,
                   tree_obj_name,
                   label,
                   oid_str))

        node_obj_namer = lambda x: "nd_%s" % id(x)
        for node in self.preorder_node_iter():
            for child in node.child_nodes():
                if node is self.seed_node:
                    nn = "%s.seed_node" % tree_obj_name
                else:
                    nn = node_obj_namer(node)
                if child.label is not None:
                    label = "'" + child.label + "'"
                else:
                    label = "None"
                if child.taxon is not None:
                    ct = taxon_obj_namer(child.taxon)
                else:
                    ct = "None"
                p.append("%s = %s.new_child(label=%s, taxon=%s, edge_length=%s%s)" %
                        (node_obj_namer(child),
                         nn,
                         label,
                         ct,
                         child.edge.length,
                         oid_str))

        return "\n".join(p)

    ###########################################################################
    ### Representation

    def as_ascii_plot(self, **kwargs):
        """
        Returns a string representation a graphic of this tree using ASCII
        characters.
        """
        ap = AsciiTreePlot(**kwargs)
        return ap.compose(self)

    def write_ascii_plot(self, stream, **kwargs):
        """
        Writes an ASCII text graphic of this tree to ``stream``.
        """
        return stream.write(self.as_ascii_plot(**kwargs))

    def print_plot(self, **kwargs):
        """
        Writes an ASCII text graphic of this tree to standard output.
        """
        import sys
        self.write_ascii_plot(sys.stdout, **kwargs)
        sys.stdout.write("\n")

    def write_as_dot(self, out, **kwargs):
        """
        Writes the tree to ``out`` as a DOT formatted digraph
        """
        if not kwargs.get("taxon_namespace"):
            kwargs["taxon_namespace"] = self.taxon_namespace
        out.write("digraph G {\n")

        nd_id_to_dot_nd = {}
        for n, nd in enumerate(self.preorder_node_iter()):
            label = _format_node(nd, **kwargs)
            if nd is self.seed_node:
                label = "root %s" % label
            dot_nd = "n%d" % n
            out.write(' %s  [label="%s"];\n' % (dot_nd, label))
            nd_id_to_dot_nd[nd] = dot_nd
        for nd, dot_nd in nd_id_to_dot_nd.iteritems():
            try:
                e = nd.edge
                par_dot_nd = nd_id_to_dot_nd[e.tail_node]
            except:
                pass
            else:
                label = _format_edge(e, **kwargs)
                s = ' %s -> %s [label="%s"];\n' % (par_dot_nd, dot_nd, label)
                out.write(s)
        out.write("}\n")

    ###########################################################################
    ### Debugging/Testing

    def _assign_node_labels_from_taxon(self):
        for nd in self.postorder_node_iter():
            if nd.label is not None:
                continue
            if nd.taxon is not None:
                nd.label = nd.taxon.label

    def _get_indented_form(self, **kwargs):
        out = StringIO()
        self._write_indented_form(out, **kwargs)
        return out.getvalue()

    def _write_indented_form(self, out, **kwargs):
        if kwargs.get("bipartitions"):
            if not kwargs.get("taxon_namespace"):
                kwargs["taxon_namespace"] = self.taxon_namespace
        self.seed_node._write_indented_form(out, **kwargs)

    def _debug_check_tree(self, logger_obj=None, **kwargs):
        import logging, inspect
        if logger_obj and logger_obj.isEnabledFor(logging.DEBUG):
            try:
                assert self._debug_tree_is_valid(logger_obj=logger_obj, **kwargs)
            except:
                calling_frame = inspect.currentframe().f_back
                co = calling_frame.f_code
                emsg = "\nCalled from file %s, line %d, in %s" % (co.co_filename, calling_frame.f_lineno, co.co_name)
                _LOG.debug("%s" % str(self))
                _LOG.debug("%s" % self._get_indented_form(**kwargs))
        assert self._debug_tree_is_valid(logger_obj=logger_obj, **kwargs)

    def _debug_tree_is_valid(self, **kwargs):
        """Performs sanity-checks of the tree data structure.

        kwargs:
            ``check_bipartitions`` if True specifies that the bipartition attributes are checked.
        """
        check_bipartitions = kwargs.get('check_bipartitions', False)
        unique_bipartition_edge_mapping = kwargs.get('unique_bipartition_edge_mapping', False)
        taxon_namespace = kwargs.get('taxon_namespace')
        if taxon_namespace is None:
            taxon_namespace = self.taxon_namespace
        if check_bipartitions:
            taxa_mask = self.seed_node.edge.bipartition._leafset_bitmask
        nodes = {}
        edges = {}
        curr_node = self.seed_node
        assert curr_node._parent_node is None, \
                "{} is seed node, but has non-'None' parent node: {}".format(curr_node, curr_node._parent_node)
        assert curr_node.edge.tail_node is None, \
                "{} is seed node, but edge has non-'None' tail node: {}".format(curr_node, curr_node.edge._parent_node)
        ancestors = []
        siblings = []
        while curr_node:
            assert curr_node not in nodes, \
                    "Node {} seen multiple times".format(curr_node)
            curr_edge = curr_node.edge
            assert curr_edge not in edges, \
                    "Edge of {}, {}, is also an edge of {}".format(curr_node, curr_node.edge, edges[curr_edge])
            edges[curr_edge] = curr_node
            nodes[curr_node] = curr_edge
            assert curr_edge.head_node is curr_node, \
                    "Head node of edge of {}, {}, is {}, not {}".format(curr_node, curr_edge, curr_edge.head_node, curr_node)
            assert curr_edge.tail_node is curr_node._parent_node, \
                    "Tail node of edge of {}, {}, is {}, but parent node is {}".format(curr_node, curr_edge, curr_edge.tail_node, curr_node._parent_node)
            if check_bipartitions:
                cm = 0
                assert (curr_edge.bipartition._leafset_bitmask | taxa_mask) == taxa_mask, \
                        "Bipartition mask error: {} | {} == {} (expecting: {})".format(
                                curr_edge.bipartition.leafset_as_bitstring(),
                                self.seed_node.edge.bipartition.leafset_as_bitstring(),
                                self.seed_node.edge.bipartition.bitmask_as_bitstring(curr_edge.bipartition._leafset_bitmask | taxa_mask),
                                self.seed_node.edge.bipartition.leafset_as_bitstring(), )
            c = curr_node._child_nodes
            if c:
                for child in c:
                    assert child._parent_node is curr_node, \
                            "Child of {}, {}, has {} as parent".format(curr_node, child, child._parent_node)
                    if check_bipartitions:
                        cm |= child.edge.bipartition._leafset_bitmask
            elif check_bipartitions:
                assert curr_node.taxon is not None, \
                        "Cannot check bipartitions: {} is a leaf node, but its 'taxon' attribute is 'None'".format(curr_node)
                cm = taxon_namespace.taxon_bitmask(curr_node.taxon)
            if check_bipartitions:
                assert (cm & taxa_mask) == curr_edge.bipartition._leafset_bitmask, \
                        "Bipartition leafset bitmask error: {} (taxa: {}, leafset: {})".format(
                                curr_edge.bipartition.bitmask_as_bitstring(cm),
                                curr_edge.bipartition.bitmask_as_bitstring(taxa_mask),
                                curr_edge.bipartition.leafset_as_bitstring())
                if unique_bipartition_edge_mapping:
                    assert self.bipartition_edge_map[curr_edge.bipartition] is curr_edge, \
                            "Expecting edge {} for bipartition {}, but instead found {}".format(curr_edge, curr_edge.bipartition, self.bipartition_edge_map[curr_edge.bipartition])
            curr_node, level = _preorder_list_manip(curr_node, siblings, ancestors)
        if check_bipartitions:
            for b in self.bipartition_encoding:
                e = self.bipartition_edge_map[b]
                assert e in edges, "{}: {} => {}".format(e, e.tail_node, e.head_node)
                if unique_bipartition_edge_mapping:
                    assert b is e.bipartition
        return True

    def _as_newick_string(self, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.
        For production purposes, use the the full-fledged 'as_string()'
        method of the object.
        """
        return self.seed_node._as_newick_string(**kwargs)

    def _print_newick(self, **kwargs):
        """
        Convenience method to newick string representation of this tree
        to the standard output stream.
        """
        import sys
        sys.stdout.write(self._as_newick_string(**kwargs))
        sys.stdout.write("\n")

    def _write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.  For
        production purposes, use the the full-fledged 'write_to_stream()'
        method of the object.
        """
        self.seed_node._write_newick(out, **kwargs)

    def _plot_bipartitions_on_tree(self,
            show_splits=True,
            show_leafsets=True,
            show_taxon_labels=False,
            is_bipartitions_updated=False,
            width=120):
        if not is_bipartitions_updated:
            self.encode_bipartitions()
        def _print_node(nd):
            d = []
            if show_splits:
                d.append(nd.bipartition.split_as_bitstring())
            if show_leafsets:
                d.append(nd.bipartition.leafset_as_bitstring())
            s = "/".join(d)
            if show_taxon_labels and nd.taxon is not None:
                s = s + " ({})".format(nd.taxon.label)
            return s
        return self.as_ascii_plot(
                show_internal_node_labels=True,
                node_label_compose_fn=_print_node,
                width=width,
                )

###############################################################################
### AsciiTreePlot

class AsciiTreePlot(object):

    class NullEdgeLengthError(ValueError):
        def __init__(self, *args, **kwargs):
            ValueError.__init__(self, *args, **kwargs)

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------

        plot_metric : str
            A string which specifies how branches should be scaled, one of:
            'age' (distance from tips), 'depth' (distance from root),
            'level' (number of branches from root) or 'length' (edge
            length/weights).
        show_internal_node_labels : bool
            Whether or not to write out internal node labels.
        leaf_spacing_factor : int
            Positive integer: number of rows between each leaf.
        width : int
            Force a particular display width, in terms of number of columns.
        node_label_compose_fn : function object
            A function that takes a Node object as an argument and returns
            the string to be used to display it.

        """
        self.plot_metric = kwargs.pop('plot_metric', 'depth')
        self.show_internal_node_labels = kwargs.pop('show_internal_node_labels', False)
        self.show_external_node_labels = kwargs.pop('show_internal_node_labels', True)
        self.leaf_spacing_factor = kwargs.pop('leaf_spacing_factor', 2)
#        self.null_edge_length = kwargs.pop('null_edge_length', 0)
        self.width = kwargs.pop('width', None)
        self.display_width = kwargs.pop('display_width', self.width) # legacy
        self.compose_node = kwargs.pop("node_label_compose_fn", None)
        if self.compose_node is None:
            self.compose_node = self.default_compose_node
        if kwargs:
            raise TypeError("Unrecognized or unsupported arguments: {}".format(kwargs))

    def default_compose_node(self, node):
        if node.taxon is not None and node.taxon.label is not None:
            return node.taxon.label
        elif node.label is not None:
            return node.label
        else:
            return "@"

    def reset(self):
        self.grid = []
        self.node_row = {}
        self.node_col = {}
        self.node_offset = {}
        self.current_leaf_row = 0
        self.node_label_map = {}

    def _calc_node_offsets(self, tree):
        if self.plot_metric == 'age' or self.plot_metric == 'depth':

            for nd in tree.postorder_node_iter():
                cnds = nd.child_nodes()
                if self.plot_metric == 'depth': # 'number of branchings from tip'
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        depths = [self.node_offset[v] for v in cnds]
                        curr_node_offset = max(depths) + 1
                elif self.plot_metric == 'age': # 'sum of edge weights from tip'
                    # note: no enforcement of ultrametricity!
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        if cnds[0].edge.length is not None:
                            curr_node_offset = self.node_offset[cnds[0]] + cnds[0].edge.length
                else:
                    raise ValueError("Unrecognized plot metric '%s' (must be one of: 'age', 'depth', 'level', or 'length')" % self.plot_metric)
                self.node_offset[nd] = curr_node_offset
            flipped_origin = max(self.node_offset.values())
            for nd in self.node_offset:
                self.node_offset[nd] = flipped_origin - self.node_offset[nd]
        else:
            for nd in tree.preorder_node_iter():
                if self.plot_metric == 'level': # 'number of branchings from root'
                    curr_edge_len = 1
                elif self.plot_metric == 'length': # 'sum of edge weights from root'
                    if nd.edge.length is not None:
                        curr_edge_len = nd.edge.length
                    else:
                        curr_edge_len = 0
                else:
                    raise ValueError("Unrecognized plot metric '%s' (must be one of: 'age', 'depth', 'level', or 'length')" % self.plot_metric)
                if nd._parent_node is None:
                    self.node_offset[nd] = curr_edge_len
                else:
                    self.node_offset[nd] =  curr_edge_len + self.node_offset[nd._parent_node]
#        print "\n".join([str(k) for k in self.node_offset.values()])

    def draw(self, tree, dest):
        dest.write(self.compose(tree))

    def get_label_for_node(self, node):
        try:
            return self.node_label_map[node]
        except KeyError:
            if node._child_nodes and self.show_internal_node_labels:
                label = self.compose_node(node)
            elif not node._child_nodes and self.show_external_node_labels:
                label = self.compose_node(node)
            else:
                label = ""
            self.node_label_map[node] = label
            return label

    def compose(self, tree):
        self.reset()
        if self.display_width is None:
            display_width = terminal.terminal_width() - 1
        else:
            display_width = self.display_width
        max_label_len = max([len(self.get_label_for_node(i)) for i in tree.leaf_node_iter()])
        if max_label_len <= 0:
            max_label_len = 0
        #effective_display_width = display_width - max_label_len - len(tree.internal_nodes) - 1
        effective_display_width = display_width - max_label_len - 1
        self._calc_node_offsets(tree)
        widths = [self.node_offset[i] for i in tree.leaf_node_iter() if self.node_offset[i] is not None]
        max_width = float(max(widths))
        if max_width == 0:
            raise AsciiTreePlot.NullEdgeLengthError("Tree cannot be plotted under metric '%s' due to zero or null edge lengths: '%s'" % (self.plot_metric, tree._as_newick_string()))
        edge_scale_factor = float(effective_display_width) / max_width
        self.calc_plot(tree.seed_node,
                       edge_scale_factor=edge_scale_factor)
        for i in range(len(tree.leaf_nodes())*self.leaf_spacing_factor + 1):
            self.grid.append([' ' for i in range(0, display_width)])
        self.draw_node(tree.seed_node)
        display = '\n'.join([''.join(i) for i in self.grid])
        return display

    def calc_plot(self, node, edge_scale_factor):
        """
        First pass through tree, post-order traversal to calculate
        coordinates of each node.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            for n in child_nodes:
                self.calc_plot(n, edge_scale_factor)
            ys = [self.node_row[n] for n in child_nodes]
            self.node_row[node] = int(float((max(ys)-min(ys)) / 2) + min(ys))
        else:
            self.node_row[node] = self.current_leaf_row
            self.current_leaf_row = self.current_leaf_row + self.leaf_spacing_factor
        if node.edge.length is None:
            self.node_col[node] = 1
        else:
            self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)
        self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)

    def draw_label(self, label, row, start_col):
        if label:
            for i in range(len(label)):
                if start_col + i < len(self.grid[row]):
                    self.grid[row][start_col+i] = label[i]

    def draw_node(self, node):
        """
        Second pass through tree, plotting nodes onto given self.grid.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            for i, child_node in enumerate(child_nodes):
                start_row = min([self.node_row[node], self.node_row[child_node]])
                end_row = max([self.node_row[node], self.node_row[child_node]])
                if i == 0:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '/'
                    start_row = start_row+1
                    edge_row = self.node_row[child_node]
                elif i == len(child_nodes)-1:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '\\'
                    edge_row = self.node_row[child_node]
                else:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '+'
                    edge_row = self.node_row[child_node]
                self.draw_node(child_node)
                for x in range(self.node_col[node]+1, self.node_col[child_node]):
                    self.grid[edge_row][x] = '-'
                for y in range(start_row, end_row):
                    self.grid[y][self.node_col[node]] = '|'
            label = []
            if self.show_internal_node_labels:
                label = self.get_label_for_node(node)
                self.draw_internal_text(label, self.node_row[node], self.node_col[node])
            else:
                self.grid[self.node_row[node]][self.node_col[node]]='+'
        else:
            label = self.get_label_for_node(node)
            self.draw_label(label, self.node_row[node], self.node_col[node]+1)

    def draw_internal_text(self, label, r, c):
        row = self.grid[r]
        try:
            for n, letter in enumerate(label):
                row[c + n] = letter
        except:
            pass

###############################################################################
### Helper Functions

def _preorder_list_manip(n, siblings, ancestors):
    """
    Helper function for recursion free preorder traversal, that does
    not rely on attributes of the node other than child_nodes() (thus it
    is useful for debuggging).

    Returns the next node (or None) and the number of levels toward the
    root the function "moved".
    """
    levels_moved = 0
    c = n.child_nodes()
    if c:
        levels_moved += 1
        ancestors.append(list(siblings))
        del siblings[:]
        siblings.extend(c[1:])
        return c[0], levels_moved
    while not siblings:
        if ancestors:
            levels_moved -= 1
            del siblings[:]
            siblings.extend(ancestors.pop())
        else:
            return None, levels_moved
    return siblings.pop(0), levels_moved

def _format_node(nd, **kwargs):
    nf = kwargs.get('node_formatter', None)
    if nf:
        return nf(nd)
    if nd.taxon is not None:
        return str(nd.taxon)
    if nd.label is not None:
        return nd.label
    return ""

def _format_edge(e, **kwargs):
    ef = kwargs.get('edge_formatter', None)
    if ef:
        return ef(e)
    return str(e)

def _format_split(split, length=None, **kwargs):
    if length is None:
        length = len(kwargs.get("taxon_namespace"))
    return bitprocessing.int_as_bitstring(split, length=length)

def _convert_node_to_root_polytomy(nd):
    """If ``nd`` has two children and at least on of them is an internal node,
    then it will be converted to an out-degree three node (with the edge length
    added as needed).

    Returns a tuple of child nodes that were detached (or() if the tree was not
    modified). This can be useful for removing the deleted node from the split_edge_map
    dictionary.
    """
    nd_children = nd.child_nodes()
    if len(nd_children) > 2:
        return ()
    try:
        left_child = nd_children[0]
    except:
        return ()
    if not left_child:
        return ()
    if len(nd_children) == 1:
        right_child = None
        dest_edge_head = nd
    else:
        right_child = nd_children[1]
        dest_edge_head = right_child
    curr_add = None
    if right_child and right_child.is_internal():
        try:
            left_child.edge.length += right_child.edge.length
        except:
            pass
        nd.remove_child(right_child)
        grand_kids = right_child.child_nodes()
        for gc in grand_kids:
            nd.add_child(gc)
        curr_add = right_child
    elif left_child.is_internal():
        try:
            dest_edge_head.edge.length += left_child.edge.length
        except:
            pass
        nd.remove_child(left_child)
        grand_kids = left_child.child_nodes()
        for gc in grand_kids:
            nd.add_child(gc)
        curr_add = left_child
    if curr_add:
        ndl = [curr_add]
        t = _convert_node_to_root_polytomy(nd)
        ndl.extend(t)
        return tuple(ndl)
    return ()

