#! /usr/bin/env python
# -*- coding: utf-8 -*-

from dendropy.utility import bitprocessing

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
    desirable property, and thus the adoption of the LSB scheme.

    Constraining the first taxon to be in group 0 (LSB-0) rather than group 1
    (LSB-1) is motivated by the fact that, in the former, we would combine
    the bitmasks of child nodes using OR (logical addition) operations when
    calculating the bitmask for a parent node, whereas, with the latter, we
    would need to use AND operations. The former strikes us as more intuitive.

    """

    @staticmethod
    def normalize_bitmask(bitmask, fill_bitmask, lowest_relevant_bit=1):
        if bitmask & lowest_relevant_bit:
            return (~bitmask) & fill_bitmask  # force least-significant bit to 0
        else:
            return bitmask & fill_bitmask  # keep least-significant bit as 0

    @staticmethod
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

    @staticmethod
    def is_trivial_leafset(leafset_bitmask):
        return bitprocessing.num_set_bits(leafset_bitmask) == 1

    @staticmethod
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

    ## Life-cycle

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------
        bitmask : integer
            A bit array representing the membership of taxa, with the
            least-significant bit corresponding to the first taxon, the next
            least-significant bit corresponding to the second taxon, and so on,
            till the last taxon corresponding to the most-significant bit.
            Taxon membership in one of two arbitrary groups, '0' or '1', is
            indicated by its corresponding bit being unset or set,
            respectively.
        leafset_bitmask : integer
            A bit array representing the presence or absence of taxa in the
            subtree descending from the child node of the edge of which this
            bipartition is associated. The least-significant bit corresponds to
            the first taxon, the next least-significant bit corresponds to the
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
        is_mutable : bool
            Specifies whether or not the tree is mutable.
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
                tree_leafset_bitmask=self._tree_leafset_bitmask,
            )
            if is_mutable is None:
                self.is_mutable = True
            else:
                self.is_mutable = is_mutable
        elif is_mutable is not None:
            self.is_mutable = is_mutable

    ## Identity

    def __hash__(self):
        assert not self.is_mutable, "Bipartition is mutable: hash is unstable"
        return self._split_bitmask or 0

    def __eq__(self, other):
        # return self._split_bitmask == other._split_bitmask
        return (
            self._split_bitmask is not None
            and self._split_bitmask == other._split_bitmask
        ) or (self._split_bitmask is other._split_bitmask)

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

    tree_leafset_bitmask = property(
        _get_tree_leafset_bitmask, _set_tree_leafset_bitmask
    )

    def _get_is_rooted(self):
        return self._is_rooted

    def _set_is_rooted(self, value):
        assert self.is_mutable, "Bipartition instance is not mutable"
        self._is_rooted = value

    is_rooted = property(_get_is_rooted, _set_is_rooted)

    ## Representation

    def __str__(self):
        return bin(self._split_bitmask)[2:].rjust(
            bitprocessing.bit_length(self._tree_leafset_bitmask), "0"
        )

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
            mask=self._split_bitmask, symbol0=symbol0, symbol1=symbol1, reverse=reverse
        )

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
            reverse=reverse,
        )

    def bitmask_as_bitstring(self, mask, symbol0=None, symbol1=None, reverse=False):
        return bitprocessing.int_as_bitstring(
            mask,
            length=bitprocessing.bit_length(self._tree_leafset_bitmask),
            symbol0=symbol0,
            symbol1=symbol1,
            reverse=reverse,
        )

    ## Calculation

    def compile_tree_leafset_bitmask(
        self, tree_leafset_bitmask, lowest_relevant_bit=None
    ):
        """
        Avoids recalculation of ``lowest_relevant_bit`` if specified.
        """
        assert self.is_mutable, "Bipartition instance is not mutable"
        self._tree_leafset_bitmask = tree_leafset_bitmask
        if lowest_relevant_bit is not None:
            self._lowest_relevant_bit = lowest_relevant_bit
        elif self._tree_leafset_bitmask:
            self._lowest_relevant_bit = bitprocessing.least_significant_set_bit(
                self._tree_leafset_bitmask
            )
        else:
            self._lowest_relevant_bit = None
        return self._tree_leafset_bitmask

    def compile_leafset_bitmask(self, leafset_bitmask=None, tree_leafset_bitmask=None):
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

    def compile_split_bitmask(
        self,
        leafset_bitmask=None,
        tree_leafset_bitmask=None,
        is_rooted=None,
        is_mutable=True,
    ):
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
                lowest_relevant_bit=self._lowest_relevant_bit,
            )
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
        self.compile_split_bitmask(
            self,
            leafset_bitmask=self._leafset_bitmask,
            tree_leafset_bitmask=self._tree_leafset_bitmask,
            is_rooted=self._is_rooted,
            is_mutable=is_mutable,
        )

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
        return (m1 & m2) == m1

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
        return (m2 & self._leafset_bitmask) == self._leafset_bitmask

    def is_trivial(self):
        """
        Returns
        -------
        bool
            |True| if this bipartition divides a leaf and the rest of the
            tree.
        """
        return Bipartition.is_trivial_bitmask(
            self._split_bitmask, self._tree_leafset_bitmask
        )

    def split_as_newick_string(
        self, taxon_namespace, preserve_spaces=False, quote_underscores=True
    ):
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
            quote_underscores=quote_underscores,
        )

    def leafset_as_newick_string(
        self, taxon_namespace, preserve_spaces=False, quote_underscores=True
    ):
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
            quote_underscores=quote_underscores,
        )

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
            bitmask=self._leafset_bitmask, index=index
        )

    def _format_bipartition(self, length=None, **kwargs):
        if length is None:
            length = len(kwargs.get("taxon_namespace"))
        return bitprocessing.int_as_bitstring(self, length=length)

    # def as_newick_string
    # def is_trivial
    # def is_non_singleton
    # def leafset_hash
    # def leafset_as_bitstring
    # def is_compatible
