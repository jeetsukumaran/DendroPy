************
Bipartitions
************

Many tree statistics and operations in DendroPy use the *bipartition encoding*
of a :class:`Tree` instance in the background, including, for example:

    -   tree statistics and metrics
    -   tree comparisons
    -   tree scoring

By default, the DendroPy functions assume that bipartitions are *not* encoded,
or are not up-to-date with respect to the current tree structure, resulting in
their recalculation *every* time. This is computationally inefficient, and you
want to avoid it if, indeed, the bipartition encoding of a tree is current. You
can control whether or not these service functions recalculate the bipartition
encoding by passing in the argument ``is_bipartitions_updated=True`` to
suppress the recalculation or ``is_bipartitions_updated=False`` to force it.

If you are doing multiple operations that require a bipartition encoding, you
should call :class:`Tree.encode_bipartitions()` *once* for each tree, and,
then, as long as the trees are *not* *modified* since the encoding, specify the
``is_bipartitions_updated=True`` argument to each of the functions that use it
to ensure that the bipartitions are not recalculated each time.

If, on the other hand, you modify a tree structure in any way, e.g., rerooting,
pruning, add/removing nodes or subtrees, you should update the bipartition
encoding of a tree yourself by calling :class:`Tree.encode_bipartitions()`, or
make sure to specify ``is_bipartitions_updated=False`` to the *first* function
that you call following the tree modification.

Modeling Bipartitions
=====================

A Bipartition is a Partitioning of Taxa Corresponding to an Edge of a Tree
--------------------------------------------------------------------------

A bipartition is the division or sorting of the leaves/tips of a tree into two
mutually-exclusive and collectively-exhaustive subsets (i.e., a *partition*, in
the set theory sense, of the leaves of the tree into exactly two non-empty
subsets; hence the term, "*bi*-partition"). Every edge on a tree corresponds to
a bipartition in the sense that if were were to split or bisect a tree at a
particular edge, the leaf sets of each of the two new trees constitute the a
bipartition of the leaf set of the original tree. In the context of
evolutionary trees like a phylogeny, the leaves typically are associated with
operational taxonomic unit concepts, or, for short, taxa. So, just as we view a
tree as a schematic representation of the relationships of taxa, we can see
bipartitions as a representation of a clustering of taxa.

For example, given a tree:

    ((a,(b,c)),(d,(e,f)));

the edge subtending the leaf node with taxon "d" corresponds to the bipartition
that splits "d" from the rest of the taxa. Similarly, the edge subtending the
most-recent common ancestor (MRCA) node of taxa "d", "e", and "f" corresponds
to the bipartition that splits "d", "e", and "f" from the rest of the taxa,
"a", "b", and "c".

A Bipartition Can Be Described by a *Bitmask*
---------------------------------------------

If we were to index the taxa of the tree, with the first taxon getting index 1,
the second index 2, the third index 3, etc. and so on until index $n$, we can
represent any possible split as sequence of symbols, such as:

    abbabbaa

where the symbol indicates membership in one arbitrarily-labeled group (e.g.,
"a") or the other (e.g., "b") of a particular taxon, based on how we relate the
taxon indexes to the position of the symbols in sequence.

If we were to use a left-to-right order, such that the first element
corresponded to the first taxon, the second to the second taxon, and so one,
the above sequence would describe the a partition of the taxa {1,2,...,8} into
the sets {1,4,7,8} and {2,3,5,6}. However, in DendroPy, we use a right-to-left
order (for reasons explained below), such that the right-most element
corresponds to the taxon with index 1, the next right-most element corresponds
to the taxon with index 2 and so on, so the sequence above represents a
partition of the taxa {1,2,...,8} into the sets {1,2,5,8} and {3,4,6,7}.

Let us say that we had a set of 8 taxa {A,B,...,H}:

    A, B, C, D, E, F, G, H

which we assign indexes {1,2,...,8} according to the following scheme:

    +-------+---+---+---+---+---+---+---+---+
    | Taxon | A | B | C | D | E | F | G | H |
    +-------+---+---+---+---+---+---+---+---+
    | Index | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
    +-------+---+---+---+---+---+---+---+---+

Then we can describe a bipartition that divides the taxa into two groups
{A,B,E,H} and {C,D,F,G}, using right-to-left ordering and symbols "0" and "1"
(instead of "a" and "b") as:

    +-------------+---+---+---+---+---+---+---+---+
    | Taxon Index | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 |
    +-------------+---+---+---+---+---+---+---+---+
    | Group       | 0 | 1 | 1 | 0 | 1 | 1 | 0 | 0 |
    +-------------+---+---+---+---+---+---+---+---+

We can succintly and usefully represent the bipartition above with an integer
given by interpreting the sequence of 0's and 1's as bits. Interpreting the
sequence above, "01101100", as a binary number or bitmask means that this
bipartition can be represented as the decimal integer "108".

This, in essence, is how bipartitions are represented in DendroPy: as integers
that are interpreted as *bitmasks* (also known as bit arrays, bit vectors, or
bit fields, though exact application of terminology varies depending on whether
or not primary operations are bitwise or dereferencing through offset indexes
or field names, etc.), where the 0's and 1's assign taxa to different subsets
of the bipartition.

As an example, consider the following tree::

    [&R] (A,(B,(C,(D,E))));

This would be encoded as::

    /----------------------------- 00001 (A)
    11111
    |      /---------------------- 00010 (B)
    \------11110
        |       /-------------- 00100 (C)
        \-------11100
                |      /------- 01000 (D)
                \------11000
                        \------- 10000 (E)

The leaves are assigned bitmasks based on the indexes of the taxa, while the
internal nodes are given by a `bitwise OR <http://en.wikipedia.org/wiki/Bitwise_operation#OR>`_-ing of the bitmasks of their children.

Modeling Bipartitions Using Leafset Bitmasks and Split Bitmasks
---------------------------------------------------------------

In DendroPy, bipartitions are modeled using bitmasks as discussed above, i.e.,
integers that, when represented as a bitarray or bitstring, specify the
assignment of taxa into one of two groups, based on whether or not the bit
corresponding to the taxon index is set or not.

In fact, each bipartition is actually modeled by *two* types of bitmasks: a
*leafset bitmask* and a *split bitmask*:

    - A leafset bitmask is a bit array in which the presence of a taxon in the
      leaves descending from the edge associated with the bipartition is
      represented by a set bit ("1"), while its absence is represented by an
      unset bit ("0"). The taxa are mapped to bit positions using a
      least-significant bit mapping scheme, in which the first taxon is represented by
      the least significant bit, the second taxon is represented by the next
      most significant bit, and so on.

    - A split bitmask is a bit array which divides or partitions taxa by assign
      each taxon to one of two arbitrarily-labeled groups, "0" or "1",
      depending on whether or not a bit is set or not in the position
      corresponding the taxon index under a least-signficant bit mapping scheme
      as described above.

        - For bipartitions of rooted trees, the split bitmask is the same value
          as the leafset bitmask.

        - For bipartitions of unrooted trees, the split bitmask is the same
          value as the leafset bitmask *if and only if* the least-signficant
          bit of the leafset bitmask is 0 (i.e., the first taxon is assigned to
          group "0"), or the *complement* of the leafset bitmask if this is the
          case. In other words, with unrooted trees we constrain the split
          bitmasks such that the first taxon and all other taxa grouped
          together with it are always placed in group "0".

Why this complication?

Consider the following unrooted tree::

    A    C    D
     \   |   /
      +--+--+
     /       \
    B         E

This could be represented by either of the following NEWICK strings::

    [&U] ((A,B),(C,(D,E)));
    [&U] (((A,B),C),(D,E));

Both the above topologies, while distinct if interpreted as rooted, represent
*identical* unrooted toplogies.

When the bipartitions are encoded as leafset bitmasks, we get the following if
the first tree statement is parsed by DendroPy::

                        /--------- 00001 (A)
    /-------------------00011
    |                   \--------- 00010 (B)
    11111
    |         /------------------- 00100 (C)
    \---------11100
            |         /--------- 01000 (D)
            \---------11000
                        \--------- 10000 (E)

and the following if the second tree statement is parsed by DendroPy::

                        /--------- 00001 (A)
            /---------00011
    /---------00111     \--------- 00010 (B)
    |         |
    11111     \------------------- 00100 (C)
    |
    |                   /--------- 01000 (D)
    \-------------------11000
                        \--------- 10000 (E)

Note that the leafset bitmask "11100" in the first tree is absent in the second
tree, while conversely, the leafset bitmask "00111" in the second tree is
absent in the first tree.

This difference is due purely to the placement of the root to one side or the
other of taxon 'C'. In rooted trees, this root is a real root, and this
difference in bipartitions as given by the leafset bitmasks is also real. In
unrooted trees, this "root" is actually an artifact of the tree structure, and
the placement is an artifact of the NEWICK string representation. In unrooted
trees, then the difference in bipartitions as given by the leafset bitmasks is,
thus, wholly artifactual. This means that it would be impossible to robustly
and reliably compare, relate, and perform any operations on bipartitions coded
using leafset bitmasks on unrooted trees: what is effectively the equivalent
bipartition of taxa maybe represented either by placing, the first taxon and all
the other taxa in the same group as it in group "0" in one representation, or
group "1" in another one, and which representation is used is arbitrary and
random and unpredictable.

Thus, to allow for robustly establishing equivalence of bipartitions across
different representations and instantiations of different unrooted trees, we
*normalize* the bit array representation of bipartitions in unrooted trees to
always ensure that the first taxon is assigned to group "0", *whether* *or*
*not* *this* *taxon* *is* *actually* *a* *descendent* *or* *a* *member* *of*
*the* *leafset* *of* *the* *edge*. [We also collapse the basal bifurcation of
unrooted trees to avoid redundant representation of artifactual bipartitions.]

As the first taxon corresponds to the least-significant bit in the DendroPy
scheme, this normalization is known as the least-significant bit 0 or "LSB-0"
normalization scheme. The choice of 0 as opposed to 1 is arbitrary, but the
reason is so ensure that we can have consistent comparisons of groups across
trees of different rotations (and "pseudo-rootings" created by the constraints
of tree representation in, e.g., the NEWICK format) by enforcing the convention
that group "0" will always be the group that includes the first taxon (i.e.,
the taxon with index 1, corresponding to the position of the least-significant
or right-most bit).

We refer to this normalized version of the leafset bitmask as a *split
bitmask*. For consistency, bipartitions of rooted trees are also assigned split
bitmasks, but here these are simply the unmodifed leafset bitmasks. For both
unrooted and rooted trees we maintain the leafset bitmask representation in
parallel for each bipartition, as this has useful information is lost when
normalized, e.g., establishing whether or not a particular subtree or taxon can
be found within bipartition.

Thus, regardless of whether the tree is rooted or unrooted, each bipartition on
is modeled by *two* bitmasks: a split bitmask and leafset bitmask. For rooted
trees, these are identical in value. For unrooted trees, the split bitmask is
the leafset bitmask normalized to constrain the least significant bit to be 0.

A split bitmask is used to establish *identity* across different trees (for
this reason it is also sometimes called a split or bipartition hash), while a
leafset bitmask is used to work with various ancestor-descendent relationships
within the same tree (it can be used to, for example, quickly assess if a taxon
descends from a particular node, or if a particular node is a common ancestor
of two taxa).

Leafset bitmasks are unstable representations of bipartitions for unrooted
trees, but remain accurate and convenient representations of the descendent
leaf-sets of nodes in both unrooted and rooted trees. Split bitmasks, on the
other hand, *are* stable representations of bipartitions for both unrooted as
well as rooted trees, but are not accurate representations of the taxa
associated with the leaves descended from the bipartition of a particular edge.

Using Bipartitions
==================

Bipartition Encoding
--------------------

The bipartition encoding of a tree is a specification of the structure of tree
in terms of the complete set of bipartitions that can be found on it. Given a
bipartition encoding of a tree, the entire topology can be reconstructed
completely and accurately. In addition, the bipartition encoding of trees can
be used to quickly and accurately compare, relate, and calculate various
statistics between different trees and within the same tree.

In DendroPy, the :meth:`Tree.encode_bipartitions()` method calculates the
bipartitions of a tree. The :attr:`Edge.bipartition` attribute of each edge
will be populated by a :class:`Bipartition` instance, each of which has the
bipartition's split bitmask stored in the :attr:`Bipartition.split_bitmask`
attribute and the leafset bitmask stored in the
:attr:`Bipartition.leaf_bitmask` attribute. In addition, each
:class:`Bipartition` also stores a reference to the edge to which it
corresponds in its :attr:`Bipartition.edge` attribute. For convenience, the
split bitmask and the leafset bitmask associated with each bipartition of an
edge can be also be accesed through the :attr:`Edge.split_bitmask` and
:attr:`Ede.leafset_bitmask` properties, respectively.

You can access these :class:`Bipartition` objects by iterating over the edges
of the tree, but it might be more convenient to access them through the
:attr:`Tree.bipartition_encoding` attribute of the :class:`Tree`. You can also
access a dictionary mapping :class:`Bipartition` instances to their
corresponding edges through the :attr:`Tree.bipartition_edge_map` attribute, or
a dictionary mapping split bitmasks to their corresponding edges through the
:attr:`Tree.split_bitmask_edge_map` attribute.

By default, the :class:`Bipartition` instances created are immutable. This is
to allow them to be used in sets or dictionary keys, and thus exploit O(1)
look-up/access performance. The hash value of a :class:`Bipartition` object is
its :attr:`Bipartition.split_bitmask` attribute; two distinct
:class:`Bipartition` objects are considered equivalent even if they refer to
different :class:`Edge` objects on different :class:`Tree` objects if their
:attr:`Bipartition.split_bitmask` values are the same. If you need to modify
the values of a :class:`Bipartition`, you need to set the
:attr:`Bipartition.is_mutable` attribute to ``True``. Note that changing any
values that modify the hash of a :class:`Bipartition` instance that is already
in a hash container such as a set or dictionary will make that instance or
possibly other members of the container inaccessible: never change the value of
a :class:`Bipartition` instance if it is in a set or dictionary.

Calculating Bipartitions on Trees
---------------------------------

A large number of DendroPy functions calculate the split and leafset bitmasks
in the background: from tree comparison approaches (e.g., calculating the
Robinson-Foulds distance), to working with within-tree operations (e.g.,
finding the most-recent common ancestor between two nodes or patrisitic
distances between taxa), to tree-set operations (e.g., building consensus trees
or scoring tree leafset credibilities and finding the maximum leafset
credibility tree).

When passing trees to these methods and functions, these functions will call
:meth:`Tree.encode_bipartitions()` automatically for you unless you explicitly
specify that this should not be done by passing in
'``is_bipartitions_updated=True``'.

The typical usage idiom in this context would be to:

    (1) Establish a common taxon namespace [i.e., creating a global
        :class:`TaxonNamespace` object and pass it in to all
        reading/parsing/input operations]
    (2) Read/load the trees, calling :meth:`Tree.encode_bipartitions()` on each one.
    (3) Perform the calculations, making sure to specify ``is_bipartitions_updated=True``.

For, example, the following snippet shows how you might count the number of trees in a bootstrap file that have the same topology as a tree of interest::

    import dendropy
    from dendropy.calculate import treecompare
    taxa = dendropy.TaxonNamespace()
    target_tree = dendropy.Tree.get_from_path(
        "mle.tre",
        "nexus",
        taxon_namespace=taxa)
    count = 0
    for sup_tree in dendropy.Tree.yield_from_files(
        files=["boots1.tre", "boots2.tre", "boostraps3.tre"],
        schema="nexus",
        taxon_namespace=taxa):
        d = treecompare.symmetric_difference(target_tree, sup_tree)
        if d == 0:
            count += 1
    print(count)

For this application, it is simpler just to let the calculations take place in the background. But, for example, if for some reason you wanted to do something more complicated, as it calculating the counts with respect to multiple trees of interest, you should try and avoid the redundant recalculation of the bitmasks::

    import dendropy

    from dendropy.calculate import treecompare
    taxa = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get_from_path(
        "mle1.tre",
        "nexus",
        taxon_namespace=taxa)
    tree1.encode_bipartitions()
    tree2 = dendropy.Tree.get_from_path(
        "mle2.tre",
        "nexus",
        taxon_namespace=taxa)
    tree2.encode_bipartitions()
    counts1 = 0
    counts2
    for sup_tree in dendropy.Tree.yield_from_files(
        files=["boots1.tre", "boots2.tre", "boostraps3.tre"],
        schema="nexus",
        taxon_namespace=taxa):
        sup_tree.encode_bipartitions()
        if treecompare.symmetric_difference(
                tree1, sup_tree, is_bipartitions_updated=True):
            count1 += 1
        if treecompare.symmetric_difference(
                tree2, sup_tree, is_bipartitions_updated=True):
            count2 += 2
    print(count1, count2)

Note that in all cases, for bipartitions to be meaningfully compared two conditions must hold:

    1. The trees must reference the *same* operational taxonomic unit namespace
       object, :class:`TaxonNamespace`.
    2. The trees must have the same rooting state (i.e., all rooted or all
       unrooted).

