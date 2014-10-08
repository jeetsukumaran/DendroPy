*******************************************
Splits, Split Bitmasks and Leafset Bitmasks
*******************************************

Splits or Bipartitions are a Division of Leaves into Two Groups
---------------------------------------------------------------

A "split" is a bipartition of a leaf or tip set of a tree, i.e., the division or sorting of the leaves/tips of a tree into two mutually-exclusive and collectively-comprehensive subsets.
The term "bipartition" is often also used to refer to the same concept. For example, the Mr. Bayes documentation says that: "taxon bipartition is defined by removing a branch on the tree, dividing the tree into those taxa to the left and right of the removed branch. This set is called a taxon bipartition"

Edges Correspond to Splits, and Splits Correspond to Partitions of Taxa
-----------------------------------------------------------------------

Every edge on a tree corresponds to a split in the sense that if were were to "split" or bisect a tree at a particular edge, the leaf sets of each of the two new trees constitute the a bipartition of the leaf set of the original tree.

In the context of evolutionary trees like a phylogeny, the leaves typically are associated with operational taxonomic unit concepts, or, for short, taxa. So, just as we view a tree as a schematic representation of the relationships of taxa, we can see splits as a representation of a (bi-)partition of taxa.

For example, given a tree:

    ((a,(b,c)),(d,(e,f)));

the edge subtending the leaf node with taxon "d" corresponds to the split that divides "d" from the rest of the taxa. Similarly, the edge subtending the most-recent common ancestor (MRCA) node of taxa "d", "e", and "f" corresponds to the split that separates "d", "e", and "f" from the rest of the taxa, "a", "b", and "c".

A Split Can Be Represented by a Sequence of Symbols Which In Turn Can be Represented By an Integer
--------------------------------------------------------------------------------------------------

If we were to index the taxa of the tree, with the first taxon getting index 1, the second index 2, the third index 3, etc. and so on until index $n$, we can represent any possible split as sequence of symbols, such as:

    abbabbaa

where the symbol indicates membership in one arbitrarily-labeled group (e.g., "a") or the other (e.g., "b") of a particular taxon, based on how we relate the taxon indexes to the position of the symbols in sequence.

If we were to use a left-to-right order, such that the first element corresponded to the first taxon, the second to the second taxon, and so one, the above sequence would describe the a partition of the taxa {1,2,...,8} into the sets {1,4,7,8} and {2,3,5,6}.

However, in DendroPy, we use a right-to-left order (for reasons explained below), such that the right-most element corresponds to the taxon with index 1, the next right-most element corresponds to the taxon with index 2 and so on, so the sequence above represents a partition of the taxa {1,2,...,8} into the sets {1,2,5,8} and {3,4,6,7}.

Splits Can be Represented as Integer Bitmasks
---------------------------------------------

Let us say that we had a set of 8 taxa {A,B,...,H}:

    A, B, C, D, E, F, G, H

which we assign indexes {1,2,...,8} according to the following scheme:

    +-------+---+---+---+---+---+---+---+---+
    | Taxon | A | B | C | D | E | F | G | H |
    +-------+---+---+---+---+---+---+---+---+
    | Index | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
    +-------+---+---+---+---+---+---+---+---+

Then we can describe a split that divides the taxa into two groups {A,B,E,H} and {C,D,F,G}, using right-to-left ordering and symbols "0" and "1" (instead of "a" and "b") as:

    +-------------+---+---+---+---+---+---+---+---+
    | Taxon Index | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 |
    +-------------+---+---+---+---+---+---+---+---+
    | Group       | 0 | 1 | 1 | 0 | 1 | 1 | 0 | 0 |
    +-------------+---+---+---+---+---+---+---+---+

We can succintly and usefully represent the split that induces the partition above with an integer given by interpreting the sequence of 0's and 1's as bits. Interpreting the sequence above, "01101100", as a binary number or bitfield means that this split can be represented as the decimal integer "108".

This, in essence, is how splits are represented in DendroPy: as integers that are interpreted as bitfields, or the terminology preferred by DendroPy, *bitmasks*, where the 0's and 1's assign taxa to different subsets of the bipartition.

Split Bitmasks and Clade Bitmasks
---------------------------------

In DendroPy, a "split bitmask" describes the split represented by that edge (i.e., the partition of the taxa induced by bisecting the tree on that edge).
A "clade bitmask", on the other hand, describes the taxa associated with the descendent leaves of that node, with a "0" indicating the absence of that taxon and a "1" indicating the presence of that taxon.
A split bitmask is used to establish *split* *identity* across different trees (for this reason it is also sometimes called a split hash), while a clade bitmask is used to work with various ancestor-descendent relationships within the same tree (it can be used to, for example, quickly assess if a taxon descends from a particular node, or if a particular node is a common ancestor of two taxa).
For rooted trees, the split bitmask and the clade bitmask of a particular edge are *always* the same value.
For unrooted trees, however, this is not necessarily the case.
This is because the split bitmasks of unrooted trees are *normalized* to always have a least-significant bit of 0. The choice of 0 as opposed to 1 is arbitrary, but the reason is so ensure that we can have consistent comparisons of groups across trees of different rotations (and "pseudo-rootings" created by the constraints of tree representation in, e.g., the NEWICK format) by enforcing the convention that group "0" will always be the group that includes the first taxon (i.e., the taxon with index 1, corresponding to the position of the least-significant or right-most bit).

Why can't we use the clade bitmask as the split bitmask for unrooted trees?
This is because we can rely on the clade bitmask to reliably establish split identity across different representation of unrooted trees.
In an urooted tree, any directionality, as given by, e.g. going from a parent node or child node of an edge, is an artifact of the tree representation.
With an unrooted tree, we can only be sure that a node is either internal or a leaf, and distinctions between parent/ancestor nodes and child/descendent nodes are arbitrary and ephemeral.
Because a tree structure requires a parent-child relationship between nodes, this is imposed on an unrooted tree when it is constructed, and for any pair of neighbor nodes, the designation of one as parent and the other as child is contigent on the vagaries of how the tree was serialized, deserialized, and constructed.
Thus, the edge corresponding to a particular split may have a particular child node in one unrooted tree representation, but in *different* representation of the *same* split in the *same* unrooted tree, the same child node may actually be a parent node.
So, 'clade bitmasks' are unstable representations of splits for unrooted trees (but remain accurate and convenient representations of the descendent leaf-sets of nodes in both unrooted and rooted trees).
By normalizing the clade bitmasks to ensure a particular consistency in membership (i.e., group "0" always includes the first taxon), we can usefully compare splits (or the groups induced by splits) across different representations, pseudo-rooting, or rotations of unrooted trees.

Calculating Splits on Trees
---------------------------

Calling :meth:`Tree.encode_splits()` on an object populates the edges of the :class:`Tree` object with two attributes, :attr:`Edge.split_bitmask` and :attr:`Edge.clade_bitmask`. The first is a hash of the split identity, and can be used to identify and compare splits across different trees. The second is a hash of the leaf-set descended from the edge, and can be used to establish ancestor-descendent relationships.

A large number of DendroPy functions calculate the split and clade bitmasks in the background: from tree comparison approaches (e.g., calculating the Robinson-Foulds distance), to working with within-tree operations (e.g., finding the most-recent common ancestor between two nodes or patrisitic distances between taxa), to tree-set operations (e.g., building consensus trees or scoring tree clade credibilities and finding the maximum clade credibility tree).
When passing trees to these methods and functions, these functions will call :meth:`Tree.encode_splits()` automatically for you unless you explicitly specify that this should not be done by passing in '``is_splits_encoded=True``'. You want to do this to avoid the overhead of unnecessarily recalculating the split and clade bitmasks every time, as these operations can simply use the values stored in the :attr:`Edge.split_bitmask` and :attr:`Edge.clade_bitmask` attributes.

The typical usage idiom in this context would be to:

    (1) Establish a common taxon namespace [i.e., creating a global
        :class:`TaxonNamespace` object and pass it in to all
        reading/parsing/input operations]
    (2) Read/load the trees, calling :meth:`Tree.encode_splits()` on each one.
    (3) Perform the calculations, making sure to specify ``is_splits_encoded=True``.

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
    tree1.encode_splits()
    tree2 = dendropy.Tree.get_from_path(
        "mle2.tre",
        "nexus",
        taxon_namespace=taxa)
    tree2.encode_splits()
    counts1 = 0
    counts2
    for sup_tree in dendropy.Tree.yield_from_files(
        files=["boots1.tre", "boots2.tre", "boostraps3.tre"],
        schema="nexus",
        taxon_namespace=taxa):
        sup_tree.encode_splits()
        if treecompare.symmetric_difference(
                tree1, sup_tree, is_splits_encoded=True):
            count1 += 1
        if treecompare.symmetric_difference(
                tree2, sup_tree, is_splits_encoded=True):
            count2 += 2
    print(count1, count2)

Once the split and clade bitmasks have been calculated, as noted above, they are stored in the tree. Some functions inspect the tree to see if the :attr:`Edge.split_bitmask` and :attr:`Edge.clade_bitmask` fields have been populated, and, if so, skip the splits encoding operations themselves.
*If* the tree structure has changed since the splits were last encoded (either explicitly, by you calling :meth:`Tree.encode_splits()` yourself, or implicitly, by a calculation operation), then you have to make sure to update he splits on the tree by re-calling :meth:`Tree.encode_splits()` to avoid errors.

Note that in all cases, for splits to be meaningfully compared two conditions must hold:

    (1) The trees must reference the *same* operational taxonomic unit namespace object, :class:`TaxonNamespace`.
    (2) The trees must have the same rooting state (i.e., all rooted or all unrooted).
