*******************************************
Splits, Split Bitmasks and Leafset Bitmasks
*******************************************

Splits or Bipartitions are a Division of Leaves into Two Groups
---------------------------------------------------------------

A "split" is a bipartition of a leaf or tip set of a tree, i.e., the division or sorting of the leaves/tips of a tree into two mutually-exclusive and collectively-comprehensive subsets.
The term "bipartition" is often also used to refer to the same concept. For example, the Mr. Bayes documentation says that: "taxon bipartition is defined by removing a branch on the tree, dividing the tree into those taxa to the left and right of the removed branch. This set is called a taxon bipartition"

Splits Correspond to Edges and Splits as Partitions of Taxa
-----------------------------------------------------------

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

Split Bitmasks and Leafset Bitmasks
-----------------------------------

In DendroPy, a "split bitmask" is an attribute of a tree edge that describes the split represented by that edge (i.e., the partition of the taxa induced by bisecting the tree on that edge).
A "leafset bitmask", on the other hand, is an attribute of a tree *node* that describes the taxa associated with the descendent leaves of that node.
For rooted trees, the split bitmask of an edge will be identical to the leafset bitmask of the node subtending that edge.
For unrooted trees, however, the split bitmask of an edge may not be the same as the leafset bitmask of the node.
This is because the split bitmasks of unrooted trees are *normalized* to always have a least-significant bit of 0. The choice of 0 as opposed to 1 is arbitrary, but the reason is so ensure that we can have consistent comparisons of groups across trees of different rotations (and "pseudo-rootings" created by the constraints of tree representation in, e.g., the NEWICK format) by enforcing the convention that group "0" will always be the group that includes the first taxon (i.e., the taxon with index 1, corresponding to the position of the least-significant or right-most bit).

Thus, the child node of particular edge may have the first three taxa as ultimate descendents (i.e., the leafset bitmask of its daughter nodes is "00111"). In a rooted tree, which has a distincting directionality from a real root to tip, that is all there is to it. However, in an urooted tree, there is no such distinct directionality, and the node that "child node" of the edge may easily become a "parent node" in another representation of the *same* unrooted topology, in which case the leafset bitmask would be "11000". If we did not normalize the split bitmasks, these two splits would be considered different, when, for unrooted topologies they are actually the same. By adopting the convention of enforcing the least-significant bit to be 0 in unrooted tree split bitmasks, i.e., that the first-indexed taxon is always in group "0", the split bitmask of "00111" is normalized to "11100": the memberships of the partition subsets remain the same, with the same taxa on either side of the bipartiton; only the arbitrary label assignment of "1" or "0" to the subsets is switched.

To see why, let us consider the following trees:


    (A,(B,(C,(D,E))));
    ((((A,B),C),E),D);

These two topologies are actually identical except for the rooting position: when interpreted as unrooted trees, the symmetric distance between them is 0.


Let us consider the split that separates {A,B,C} from {D,E}, i.e. "00111".


