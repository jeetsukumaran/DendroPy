Release 4.5.2
-------------

-   Support for user-specified random seed in RaXML wrapper (thanks @NoahAmsel)
-   *MUCH* faster label lookup (thanks Sam Nicholls / @SamStudio8 !)
-   Faster birth-death tree generation (thanks @NicolaDM !)
-   Storage of supplemental NEXUS blocks
-   Fix type: "PhylogeneticIndependentConstrasts" => "PhylogeneticIndependentContrasts"

Release 4.4.0
-------------

-   Calculation of birth-death likelihoods.
-   Bipartitions inherit rooting state of trees.
-   Patristic paths between tips can be tracked in ``PatristicDistanceMatrix``.
-   Character column metadata annotations now actually possible.
-   Standard character matrix defaults to 0-9 alphabet instead of just 01.
-   Reorganization of package directory: from "$HOME/dendropy" and "$HOME/dendropy/test" to more modern "$HOME/src/dendropy" and "$HOME/tests" respectively.

Release 4.3.0
-------------

-   [SumTrees]: Important bugfix in tracking splits on consensus tree from rooted trees: previously it was possible for a split on the consensus tree to be ignored, resulting in a null (0) edge length and no (0.0) support.
-   Added ``sumlabels.py`` application.
-   Birth-death tree (``dendropy.model.birth_death_tree``) now allows for preservation of extinct tips.
-   Improved performance of character subset export

Release 4.2.0
-------------

-   0 branch lengths assigned to randomly resolved polytomies.
-   Explicitly set rooting for NJ and UPGMA trees.
-   Faster pruning (kyungtaekLIM)
-   Fix nesting bug in raised KeyError in basemodel.AnnotationSet.__deepcopy__ (Steve Bond)
-   Catch edge case during deepcopy when Edge object has no _annotations (Steve Bond)
-   Optimizations and fixes for various population genetic calculations (Andrew Guy)
-   newickreader: Parse jplace style edge numbering. (Ben J Woodcroft)
-   Calculate probability of gene tree(s) in species trees under the Multispecies Coalescent model.
-   New approaches to calculate distances between unlabeled trees of different sizes: ``dendropy.profiledistance`` and ``dendropy.calculate.treecompare.TreeShapeKernel``.
-   When parsing Newick/NEXUS, allow for internal node labels to be associated with either nodes or edges.

Release 4.1.0
-------------

New or Updated Features
^^^^^^^^^^^^^^^^^^^^^^^

    -   [SumTrees]: tip-dating/non-contemporaneous tip age assignment using the "``--tip-ages``" argument (http://dendropy.org/programs/sumtrees.html#setting-the-node-ages-of-the-summary-trees).
    -   [SumTrees]: "``--min-clade-freq``" applies to all summary targets (i.e., not just consensus trees, but user-specified as well as, e.g. MCCT trees).
    -   Fast, flexible, and powerful tree and subtree cloning, extracting only nodes/taxa of interest (http://dendropy.org/primer/treemanips.html#extracting-trees-and-subtrees-from-an-existing-tree).
    -   Neighbor-joining and UPGMA trees (http://dendropy.org/primer/phylogenetic_distances.html#generating-distance-trees-from-a-phylogeneticdistancematrix-object).
    -   The new (actually, warmed-over) PhylogeneticDistanceMatrix to manage various "within-tree" distances, such patristic distances, or the ecological statistics described below (http://dendropy.org/primer/phylogenetic_distances.html#creating-a-phylogeneticdistancematrix-object).
    -   Added phylogenetic community ecology statistic calculations: Mean Pairwise Distance (MPD), Mean Nearest Taxon Distance (MNTD), Standardized Effect Size MPD and MNTD, equivalent to -1 * NRI and -1 * NTI (http://dendropy.org/primer/phylogenetic_distances.html#phylogenetic-community-statistics).
    -   Added DataTable class to manage community ecology (as well as more general classes of) data.
    -   Implementation of the Protracted Speciation model: a Birth-Death process with explicit modeling of speciation-as-a-process rather than speciation-as-an-event by incorporating the lag between speciation initiation and speciation completion.
    -   NEWICK terminating semicolon requirement relaxation.
    -   Some more refined node filtering/dropping.
    -   Return list of nodes dropped when filtering out leaves.
    -   Force max/min ages when calculating node ages; and beginning of support for setting node ages by function.
    -   Implementation of Tree.find_nodes() to return collection of nodes that match instaed of just the first one.

Bug Fixes
^^^^^^^^^

    -   Handle sequence comparison where there are no non-ignored sites in common.
    -   Update string type checking to handle unicode etc. under Python 2.
    -   Exclusion of trees from data set reads actually works.
    -   Actually implement symbol to state (alphabet) identity coercion in derived classes.
    -   Pop out inner classes to enable pickling.
    -   Several bugs, mostly caused by leftovers of DendroPy3 code.
    -   Made group_ranges work properly with unordered iterables.
    -   Make PHYLIP writing work correctly with missing taxa.


Release 4.0.3
-------------

Bug Fixes
^^^^^^^^^

    -   [SumTrees]: propagate ``-f``/``--frequency`` option to underlying summarization engine.
    -   [SumTrees]: ``-v``/``--ultrametricity-precision`` option takes numeric value.
    -   Exporting of characters from matrix suppresses cloned character subset definitions.

Release 4.0.2
-------------

Bug Fixes
^^^^^^^^^

    -   Adjustment of child edge lengths when collapsing basal bifurcations.

Release 4.0.1
-------------

Bug Fixes
^^^^^^^^^

    -   Fix for installing using in virtual environments under ``virtualenv``.

