#############
DendroPy News
#############

2018-05-07
==========

DendroPy 4.4.0 has been released!

-   Calculation of birth-death likelihoods.
-   Bipartitions inherit rooting state of trees.
-   Patristic paths between tips can be tracked in ``PatristicDistanceMatrix``.
-   Character column metadata annotations now actually possible.
-   Standard character matrix defaults to 0-9 alphabet instead of just 01.
-   Reorganization of package directory: from "$HOME/dendropy" and "$HOME/dendropy/test" to more modern "$HOME/src/dendropy" and "$HOME/tests" respectively.

For more information:

    http://dendropy.org

Get it now with:

    $ sudo pip install --upgrade dendropy

or:
    $ sudo pip install -U dendropy

You can install directly from the main GitHub repository using:

    $ pip install git+https://github.com/jeetsukumaran/DendroPy.git

or:

    $ pip install git+git://github.com/jeetsukumaran/DendroPy.git

Or, if you are working out of a local Git repo:

    $ git pull origin
    $ python setup.py build
    $ sudo python setup.py develop

(note: ``sudo`` should be used only if you are installing system-wide using your system Python).

2017-06-17
==========

DendroPy 4.3.0 has been released!

Folks, this one includes a VERY important bugfix for SumTrees, so I highly recommend everyone upgrade to this.

Even if you do not use SumTrees, you probably want to grab this for all the goodies that it offers. Apart from the usualy minor tweaks and fixes, you are going to be getting:

-   A new ``sumlabels.py`` application: (string) concatenation of labels for corresponding edges across multiple trees.
-   The birth-death tree simulation (``dendropy.model.birth_death_tree``) now allows for preservation of extinct tips.
-   Improved performance of character subset export.

For more information:

    http://dendropy.org

Get it now with:

    $ sudo pip install --upgrade dendropy

or:
    $ sudo pip install -U dendropy

You can install directly from the main GitHub repository using:

    $ pip install git+https://github.com/jeetsukumaran/DendroPy.git

or:

    $ pip install git+git://github.com/jeetsukumaran/DendroPy.git

Or, if you are working out of a local Git repo:

    $ git pull origin
    $ python setup.py build
    $ sudo python setup.py develop

(note: ``sudo`` should be used only if you are installing system-wide using your system Python).

2016-12-28
==========

DendroPy 4.2.0 has been released!!
Highlights:

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

For more information:

    http://dendropy.org

Get it now with:

    $ sudo pip install --upgrade dendropy

or:
    $ sudo pip install -U dendropy

You can install directly from the main GitHub repository using:

    $ pip install git+https://github.com/jeetsukumaran/DendroPy.git

or:

    $ pip install git+git://github.com/jeetsukumaran/DendroPy.git

Or, if you are working out of a local Git repo:

    $ git pull origin
    $ python setup.py build
    $ sudo python setup.py develop

(note: ``sudo`` should be used only if you are installing system-wide using your system Python).


2016-12-28
==========

DendroPy 4.2.0 has been released!!
Highlights:

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

For more information:

    http://dendropy.org

Get it now with:

    $ sudo pip install --upgrade dendropy

or:
    $ sudo pip install -U dendropy

You can install directly from the main GitHub repository using:

    $ pip install git+https://github.com/jeetsukumaran/DendroPy.git

or:

    $ pip install git+git://github.com/jeetsukumaran/DendroPy.git

Or, if you are working out of a local Git repo:

    $ git pull origin
    $ python setup.py build
    $ sudo python setup.py develop

(note: ``sudo`` should be used only if you are installing system-wide using your system Python).


2016-03-16
==========

DendroPy 4.1.0 has been released!!

Lots of new features, improvements, and bug fixes!

Major new features include:

    -   [SumTrees]: Tip dating! Tip-dating/non-contemporaneous tip age assignment using the "``--tip-ages``" argument (http://dendropy.org/programs/sumtrees.html#setting-the-node-ages-of-the-summary-trees).
    -   [SumTrees]: Collapse clades! "``--min-clade-freq``" applies to all summary targets (i.e., not just consensus trees, but user-specified as well as, e.g. MCCT trees).
    -   Distance trees! Neighbor-joining and UPGMA trees (http://dendropy.org/primer/phylogenetic_distances.html#generating-distance-trees-from-a-phylogeneticdistancematrix-object).
    -   Phylogenetic community ecology statistics! Mean Pairwise Distance (MPD), Mean Nearest Taxon Distance (MNTD), Standardized Effect Size MPD and MNTD, equivalent to -1 * NRI and -1 * NTI (http://dendropy.org/primer/phylogenetic_distances.html#phylogenetic-community-statistics).
    -   The Protracted Speciation birth-death model: a Birth-Death process with explicit modeling of speciation-as-a-process rather than speciation-as-an-event by incorporating the lag between speciation initiation and speciation completion.
    -   Fast, flexible, and powerful tree and subtree cloning, extracting only nodes/taxa of interest (http://dendropy.org/primer/treemanips.html#extracting-trees-and-subtrees-from-an-existing-tree).

*LOTS* more happening under the hood! For more information:

    http://dendropy.org

Get it now with:

    $ sudo pip install --upgrade dendropy

or:
    $ sudo pip install -U dendropy


If you are working out of a Git repo:

    $ git pull origin
    $ python setup.py build
    $ sudo python setup.py develop

(note: ``sudo`` should be used only if you are installing system-wide using your system Python).


2015-04-06
==========

The fourth major version series of the DendroPy Phylogenetic Computing Library has been released!

    http://dendropy.org

Get it now with:

    $ sudo pip install -U dendropy

-   DendroPy 4 runs under Python 2.7 and Python 3.x
-   Re-architectured and re-engineered from the ground up, yet preserving (as much as possible, though certainly not all) the public API of DendroPy 3.x.
-   MAJOR, MAJOR, MAJOR performance improvements in data file reading and processing! Newick and Nexus tree file parsing crazily optimized, with performance scaling at O(1) rather than O(N) or O(n^2) (i.e., in practical terms, you will see better performance improvements with bigger trees when comparing DendroPy 4 vs. DendroPy 3). A thousand-tip tree can be parsed in 0.1 seconds with DendroPy 4 vs. 0.2 seconds with DendroPy 3, while a one million-tip tree can be parsed in under two minutes with DendroPy 4, vs. over 4 days with DendroPy 3. These performance improvements will percolate down to all applications based on DendroPy, including, for example, SumTrees.
-   Tests, tests, tests, tests, and more tests! The core library has a stupendous amount of new tests added, and with each one the ability to zero in and identify, isolate, and deal with bugs is improved.
-   Related to above: dozens of nasty bugs have been dealt with. No, not killed, because we are not that kind of organization. Rather, they have been taken to the big testing farm in the quarantine zone where they can lead healthy lives munching on mock constructs and helping us test the the library to ensure that it works as advertised so that <em>your</em> code works as advertised.
-   Documentation, documentation, documentation! The goal is to have <em>every</em> public method, function, or class fully-documented.
-   Many, many, many, many new features: e.g., a high-performance TreeArray class, calculation of MCCT topologies, new simulation models, new tree statistics, new tree manipulation routines.
-   SumTrees works faster than ever before thanks to the above improvements, and also allows for many new operations such as rerooting the target tree, using an MCCT tree as the target topology, extensive extra information summarized, auto-detection of number of parallel processors etc.: http://dendropy.org/programs/sumtrees.html .
-   The newly rewritten DendroPy primer is just full of information to get you started: http://dendropy.org/primer/index.html .
-   The "work-in-progress" migration primer will help ease the transition from 3 to 4: http://dendropy.org/migration.html .
-   Comprehensive documentation of all the data formats supported, plus all the keyword arguments you can use to control and customize reading and writing in all these different formats: http://dendropy.org/schemas/index.html .
-   A glossary of terms, to clarify the simultaneously redundant and oversubscribed/conflicting terminological soup that characterizes a lot of phylogenetics: http://dendropy.org/glossary.html .

