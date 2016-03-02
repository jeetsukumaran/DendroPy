**************************************
Phylogenetic Distance-Based Operations
**************************************

The |PhylogeneticDistanceMatrix| class is a comprehensive class for tracking taxon-to-taxon distances and operations associated with these values.


Creating a |PhylogeneticDistanceMatrix| Object
===============================================

There are basically two ways to create a |PhylogeneticDistanceMatrix|
instance:

1.  Create it from an existing tree based on the patristic distances on the tree
2.  Create it from an external data source specifying the distances between taxa.

Creating a |PhylogeneticDistanceMatrix| Object From a |Tree|
------------------------------------------------------------

If you have an existing |Tree| instance, e.g. ``tree``, the :meth:`~dendropy.datamodel.treemodel.Tree.phylogenetic_distance_matrix` method returns a "snapshot" of the tip taxon-to-taxon distances given on the tree:

.. literalinclude:: /examples/pdm.py

The new |PhylogeneticDistanceMatrix| object will reference the same |TaxonNamespace| and member |Taxon| objects as the tree, ``tree``.

Creating a |PhylogeneticDistanceMatrix| Object From an External Data Source
---------------------------------------------------------------------------

The :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.from_csv` method reads a token-delimited external data source specifying taxon-to-taxon distances and creates and returns a corresponding |PhylogeneticDistanceMatrix|.
This data source is expected to provide a table where each row is a separate line and each column is separated from the preceding by a token (typically a comma or a tab character).
The cells of the table are numeric (typically real) values that indicate the
distance between the taxa of the current row and column. Note that *only* the
upper right section of the table is considered. The diagonals values are
typically zeroes and, in either case, ignored along with the lower diagonal.
Despite being ignored by the PhylogeneticDistanceMatrix object, the values are
parsed by the underlying reader and thus have to be valid numerical values.

For example,  the data from Table 1 of Saitou and Nei (1987) can be represented by the following comma-separated value (CSV) file::

    ,a,b,c,d,e,f,g,h
    a,0,7,8,11,13,16,13,17
    b,7,0,5,8,10,13,10,14
    c,8,5,0,5,7,10,7,11
    d,11,8,5,0,8,11,8,12
    e,13,10,7,8,0,5,6,10
    f,16,13,10,11,5,0,9,13
    g,13,10,7,8,6,9,0,8
    h,17,14,11,12,10,13,8,0

(Note the empty cell of the first column of the first row).

This file can be instantiated into a |PhylogeneticDistanceMatrix| by the following, which specifies the keyword argument ``delimiter=","`` for correct parsing::

    import dendropy
    pdm1 = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src=open("data.csv"),
            delimiter=",")

If tabs were used as a delimiter instead, then ``delimiter="\t"`` would be specified::

    import dendropy
    pdm1 = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src=open("data.tsv"),
            delimiter="\t")

In the above examples, the taxon namespace into which the data was read was a new one, created by default.
In many cases, e.g. if you want to compare the data to data from other sources, you will want to ensure that the |TaxonNamespace| is shared. You do this by passing in a |TaxonNamespace| instance using the ``taxon_namespace`` argument:

.. literalinclude:: /examples/pdm_tns1.py

Note that if you do pass in a non-empty |TaxonNamespace| value using the ``taxon_namespace`` argument, by default the reading process protects any new |Taxon| objects from being created in this taxon namespace. This is to protect the integrity of the taxon namespace from errors due to inadverent label mismatches (see note below).
You can specify ``is_allow_new_taxa=True`` to relax this restriction.

.. note::

    It is important that the taxon labels specified in the data source match the taxon labels of the |Taxon| objects with which you want to associate the data *exactly*. This includes taking into account format-based transformations. For example, by default, unless ``preserve_underscores=True`` is specified, when reading Newick and NEXUS format data underscores not protected by quotes will be translated into spaces. If the taxon labels in the distance file have underscores in them, however, by default these will *not* be translated into spaces. Thus the taxon "Python_regius" in a Newick/NEXUS tree file or a NEXUS character data file will be represented by "Python regius" internally, and this will not match the taxon label "Python_regius" given in the distance table. The solution is to ensure that the distance table file has spaces in place of underscores or use ``preserve_underscores=True`` when reading Newick/NEXUS data. For the former, in many cases you might be able to  pass in a function to transform or translate the distance data table labels using the ``label_transform_fn`` argument of the :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.from_csv` method. E.g.::

        label_transform_fn = lambda x: x.replace("_", " ")
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
                src=open("pythonidae.mle.weighted.pdm.csv"),
                delimiter=",",
                label_transform_fn=label_transform_fn)

Calculating Patristic Distances and Most-Recent Common Ancestors (MRCA)
=======================================================================

.. literalinclude:: /examples/mrca2.py

.. note:: "Weighted" distances (or "weighted edge" distances) refers to distances taking the edge weights or branch lengths into account. "Unweighted" distances (or "unweighted edge") distances refers to distances considering only the total *number* of edges connecting two taxa, rather than the total branch length.

Generating Distance Trees from a |PhylogeneticDistanceMatrix| Object
====================================================================

Neighbor-Joining Trees
----------------------

The :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.nj_tree` method returns a |Tree| representing the neighbor-joining tree calculated on the distances in the matrix:

.. literalinclude:: /examples/pdm_nj_tree.py

UPGMA Trees
-----------

The :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.upgma_tree` method returns a |Tree| representing the UPGMA tree calculated on the distances in the matrix:

.. literalinclude:: /examples/pdm_upgma_tree.py

Phylogenetic Community Statistics
=================================

Basic Phylogenetic Community Statistics
---------------------------------------

Various phylogenetic community statistics can be calculated for one or more definitions of "community".


    -   The Mean Pairwise Distance (MPD) is returned by the :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.mean_pairwise_distance` method, which calculates:

        .. math::

            mpd = \frac{ \sum_{i}^{n} \sum_{j}^{n} \delta_{i,j} }{n \choose 2},

        where :math:`i \neq j`, :math:`\delta_{i,j}` is the phylogenetic
        distance between species :math:`i` and :math:`j`, and :math:`n` is the number
        of species in the sample.


    -   The Mean Nearest Taxon Distance (MNTD) is returned by the :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.mean_nearest_taxon_distance` method, which calculates:

            .. math::
                mntd = \frac{ \sum_{i}^{n} min(\delta_{i,j}) }{n},

        where :math:`i \neq j`, :math:`\delta_{i,j}` is the phylogenetic
        distance between species :math:`i` and :math:`j`, and :math:`n` is the number
        of species in the sample.

Each of these methods takes a function object as a ``filter_fn`` argument.
This function object serves to filter the taxa of the tree, reducing it to so that the tips are restricted to the community or assemblage of interest.
If not specified, then all leaves are considered in the calculation.
If specified, the function object should take a |Taxon| object as its only argument and return ``True`` if the |Taxon| is considered part of the assembalge or community or |False| if not.
For example, to calculate the MPD of the entire tree and then of some (highly artificial) communities:

.. literalinclude:: /examples/pdm_mpd0.py

A more realistic example is where a tree is sampled across multiple communities, with the data read from a tab-delimited source:

.. literalinclude:: /examples/pdm_mpd1.py

which results in::

    Assemblage Memberships:
    C1: ['spA', 'spC', 'spE', 'spG', 'spM', 'spN', 'spO']
    C2: ['spB', 'spD', 'spE', 'spG', 'spJ', 'spK']
    C3: ['spA', 'spC', 'spF', 'spG', 'spH', 'spN']
    C4: ['spH', 'spI', 'spJ', 'spK', 'spL', 'spM']
    C5: ['spH', 'spI', 'spJ']
    Phylogenetic Community Statistics:
    C1: MPD=3.19428571429, MNTD=1.61142857143
    C2: MPD=1.88666666667, MNTD=1.06666666667
    C3: MPD=1.88666666667, MNTD=1.06166666667
    C4: MPD=1.91066666667, MNTD=0.108333333333
    C5: MPD=0.18, MNTD=0.13

Standardized Effect Size Statistics
-----------------------------------

The Standardized Effect Size (S.E.S.) of these statistics are useful to remove any bias associated with the decrease in variance as species richness increases to the point where assemblages are saturated.
The S.E.S. is calculated under a null model, given here by random shuffling of the tip labels of the tree.
The statistic is calculated for each randomization of the tip labels, and the mean and standard deviation of the collection of values across multiple randomization replicates is used to calculate the S.E.S
For a particular statistic, the S.E.S. is obtained by dividing the difference between the value of the statistic as given by the original data and the mean of the statistic across the replicates generated under the null model, divided by the standard deviation of the statistic across the replicates generated under the null model:

.. math::
    SES(statistic) = \frac{observed - mean(model_{null})}{sd(model_{null})}

The S.E.S. of the MPD is calculated by the :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.standardized_effect_size_mean_pairwise_distance` method, while the S.E.S of the MNTD is calculated by the :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.standardized_effect_size_mean_nearest_taxon_distance` method.

Instead of a filter function, as with the basic statistics above, these methods take an ``assemblage_memberships`` argument, the value of which should be an iterable of iterable of |Taxon| objects. For e.g., a list of sets of |Taxon| objects, where each set in the list specifies the membership of a single assemblage.
The return value of these methods is a list of ``namedtuple`` objects, with each element in the list the result of associated with community/assemblage definition in the corresponding position of the input ``assemblage_memberships`` list.
The result ``namedtuple`` objects have the following fields:

    ``obs``
        the observed value of the statistic
    ``null_model_mean``
        the mean value of the statistic under the null model
    ``null_model_sd``
        the standard deviation of the statistic under the null model
    ``z``
        the standardized effect (S.E.S.) value of the statistic
    ``p``
        the p-value of the observed value of the statistic

As an example:

.. literalinclude:: /examples/pdm_ses1.py

which results in::

    Phylogenetic Community Standardized Effect Size Statistics:
    # Assemblage 'C1' (['spA', 'spC', 'spE', 'spG', 'spM', 'spN', 'spO'])
    -     MPD: 3.19428571429
    - SES MPD: 1.42344634503
    - p-value: 0.982
    # Assemblage 'C2' (['spB', 'spD', 'spE', 'spG', 'spJ', 'spK'])
    -     MPD: 1.88666666667
    - SES MPD: -0.917064662164
    - p-value: 0.21
    # Assemblage 'C3' (['spA', 'spC', 'spF', 'spG', 'spH', 'spN'])
    -     MPD: 1.88666666667
    - SES MPD: -0.769722690565
    - p-value: 0.24
    # Assemblage 'C4' (['spH', 'spI', 'spJ', 'spK', 'spL', 'spM'])
    -     MPD: 1.91066666667
    - SES MPD: -0.87070720087
    - p-value: 0.229
    # Assemblage 'C5' (['spH', 'spI', 'spJ'])
    -     MPD: 0.18
    - SES MPD: -1.99810309955
    - p-value: 0.0

If you are unhappy with the extra book-keeping involved with co-ordinating all these different lists (``assemblage_memberships``, ``assemblage_names``, etc.), you can associate at least some of these lists in an ``OrderedDict`` object:

.. literalinclude:: /examples/pdm_ses2.py

A convenience method is available to read community data from a delimited source, :meth:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.assemblage_membership_definitions_from_csv`, which makes the process somewhat easier:

.. literalinclude:: /examples/pdm_ses3.py

