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
