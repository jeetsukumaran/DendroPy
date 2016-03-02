**************************************
Phylogenetic Distance-Based Operations
**************************************

The |PhylogeneticDistanceMatrix| class is a comprehensive class for tracking taxon-to-taxon distances and operations associated with these values.


Generating a |PhylogeneticDistanceMatrix| Object
================================================

There are basically two ways to create a |PhylogeneticDistanceMatrix|
instance:

1.  Create it from an existing tree based on the patristic distances on the tree
2.  Create it from an external data source specifying the distances between taxa.

Creating a |PhylogeneticDistanceMatrix| Object From a |Tree|
------------------------------------------------------------

If you have an existing |Tree| instance, e.g. ``tree``, the :meth:`~dendropy.datamodel.treemodel.Tree.phylogenetic_distance_matrix` method returns a "snapshot" of the tip taxon-to-taxon distances given on the tree:

.. literalinclude:: /examples/pdm.py

Reading a |PhylogeneticDistanceMatrix| Object From an External Data Source
--------------------------------------------------------------------------

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


