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

If you have an existing |Tree| instance, e.g. ``tree``, the :meth:`dendropy.datamodel.treemodel.Tree.phylogenetic_distance_matrix` method returns a "snapshot" of the tip taxon-to-taxon distances given on the tree:

.. literalinclude:: /examples/pdm.py





