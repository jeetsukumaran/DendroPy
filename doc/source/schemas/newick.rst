******
NEWICK
******

.. contents::
    :local:
    :backlinks: none

Reading
=======

.. _schema_specific_keyword_arguments_reading_newick:

Supported Methods
-----------------

:meth:`Tree.get <dendropy.datamodel.treemodel.Tree.get>`
........................................................
.. literalinclude:: /schemas/interfaces/newick_tree_get.py

:meth:`TreeList.get <dendropy.datamodel.treemodel.TreeList.get>`
................................................................
.. literalinclude:: /schemas/interfaces/newick_treelist_get.py

:meth:`TreeList.read <dendropy.datamodel.treemodel.TreeList.read>`
..................................................................
.. literalinclude:: /schemas/interfaces/newick_treelist_read.py

:meth:`TreeArray.read <dendropy.datamodel.treemodel.TreeArray.read>`
....................................................................
.. literalinclude:: /schemas/interfaces/newick_treearray_read.py

:meth:`DataSet.get <dendropy.datamodel.treemodel.DataSet.get>`
..............................................................
.. literalinclude:: /schemas/interfaces/newick_dataset_get.py

:meth:`DataSet.read <dendropy.datamodel.treemodel.DataSet.read>`
................................................................
.. literalinclude:: /schemas/interfaces/newick_dataset_read.py

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.newickreader.NewickReader.__init__
