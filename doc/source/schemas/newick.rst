******
Newick
******

.. contents::
    :local:
    :backlinks: none

Description
===========

    * http://evolution.genetics.washington.edu/phylip/newicktree.html
    * http://en.wikipedia.org/wiki/Newick_format
    * http://evolution.genetics.washington.edu/phylip/newick_doc.html

Reading
=======

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

:meth:`TreeArray.read <dendropy.datamodel.treecollectionmodel.TreeArray.read>`
..............................................................................
.. literalinclude:: /schemas/interfaces/newick_treearray_read.py

:meth:`DataSet.get <dendropy.datamodel.datasetmodel.DataSet.get>`
.................................................................
.. literalinclude:: /schemas/interfaces/newick_dataset_get.py

:meth:`DataSet.read <dendropy.datamodel.datasetmodel.DataSet.read>`
...................................................................
.. literalinclude:: /schemas/interfaces/newick_dataset_read.py

.. _schema_specific_keyword_arguments_reading_newick:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.newickreader.NewickReader.__init__

Writing
=======

Supported Methods
-----------------

:meth:`Tree.write <dendropy.datamodel.treemodel.Tree.write>`
............................................................
.. literalinclude:: /schemas/interfaces/newick_write.py

:meth:`TreeList.write <dendropy.datamodel.treecollectionmodel.TreeList.write>`
..............................................................................
.. literalinclude:: /schemas/interfaces/newick_write.py

:meth:`DataSet.write <dendropy.datamodel.datasetmodel.DataSet.write>`
.....................................................................
.. literalinclude:: /schemas/interfaces/newick_write.py

.. _schema_specific_keyword_arguments_writing_newick:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.newickwriter.NewickWriter.__init__

