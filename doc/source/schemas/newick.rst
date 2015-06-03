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

.. _schema_specific_keyword_arguments_reading_newick:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.newickreader.NewickReader.__init__

Supported Methods
-----------------

``Tree.get``
............
(:meth:`method reference <dendropy.datamodel.treemodel.Tree.get>`)

.. literalinclude:: /schemas/interfaces/newick_tree_get.py

``TreeList.get``
................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeList.get>`)

.. literalinclude:: /schemas/interfaces/newick_treelist_get.py

``TreeList.read``
.................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeList.read>`)

.. literalinclude:: /schemas/interfaces/newick_treelist_read.py

``TreeArray.read``
..................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeArray.read>`)

.. literalinclude:: /schemas/interfaces/newick_treearray_read.py

``DataSet.get``
...............
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.get>`)

.. literalinclude:: /schemas/interfaces/newick_dataset_get.py

``DataSet.read``
................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.read>`)

.. literalinclude:: /schemas/interfaces/newick_dataset_read.py

Writing
=======

.. _schema_specific_keyword_arguments_writing_newick:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.newickwriter.NewickWriter.__init__

Supported Methods
-----------------

``Tree.write``
..............
(:meth:`method reference <dendropy.datamodel.treemodel.Tree.write>`)

.. literalinclude:: /schemas/interfaces/newick_write.py

``Tree.as_string``
..................
(:meth:`method reference <dendropy.datamodel.treemodel.Tree.as_string>`)

.. literalinclude:: /schemas/interfaces/newick_as_string.py

``TreeList.write``
..................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeList.write>`)

.. literalinclude:: /schemas/interfaces/newick_write.py

``TreeList.as_string``
......................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeList.as_string>`)

.. literalinclude:: /schemas/interfaces/newick_as_string.py

``DataSet.write``
.................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.write>`)

.. literalinclude:: /schemas/interfaces/newick_write.py

``DataSet.as_string``
.....................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.as_string>`)

.. literalinclude:: /schemas/interfaces/newick_as_string.py


