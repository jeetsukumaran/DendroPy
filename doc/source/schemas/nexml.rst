*****
NeXML
*****

.. contents::
    :local:
    :backlinks: none

Description
===========

Reading
=======

.. _schema_specific_keyword_arguments_reading_nexml:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.nexmlreader.NexmlReader.__init__

Supported Methods
-----------------

:meth:`Tree.get <dendropy.datamodel.treemodel.Tree.get>`
........................................................
.. literalinclude:: /schemas/interfaces/nexml_tree_get.py

:meth:`TreeList.get <dendropy.datamodel.treecollectionmodel.TreeList.get>`
..........................................................................
.. literalinclude:: /schemas/interfaces/nexml_treelist_get.py

:meth:`TreeList.read <dendropy.datamodel.treecollectionmodel.TreeList.read>`
............................................................................
.. literalinclude:: /schemas/interfaces/nexml_treelist_read.py

:meth:`TreeArray.read <dendropy.datamodel.treecollectionmodel.TreeArray.read>`
..............................................................................
.. literalinclude:: /schemas/interfaces/nexml_treearray_read.py

:meth:`DnaCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.get>`
..........................................................................................

.. literalinclude:: /schemas/interfaces/nexml_dnacharactermatrix_get.py

:meth:`RnaCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.get>`
..........................................................................................

.. literalinclude:: /schemas/interfaces/nexml_rnacharactermatrix_get.py

:meth:`ProteinCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.get>`
..................................................................................................

.. literalinclude:: /schemas/interfaces/nexml_proteincharactermatrix_get.py

:meth:`RestrictionSitesCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.get>`
....................................................................................................................

.. literalinclude:: /schemas/interfaces/nexml_restrictionsitescharactermatrix_get.py

:meth:`InfiniteSitesCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.get>`
...............................................................................................................

.. literalinclude:: /schemas/interfaces/nexml_infinitesitescharactermatrix_get.py

:meth:`StandardCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.get>`
....................................................................................................

.. literalinclude:: /schemas/interfaces/nexml_standardcharactermatrix_get.py

:meth:`DataSet.get <dendropy.datamodel.datasetmodel.DataSet.get>`
.................................................................
.. literalinclude:: /schemas/interfaces/nexml_dataset_get.py

:meth:`DataSet.read <dendropy.datamodel.datasetmodel.DataSet.read>`
...................................................................
.. literalinclude:: /schemas/interfaces/nexml_dataset_read.py

Writing
=======

.. _schema_specific_keyword_arguments_writing_nexml:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.nexmlwriter.NexmlWriter.__init__

Supported Methods
-----------------

:meth:`Tree.write <dendropy.datamodel.treemodel.Tree.write>`
............................................................
.. literalinclude:: /schemas/interfaces/nexml_trees_write.py

:meth:`TreeList.write <dendropy.datamodel.treecollectionmodel.TreeList.write>`
..............................................................................
.. literalinclude:: /schemas/interfaces/nexml_trees_write.py

:meth:`DnaCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.write>`
..............................................................................................

.. literalinclude:: /schemas/interfaces/nexml_chars_write.py

:meth:`RnaCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.write>`
..............................................................................................

.. literalinclude:: /schemas/interfaces/nexml_chars_write.py

:meth:`ProteinCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.write>`
......................................................................................................

.. literalinclude:: /schemas/interfaces/nexml_chars_write.py

:meth:`RestrictionSitesCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.write>`
........................................................................................................................

.. literalinclude:: /schemas/interfaces/nexml_chars_write.py

:meth:`InfiniteSitesCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.write>`
..................................................................................................................

.. literalinclude:: /schemas/interfaces/nexml_chars_write.py

:meth:`StandardCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.write>`
........................................................................................................

.. literalinclude:: /schemas/interfaces/nexml_chars_write.py

:meth:`DataSet.write <dendropy.datamodel.datasetmodel.DataSet.write>`
.....................................................................

.. literalinclude:: /schemas/interfaces/nexml_dataset_write.py




