*****
NEXUS
*****

.. contents::
    :local:
    :backlinks: none

Description
===========

Reading
=======

.. _schema_specific_keyword_arguments_reading_nexus:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.nexusreader.NexusReader.__init__

Supported Methods
-----------------

``Tree.get``
............
(:meth:`method reference <dendropy.datamodel.treemodel.Tree.get>`)

.. literalinclude:: /schemas/interfaces/nexus_tree_get.py

``TreeList.get``
................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeList.get>`)

.. literalinclude:: /schemas/interfaces/nexus_treelist_get.py

``TreeList.read``
.................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeList.read>`)

.. literalinclude:: /schemas/interfaces/nexus_treelist_read.py

``TreeArray.read``
..................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeArray.read>`)

.. literalinclude:: /schemas/interfaces/nexus_treearray_read.py

``DnaCharacterMatrix.get``
..........................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/nexus_dnacharactermatrix_get.py

``RnaCharacterMatrix.get``
..........................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/nexus_rnacharactermatrix_get.py

``ProteinCharacterMatrix.get``
..............................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/nexus_proteincharactermatrix_get.py

``RestrictionSitesCharacterMatrix.get``
.......................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/nexus_restrictionsitescharactermatrix_get.py

``InfiniteSitesCharacterMatrix.get``
....................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/nexus_infinitesitescharactermatrix_get.py

``StandardCharacterMatrix.get``
...............................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/nexus_standardcharactermatrix_get.py

``DataSet.get``
...............
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.get>`)

.. literalinclude:: /schemas/interfaces/nexus_dataset_get.py

``DataSet.read``
................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.read>`)

.. literalinclude:: /schemas/interfaces/nexus_dataset_read.py

Writing
=======

.. _schema_specific_keyword_arguments_writing_nexus:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.nexuswriter.NexusWriter.__init__

Supported Methods
-----------------

``Tree.write``
..............
(:meth:`method reference <dendropy.datamodel.treemodel.Tree.write>`)

.. literalinclude:: /schemas/interfaces/nexus_trees_write.py

``Tree.as_string``
..................
(:meth:`method reference <dendropy.datamodel.treemodel.Tree.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_trees_as_string.py

``TreeList.write``
..................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeList.write>`)

.. literalinclude:: /schemas/interfaces/nexus_trees_write.py

``TreeList.as_string``
......................
(:meth:`method reference <dendropy.datamodel.treecollectionmodel.TreeList.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_trees_as_string.py

``DnaCharacterMatrix.write``
.............................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_write.py

``DnaCharacterMatrix.as_string``
................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_as_string.py


``RnaCharacterMatrix.write``
.............................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_write.py

``RnaCharacterMatrix.as_string``
................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_as_string.py


``ProteinCharacterMatrix.write``
.................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_write.py

``ProteinCharacterMatrix.as_string``
.....................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_as_string.py


``RestrictionSitesCharacterMatrix.write``
.........................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_write.py

``RestrictionSitesCharacterMatrix.as_string``
.............................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_as_string.py


``InfiniteSitesCharacterMatrix.write``
......................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_write.py

``InfiniteSitesCharacterMatrix.as_string``
..........................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_as_string.py


``StandardCharacterMatrix.write``
.................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_write.py

``StandardCharacterMatrix.as_string``
.....................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_chars_as_string.py


``DataSet.write``
.................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.write>`)

.. literalinclude:: /schemas/interfaces/nexus_dataset_write.py

``DataSet.as_string``
.....................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.as_string>`)

.. literalinclude:: /schemas/interfaces/nexus_dataset_as_string.py



