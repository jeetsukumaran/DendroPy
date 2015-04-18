******
PHYLIP
******

.. contents::
    :local:
    :backlinks: none

Description
===========

Reading
=======

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.phylipreader.PhylipReader.__init__

Supported Methods
-----------------

:meth:`DnaCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.get>`
..........................................................................................

.. literalinclude:: /schemas/interfaces/phylip_dnacharactermatrix_get.py

:meth:`RnaCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.get>`
..........................................................................................

.. literalinclude:: /schemas/interfaces/phylip_rnacharactermatrix_get.py

:meth:`ProteinCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.get>`
..................................................................................................

.. literalinclude:: /schemas/interfaces/phylip_proteincharactermatrix_get.py

:meth:`RestrictionSitesCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.get>`
....................................................................................................................

.. literalinclude:: /schemas/interfaces/phylip_restrictionsitescharactermatrix_get.py

:meth:`InfiniteSitesCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.get>`
...............................................................................................................

.. literalinclude:: /schemas/interfaces/phylip_infinitesitescharactermatrix_get.py

:meth:`StandardCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.get>`
....................................................................................................

.. literalinclude:: /schemas/interfaces/phylip_standardcharactermatrix_get.py

:meth:`DataSet.get <dendropy.datamodel.datasetmodel.DataSet.get>`
.................................................................

Note that the type of data needs to be specified using the ``data_type``
keyword argument.

.. literalinclude:: /schemas/interfaces/phylip_dataset_get.py

:meth:`DataSet.read <dendropy.datamodel.datasetmodel.DataSet.read>`
...................................................................

Note that the type of data needs to be specified using the ``data_type``
keyword argument.

.. literalinclude:: /schemas/interfaces/phylip_dataset_read.py

Writing
=======

.. _schema_specific_keyword_arguments_writing_phylip:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.phylipwriter.PhylipWriter.__init__

Supported Methods
-----------------

:meth:`DnaCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.write>`
..............................................................................................

.. literalinclude:: /schemas/interfaces/phylip_write.py

:meth:`RnaCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.write>`
..............................................................................................

.. literalinclude:: /schemas/interfaces/phylip_write.py

:meth:`ProteinCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.write>`
......................................................................................................

.. literalinclude:: /schemas/interfaces/phylip_write.py

:meth:`RestrictionSitesCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.write>`
........................................................................................................................

.. literalinclude:: /schemas/interfaces/phylip_write.py

:meth:`InfiniteSitesCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.write>`
..................................................................................................................

.. literalinclude:: /schemas/interfaces/phylip_write.py

:meth:`StandardCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.write>`
........................................................................................................

.. literalinclude:: /schemas/interfaces/phylip_write.py

:meth:`DataSet.write <dendropy.datamodel.datasetmodel.DataSet.write>`
.....................................................................

.. literalinclude:: /schemas/interfaces/phylip_write.py


