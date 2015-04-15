*****
FASTA
*****

.. contents::
    :local:
    :backlinks: none

Description
===========

    * http://en.wikipedia.org/wiki/FASTA_format
    * Lipman, DJ; Pearson, WR (1985). "Rapid and sensitive protein similarity searches". Science 227 (4693): 1435–41. `doi:10.1126/science.2983426 <http://www.sciencemag.org/content/227/4693/1435>`_.

Reading
=======

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.fastareader.FastaReader.__init__

Supported Methods
-----------------

:meth:`DnaCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.get>`
..........................................................................................

.. literalinclude:: /schemas/interfaces/fasta_dnacharactermatrix_get.py

:meth:`RnaCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.get>`
..........................................................................................

.. literalinclude:: /schemas/interfaces/fasta_rnacharactermatrix_get.py

:meth:`ProteinCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.get>`
..................................................................................................

.. literalinclude:: /schemas/interfaces/fasta_proteincharactermatrix_get.py

:meth:`RestrictionSitesCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.get>`
....................................................................................................................

.. literalinclude:: /schemas/interfaces/fasta_restrictionsitescharactermatrix_get.py

:meth:`InfiniteSitesCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.get>`
...............................................................................................................

.. literalinclude:: /schemas/interfaces/fasta_infinitesitescharactermatrix_get.py

:meth:`StandardCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.get>`
....................................................................................................

.. literalinclude:: /schemas/interfaces/fasta_standardcharactermatrix_get.py

:meth:`DataSet.get <dendropy.datamodel.datasetmodel.DataSet.get>`
.................................................................

Note that the type of data needs to be specified using the ``data_type``
keyword argument.

.. literalinclude:: /schemas/interfaces/fasta_dataset_get.py

:meth:`DataSet.read <dendropy.datamodel.datasetmodel.DataSet.read>`
...................................................................

Note that the type of data needs to be specified using the ``data_type``
keyword argument.

.. literalinclude:: /schemas/interfaces/fasta_dataset_read.py

Writing
=======

.. _schema_specific_keyword_arguments_writing_fasta:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.fastawriter.FastaWriter.__init__

Supported Methods
-----------------

:meth:`DnaCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.write>`
..............................................................................................

.. literalinclude:: /schemas/interfaces/fasta_write.py

:meth:`RnaCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.write>`
..............................................................................................

.. literalinclude:: /schemas/interfaces/fasta_write.py

:meth:`ProteinCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.write>`
......................................................................................................

.. literalinclude:: /schemas/interfaces/fasta_write.py

:meth:`RestrictionSitesCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.write>`
........................................................................................................................

.. literalinclude:: /schemas/interfaces/fasta_write.py

:meth:`InfiniteSitesCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.write>`
..................................................................................................................

.. literalinclude:: /schemas/interfaces/fasta_write.py

:meth:`StandardCharacterMatrix.write <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.write>`
........................................................................................................

.. literalinclude:: /schemas/interfaces/fasta_write.py

:meth:`DataSet.write <dendropy.datamodel.datasetmodel.DataSet.write>`
.....................................................................

.. literalinclude:: /schemas/interfaces/fasta_write.py


