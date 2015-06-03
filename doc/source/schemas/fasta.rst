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

``DnaCharacterMatrix.get``
..........................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/fasta_dnacharactermatrix_get.py

``RnaCharacterMatrix.get``
..........................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/fasta_rnacharactermatrix_get.py

``ProteinCharacterMatrix.get``
..............................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/fasta_proteincharactermatrix_get.py

``RestrictionSitesCharacterMatrix.get``
.......................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/fasta_restrictionsitescharactermatrix_get.py

``InfiniteSitesCharacterMatrix.get``
....................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/fasta_infinitesitescharactermatrix_get.py

``StandardCharacterMatrix.get``
...............................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.get>`)

.. literalinclude:: /schemas/interfaces/fasta_standardcharactermatrix_get.py

``DataSet.get``
...............
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.get>`)

Note that the type of data needs to be specified using the ``data_type``
keyword argument.

.. literalinclude:: /schemas/interfaces/fasta_dataset_get.py

``DataSet.read``
................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.read>`)

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

``DnaCharacterMatrix.write``
.............................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/fasta_write.py

``DnaCharacterMatrix.as_string``
................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/fasta_as_string.py


``RnaCharacterMatrix.write``
.............................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/fasta_write.py

``RnaCharacterMatrix.as_string``
................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/fasta_as_string.py


``ProteinCharacterMatrix.write``
.................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/fasta_write.py

``ProteinCharacterMatrix.as_string``
.....................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/fasta_as_string.py


``RestrictionSitesCharacterMatrix.write``
.........................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/fasta_write.py

``RestrictionSitesCharacterMatrix.as_string``
.............................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/fasta_as_string.py


``InfiniteSitesCharacterMatrix.write``
......................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/fasta_write.py

``InfiniteSitesCharacterMatrix.as_string``
..........................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/fasta_as_string.py


``StandardCharacterMatrix.write``
.................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.write>`)

.. literalinclude:: /schemas/interfaces/fasta_write.py

``StandardCharacterMatrix.as_string``
.....................................
(:meth:`method reference <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.as_string>`)

.. literalinclude:: /schemas/interfaces/fasta_as_string.py


``DataSet.write``
.................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.write>`)

.. literalinclude:: /schemas/interfaces/fasta_write.py

``DataSet.as_string``
.....................
(:meth:`method reference <dendropy.datamodel.datasetmodel.DataSet.as_string>`)

.. literalinclude:: /schemas/interfaces/fasta_as_string.py

