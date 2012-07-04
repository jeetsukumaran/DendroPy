**********************************
GenBank Genetic Sequence Databases
**********************************

The :mod:`~dendropy.interop.genbank` module provides the classes and methods to download sequences from |GenBank| and instantiate them into |DendroPy| phylogenetic data objects.
Three classes are provided, all of which have an identical interface, varying only in that type of data retrieved:

   :class:`~dendropy.interop.genbank.GenBankDna`

        Acquire and manage DNA sequence data from the |GenBank| Nucleotide database.

   :class:`~dendropy.interop.genbank.GenBankRna`

        Acquire and manage RNA sequence data from the |GenBank| Nucleotide database.

   :class:`~dendropy.interop.genbank.GenBankProtein`

        Acquire and manage AA sequence data from the |GenBank| Protein database.


Quick Start
===========

The basic way to retrieve sequence data is create a
:class:`~dendropy.interop.genbank.GenBankDna`,
:class:`~dendropy.interop.genbank.GenBankRna`, or
:class:`~dendropy.interop.genbank.GenBankProtein` object, and pass in a list of
identifiers to be retrieved using the "``ids``"  argument::


    >>> from dendropy.interop import genbank
    >>> gb_dna = genbank.GenBankDna(ids=['EU105474', 'EU105475'])
    >>> for gb in gb_dna:
    ...     print gb
    gi|158930545|gb|EU105474.1| Homo sapiens Ache non-coding region T864 genomic sequence
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence

The records are stored as :class:`~dendropy.interop.genbank.GenBankAccessionRecord` objects.
To generate a |CharacterMatrix| object from the collection of sequences, call the :meth:`~dendropy.interop.genbank.GenBankDna.generate_char_matrix`  method::

    >>> from dendropy.interop import genbank
    >>> gb_dna = genbank.GenBankDna(ids=['EU105474', 'EU105475'])
    >>> char_matrix = gb_dna.generate_char_matrix()
    >>> print char_matrix.as_string("nexus")
    #NEXUS
    BEGIN TAXA;

        DIMENSIONS NTAX=2;
        TAXLABELS
            EU105474
            EU105475
    ;
    END;
    BEGIN CHARACTERS;
        DIMENSIONS  NCHAR=494;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    EU105474    TCTCTTATCA...
    EU105475    TCTCTTATCA...
    ;
    END;


