**********************************
GenBank Genetic Sequence Databases
**********************************

The :mod:`~dendropy.interop.genbank` module provides the classes and methods to download sequences from |GenBank| and instantiate them into |DendroPy| phylogenetic data objects.
Three classes are provided, all of which have an identical interface, varying only in the type of data retrieved:

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
identifiers to be retrieved using the "``ids``"  argument.
The value of this argument should be a container with either GenBank accession identifiers or GI numbers::


    >>> from dendropy.interop import genbank
    >>> gb_dna = genbank.GenBankDna(ids=['EU105474', 'EU105475'])
    >>> for gb in gb_dna:
    ...     print gb
    gi|158930545|gb|EU105474.1| Homo sapiens Ache non-coding region T864 genomic sequence
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence

The records are stored as :class:`~dendropy.interop.genbank.GenBankAccessionRecord` objects.
These records store the *full* information available in a |GenBank| record, including the references, feature table, qualifiers, and other details, and these are available as attributes of the :class:`~dendropy.interop.genbank.GenBankAccessionRecord` objects (e.g., "``primary_accession``", "``taxonomy``", "``feature_table``" and so on).

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

As can be seen, by default the taxon labels assigned to the sequences are set to the identifier used to request the sequences. This, and many other aspects of the character matrix generation, including annotation of taxa and sequences, can be customized, as discussed in detail below.

Acquiring Data from GeneBank
============================

The :class:`~dendropy.interop.genbank.GenBankDna`, :class:`~dendropy.interop.genbank.GenBankRna`, and :class:`~dendropy.interop.genbank.GenBankProtein` classes provide for the downloading and management of DNA, RNA, and protein (AA) sequences from |GenBank|.
The first two of these query the "nucleotide" or "nuccore" database, while the last queries the "protein" database.
The constructors of these classes accept the following arguments:

    ``ids``

        A list of accession identifiers of GI numbers of the records to be downloaded. E.g. "``ids=['EU105474', 'EU105475']``",  "``ids=['158930545', 'EU105475']``", or  "``ids=['158930545', '158930546']``".
        If "``prefix``" is specified, this string will be pre-pended to all values in the list.

    ``id_range``
        A tuple of *integers* that specify the first and last values (inclusive) of accession or GI numbers of the records to be downloaded. If "``prefix``" is specified, this string will be prepended to all numbers in this range.
        Thus specifying "``id_range=(158930545, 158930550)``" is exactly equivalent to specifying "``ids=[158930545, 158930546, 158930547, 158930548, 158930549, 158930550]``", while specifying "``id_range=(105474, 105479), prefix="EU"``" is exactly equivalent tp specifying "``ids=["EU105474", "EU105475", "EU105476", "EU105477", "EU105478", "EU105479"]``".


    ``prefix``
        This string will be prepended to all values resulting from the "``ids``" and "``id_range``".


    ``verify``
        By default, the results of the download are checked to make sure there is a one-to-one correspondence between requested id's and retrieved records. Setting "``verify=False``" skips this checking.

So, for example, the following are all different ways of instantiating |GenBank| resource data store::

    >>> from dendropy.interop import genbank
    >>> gb_dna = genbank.GenBankDna(ids=['EU105474', 'EU105475'])
    >>> gb_dna = genbank.GenBankDna(ids=['158930545', 'EU105475'])
    >>> gb_dna = genbank.GenBankDna(ids=['158930545', '158930546'])
    >>> gb_dna = genbank.GenBankDna(ids=['105474', '105475'], prefix="EU")
    >>> gb_dna = genbank.GenBankDna(id_range=(105474, 105478), prefix="EU")
    >>> gb_dna = genbank.GenBankDna(id_range=(158930545, 158930546))

You can add more records to an existing instance of :class:`~dendropy.interop.genbank.GenBankDna`, :class:`~dendropy.interop.genbank.GenBankRna`, or :class:`~dendropy.interop.genbank.GenBankProtein` objects by using the "``acquire``" or "``acquire_range``" methods.
The "``acquire``" method takes a sequence of accession identifiers or GI numbers for the first argument ("``ids``"), and, in addition an optional string prefix to be prepended can be supplied using the second argument, "``prefix``", while verification can be disabled by specifying |False| for the third argument, "``verify``".
The "``acquire_range``" method takes two mandatory *integer* arguments: the first and last value of the range of accession or GI numbers of the records to be downloaded.
As with the other method, a string prefix to be prepended can be optionally supplied using the argument "``prefix``", while verification can be disabled by specifying "``verify=|False|``".
For example::


    >>> from dendropy.interop import genbank
    >>> gb_dna = genbank.GenBankDna(['EU105474', 'EU105475'])
    >>> print len(gb_dna)
    >>> gb_dna.acquire([158930547, 158930548])
    >>> print len(gb_dna)
    >>> gb_dna.acquire_range(105479, 105480, prefix="EU")
    >>> print len(gb_dna)
    2
    4
    6

Accessing GenBank Records
=========================

The |GenBank| records accumulated in :class:`~dendropy.interop.genbank.GenBankDna`, :class:`~dendropy.interop.genbank.GenBankRna`, and :class:`~dendropy.interop.genbank.GenBankProtein` objects are represented by collections of :class:`~dendropy.interop.genbank.GenBankAccessionRecord` objects.
Each of these :class:`~dendropy.interop.genbank.GenBankAccessionRecord` objects represent the full information from the |GenBank| source as a rich Python object.

    >>> from dendropy.interop import genbank
    >>> gb_dna = genbank.GenBankDna(['EU105474', 'EU105475'])
    >>> for gb_rec in gb_dna:
    ...    print gb_rec.gi
    ...    print gb_rec.locus
    ...    print gb_rec.length
    ...    print gb_rec.moltype
    ...    print gb_rec.topology
    ...    print gb_rec.strandedness
    ...    print gb_rec.division
    ...    print gb_rec.update_date
    ...    print gb_rec.create_date
    ...    print gb_rec.definition
    ...    print gb_rec.primary_accession
    ...    print gb_rec.accession_version
    ...    print "(other seq ids)"
    ...    for osi_key, osi_value in gb_rec.other_seq_ids.items():
    ...        print "    ", osi_key, osi_value
    ...    print gb_rec.source
    ...    print gb_rec.organism
    ...    print gb_rec.taxonomy
    ...    print "(references)"
    ...    for ref in gb_rec.references:
    ...        print "    ", ref.number , ref.position , ref.authors , ref.consrtm , ref.title , ref.journal , ref.medline_id , ref.pubmed_id , ref.remark
    ...    print "(feature_table)"
    ...    for feature in gb_rec.feature_table:
    ...        print "    ", feature.key, feature.location
    ...        for qualifier in feature.qualifiers:
    ...            print "        ", qualifier.name, qualifier.value
    ...
    158930545
    EU105474
    494
    DNA
    linear
    double
    PRI
    27-NOV-2007
    27-NOV-2007
    Homo sapiens Ache non-coding region T864 genomic sequence
    EU105474
    EU105474.1
    (other seq ids)
        gb EU105474.1
        gi 158930545
    Homo sapiens (human)
    Homo sapiens
    Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Eutel...
    (references)
        1 1..494 [] None Statistical evaluation of alternativ...
        2 1..494 [] None Direct Submission Submitted (17-AUG-...
    (feature_table)
        source 1..494
            organism Homo sapiens
            mol_type genomic DNA
            db_xref taxon:9606
            chromosome 18
            note Ache
        misc_feature 1..494
            note non-coding region T864
    .
    .
    .
    (etc.)

Generating Character Matrix Objects from GenBank Data
=====================================================

The "``generate_char_matrix()``" method of :class:`~dendropy.interop.genbank.GenBankDna`, :class:`~dendropy.interop.genbank.GenBankRna`, and :class:`~dendropy.interop.genbank.GenBankProtein` objects creates and returns a |CharacterMatrix| object of the appropriate type out of the data collected in them.
When called without any arguments, it generates a new |TaxonSet| block, creating one new |Taxon| object for every sequence in the collection with a label corresponding to the identifier used to request the sequence::

    >>> from dendropy.interop import genbank
    >>> gb_dna = genbank.GenBankDna(ids=[158930545, 'EU105475'])
    >>> char_matrix = gb_dna.generate_char_matrix()
    >>> print char_matrix.as_string("nexus")
    #NEXUS

    BEGIN TAXA;

        DIMENSIONS NTAX=2;
        TAXLABELS
            158930545
            EU105475
    ;
    END;

    BEGIN CHARACTERS;
        DIMENSIONS  NCHAR=494;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    158930545    TCTCTTATCAAACTA...
    EU105475     TCTCTTATCAAACTA...
        ;
    END;


    BEGIN SETS;
    END;


Customizing/Controlling Sequences Taxa
--------------------------------------

The taxon assignment can be controlled in one of two ways:

    1. Using the "``label_components``" and optionally the "``label_component_separator``" arguments.
    2. Specifying a custom function using the "``gb_to_taxon_func``" argument that takes a :class:`~dendropy.interop.genbank.GenBankAccessionRecord` object and returns the |Taxon| object to be assigned to the sequence; this approach requires specification of a |TaxonSet| object passed using the "``taxon_set``" argument.

Specifying a Custom Label for Sequence Taxat
............................................

The "``label_components``" and the "``label_component_separator``" arguments allow for customization of the taxon labels of the |Taxon| objects created for each sequence.
The "``label_components``" argument should be assigned an ordered container (e.g., a list) of strings that correspond to attributes of objects of the :class:`~dendropy.interop.genbank.GenBankAccessionRecord` class.
The values of these attributes will be concatenated to compose the |Taxon| object label.
By default, the components will be separated by spaces, but you can override this by passing the string to be used by the "``label_component_separator``" argument.
For example::


    >>> from dendropy.interop import genbank
    >>> gb_dna = genbank.GenBankDna(ids=[158930545, 'EU105475'])
    >>> char_matrix = gb_dna.generate_char_matrix(
    ... label_components=["accession", "organism", ],
    ... label_component_separator="_")
    >>> print [t.label for t in char_matrix.taxon_set]
    ['EU105474_Homo_sapiens', 'EU105475_Homo_sapiens']
    >>> char_matrix = gb_dna.generate_char_matrix(
    ... label_components=["organism", "moltype", "gi"],
    ... label_component_separator=".")
    >>> print [t.label for t in char_matrix.taxon_set]
    ['Homo.sapiens.DNA.158930545', 'Homo.sapiens.DNA.158930546']

Specifying a Custom Taxon-Discovery Function
............................................

Full control over the |Taxon| object assignment process is given by using the "``gb_to_taxon_func``" argument.
This should be used to specify a function that takes a :class:`~dendropy.interop.genbank.GenBankAccessionRecord` object and returns the |Taxon| object to be assigned to the sequence.
The specification of a |TaxonSet| object passed using the "``taxon_set``" argument is also required, so that this can be assigned to the |CharacterMatrix| object.
