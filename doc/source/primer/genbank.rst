*************************************************
Working with GenBank Molecular Sequence Databases
*************************************************

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
    >>> print(char_matrix.as_string("nexus"))
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
When called without any arguments, it generates a new |TaxonNamespace| block, creating one new |Taxon| object for every sequence in the collection with a label corresponding to the identifier used to request the sequence::

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


Customizing/Controlling Sequence Taxa
-------------------------------------

The taxon assignment can be controlled in one of two ways:

    1. Using the "``label_components``" and optionally the "``label_component_separator``" arguments.
    2. Specifying a custom function using the "``gb_to_taxon_func``" argument that takes a :class:`~dendropy.interop.genbank.GenBankAccessionRecord` object and returns the |Taxon| object to be assigned to the sequence; this approach requires specification of a |TaxonNamespace| object passed using the "``taxon_namespace``" argument.

Specifying a Custom Label for Sequence Taxa
...........................................

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
    >>> print [t.label for t in char_matrix.taxon_namespace]
    ['EU105474_Homo_sapiens', 'EU105475_Homo_sapiens']
    >>> char_matrix = gb_dna.generate_char_matrix(
    ... label_components=["organism", "moltype", "gi"],
    ... label_component_separator=".")
    >>> print [t.label for t in char_matrix.taxon_namespace]
    ['Homo.sapiens.DNA.158930545', 'Homo.sapiens.DNA.158930546']

Specifying a Custom Taxon-Discovery Function
............................................

Full control over the |Taxon| object assignment process is given by using the "``gb_to_taxon_func``" argument.
This should be used to specify a function that takes a :class:`~dendropy.interop.genbank.GenBankAccessionRecord` object and returns the |Taxon| object to be assigned to the sequence.
The specification of a |TaxonNamespace| object passed using the "``taxon_namespace``" argument is also required, so that this can be assigned to the |CharacterMatrix| object.

A simple example that illustrates the usage of the "``gb_to_taxon_func``" argument by creating a custom label::

    #! /usr/bin/env python

    import dendropy
    from dendropy.interop import genbank

    def gb_to_taxon(gb):
        locality = gb.feature_table.find("source").qualifiers.find("note").value
        label = "GI" + gb.gi + "." + locality
        taxon = dendropy.Taxon(label=label)
        return taxon

    taxon_namespace = dendropy.TaxonNamespace()

    gb_dna = genbank.GenBankDna(ids=[158930545, 'EU105475'])
    char_matrix = gb_dna.generate_char_matrix(
        taxon_namespace=taxon_namespace,
        gb_to_taxon_func=gb_to_taxon)
    print [t.label for t in char_matrix.taxon_namespace]

which results in::

    ['GI158930545.Ache', 'GI158930546.Arara']

A more complex case might be where you may already have a |TaxonNamespace| with existing |Taxon| objects that you may want to associate with the sequences.
The following illustrates how to do this::


    #! /usr/bin/env python

    import dendropy
    from dendropy.interop import genbank

    tree = dendropy.Tree.get_from_string(
        "(Ache, (Arara, (Bribri, (Guatuso, Guaymi))))",
        "newick")
    def gb_to_taxon(gb):
        locality = gb.feature_table.find("source").qualifiers.find("note").value
        taxon = tree.taxon_namespace.get_taxon(label=locality)
        assert taxon is not None
        return taxon

    gb_ids = [158930545, 158930546, 158930547, 158930548, 158930549]

    gb_dna = genbank.GenBankDna(ids=gb_ids)
    char_matrix = gb_dna.generate_char_matrix(
        taxon_namespace=tree.taxon_namespace,
        gb_to_taxon_func=gb_to_taxon)
    print [t.label for t in char_matrix.taxon_namespace]
    print tree.taxon_namespace is char_matrix.taxon_namespace
    for taxon in tree.taxon_namespace:
        print "{}: {}".format(
            taxon.label,
            char_matrix[taxon].symbols_as_string()[:10])

which results in::

    True
    Ache: TCTCTTATCA
    Arara: TCTCTTATCA
    Bribri: TCTCTTATCA
    Guatuso: TCTCTTATCA
    Guaymi: TCTCTTATCA
    ['Ache', 'Arara', 'Bribri', 'Guatuso', 'Guaymi']

The important thing to note here is the the |Taxon| objects in the |DnaCharacterMatrix| do not just have the same labels as the |Taxon| object in the |Tree|, "``tree``", but actually *are* the same objects (i.e., reference the same operational taxonomic units within |DendroPy|).

Adding the GenBank Record as an Attribute
-----------------------------------------

It is sometimes useful to maintain a handle on the original GenBank record in the |CharacterMatrix| resulting from "``generate_char_matrix()``".
The "``set_taxon_attr``"  and "``set_seq_attr``" arguments of the "``generate_char_matrix()``" method allow you to this.
The values supplied to these arguments should be strings that specify the name of the attribute that will be created on the |Taxon| or |CharacterDataSequence| objects, respectively.
The value of this attribute will be the :class:`~dendropy.interop.genbank.GenBankAccessionRecord` that underlies the |Taxon| or |CharacterDataSequence| sequence.
For example::

    #! /usr/bin/env python

    import dendropy
    from dendropy.interop import genbank
    gb_dna = genbank.GenBankDna(ids=[158930545, 'EU105475'])
    char_matrix = gb_dna.generate_char_matrix(set_taxon_attr="gb_rec")
    for taxon in char_matrix.taxon_namespace:
        print "Data for taxon '{}' is based on GenBank record: {}".format(
            taxon.label,
            taxon.gb_rec.definition)

will result in::

    Data for taxon '158930545' is based on GenBank record: Homo sapiens Ache non-coding region T864 genomic sequence
    Data for taxon 'EU105475' is based on GenBank record: Homo sapiens Arara non-coding region T864 genomic sequence

Alternatively, the following::

    #! /usr/bin/env python

    import dendropy
    from dendropy.interop import genbank
    gb_dna = genbank.GenBankDna(ids=[158930545, 'EU105475'])
    char_matrix = gb_dna.generate_char_matrix(set_seq_attr="gb_rec")
    for sidx, sequence in enumerate(char_matrix.vectors()):
        print "Sequence {} ('{}') is based on GenBank record: {}".format(
            sidx+1,
            char_matrix.taxon_namespace[sidx].label,
            sequence.gb_rec.defline)

will result in::

    Sequence 1 ('158930545') is based on GenBank record: gi|158930545|gb|EU105474.1| Homo sapiens Ache non-coding region T864 genomic sequence
    Sequence 2 ('EU105475') is based on GenBank record: gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence

Annotating with GenBank Data and Metadata
-----------------------------------------

To persist the information in a the :class:`~dendropy.interop.genbank.GenBankAccessionRecord` object through serialization and deserialization, you can request that this information gets added as an  :class:`~dendropy.datamodel.basemodel.Annotation` (see ":doc:`Working with Metadata Annotations </primer/working_with_metadata_annotations>`") to the corresponding |Taxon| or |CharacterDataSequence| object.

Reference Annotation
....................

Specifying "``add_ref_annotation_to_taxa=True``" will result in a reference-style metadata annotation added to the |Taxon| object, while specifying "``add_ref_annotation_to_seqs=True``" will result in a reference-style metadata annotation added to the sequence.
The reference-style annotation is brief, single annotation that points to the URL of the original record.
As with metadata annotations in general, you really need to be using the NeXML format for full functionality.

So, for example::

    #! /usr/bin/env python

    import dendropy
    from dendropy.interop import genbank
    gb_dna = genbank.GenBankDna(ids=[158930545, 'EU105475'])
    char_matrix = gb_dna.generate_char_matrix(add_ref_annotation_to_taxa=True)
    print char_matrix.as_string("nexml")


will result in::

    <?xml version="1.0" encoding="ISO-8859-1"?>
    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009 ../xsd/nexml.xsd"
        xmlns:dcterms="http://purl.org/dc/terms/"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
        xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
    >
        <otus id="d4320533416">
            <otu id="d4323884688" label="158930545">
                <meta xsi:type="nex:ResourceMeta" rel="dcterms:source" href="http://www.ncbi.nlm.nih.gov/nucleotide/158930545" id="d4323884752" />
            </otu>
            <otu id="d4323884816" label="EU105475">
                <meta xsi:type="nex:ResourceMeta" rel="dcterms:source" href="http://www.ncbi.nlm.nih.gov/nucleotide/EU105475" id="d4323990736" />
            </otu>
        </otus>
        .
        .
        .
    </nex:nexml>

Alternatively::

    #! /usr/bin/env python

    import dendropy
    from dendropy.interop import genbank
    gb_dna = genbank.GenBankDna(ids=[158930545, 'EU105475'])
    char_matrix = gb_dna.generate_char_matrix(add_ref_annotation_to_seqs=True)
    print char_matrix.as_string("nexml")

will result in::

    <?xml version="1.0" encoding="ISO-8859-1"?>
    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009 ../xsd/nexml.xsd"
        xmlns:dcterms="http://purl.org/dc/terms/"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
        xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
    >
            <matrix>
                <row id="d4320533856" otu="d4322811536">
                    <meta xsi:type="nex:ResourceMeta" rel="dcterms:source" href="http://www.ncbi.nlm.nih.gov/nucleotide/158930545" id="d4322811600" />
                    <seq>TCTCTTATCAAAC.../seq>
                </row>
                <row id="d4320534384" otu="d4322811664">
                    <meta xsi:type="nex:ResourceMeta" rel="dcterms:source" href="http://www.ncbi.nlm.nih.gov/nucleotide/EU105475" id="d4322917584" />
                    <seq>TCTCTTATCAAAC...</seq>
                </row>
            </matrix>
        </characters>
    </nex:nexml>

Full Annotation
...............

Specifying "``add_full_annotation_to_taxa=True``" or "``add_full_annotation_to_seqs=True``" will result in the entire GenBank record being added as a set of annotations to the |Taxon| or |CharacterDataSequence| object, respectively.

For example::

    #! /usr/bin/env python

    import dendropy
    from dendropy.interop import genbank
    gb_dna = genbank.GenBankDna(ids=[158930545, 'EU105475'])
    char_matrix = gb_dna.generate_char_matrix(add_full_annotation_to_taxa=True)
    print char_matrix.as_string("nexml")

will result in the following::

    <?xml version="1.0" encoding="ISO-8859-1"?>
    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009 ../xsd/nexml.xsd"
        xmlns:genbank="http://www.ncbi.nlm.nih.gov/dtd/INSD_INSDSeq.mod.dtd"
        xmlns:dcterms="http://purl.org/dc/terms/"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
        xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
    >
        <otus id="d4320533416">
            <otu id="d4323884688" label="158930545">
                <meta xsi:type="nex:ResourceMeta" rel="dcterms:source" href="http://www.ncbi.nlm.nih.gov/nucleotide/158930545" id="d4323884752" >
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_locus" content="EU105474" id="d4323884880" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_length" content="494" id="d4323884944" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_moltype" content="DNA" id="d4323885008" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_topology" content="linear" id="d4323901520" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_strandedness" content="double" id="d4323901584" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_division" content="PRI" id="d4323901648" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_update-date" content="27-NOV-2007" id="d4323901712" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_create-date" content="27-NOV-2007" id="d4323901776" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_definition" content="Homo sapiens Ache non-coding region T864 genomic sequence" id="d4323901840" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_primary-accesison" content="EU105474" id="d4323901904" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_accession-version" content="EU105474.1" id="d4323901968" />
                    <meta xsi:type="nex:ResourceMeta" rel="genbank:otherSeqIds" id="d4323902032" >
                        <meta xsi:type="nex:LiteralMeta" property="genbank:gb" content="EU105474.1" id="d4323902160" />
                        <meta xsi:type="nex:LiteralMeta" property="genbank:gi" content="158930545" id="d4323902224" />
                    </meta>
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_source" content="Homo sapiens (human)" id="d4323902096" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_organism" content="Homo sapiens" id="d4323902288" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_taxonomy" content="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo" id="d4323902352" />
                    <meta xsi:type="nex:ResourceMeta" rel="genbank:INSDSeq_references" id="d4323902416" >
                        <meta xsi:type="nex:ResourceMeta" rel="genbank:INSDReference_reference" id="d4323902544" >
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_reference" content="1" id="d4323902672" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_position" content="1..494" id="d4323902736" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_title" content="Statistical evaluation of alternative models of human evolution" id="d4323902800" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_journal" content="Proc. Natl. Acad. Sci. U.S.A. 104 (45), 17614-17619 (2007)" id="d4323902864" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_pubmed" content="17978179" id="d4323902928" />
                        </meta>
                        <meta xsi:type="nex:ResourceMeta" rel="genbank:INSDReference_reference" id="d4323902608" >
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_reference" content="2" id="d4323903056" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_position" content="1..494" id="d4323903120" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_title" content="Direct Submission" id="d4323903184" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDReference_journal" content="Submitted (17-AUG-2007) Laboratorio de Biologia Genomica e Molecular, Pontificia Universidade Catolica do Rio Grande do Sul, Av Ipiranga 6681, Predio 12C, Sala 172, Porto Alegre, RS 90619-900, Brazil" id="d4323903248" />
                        </meta>
                    </meta>
                    <meta xsi:type="nex:ResourceMeta" rel="genbank:INSDSeq_feature-table" id="d4323902480" >
                        <meta xsi:type="nex:ResourceMeta" rel="genbank:INSDSeq_feature" id="d4323903312" >
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDFeature_key" content="source" id="d4323903440" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDFeature_location" content="1..494" id="d4323903504" />
                            <meta xsi:type="nex:ResourceMeta" rel="genbank:INSDFeature_quals" id="d4323903376" >
                                <meta xsi:type="nex:LiteralMeta" property="genbank:organism" content="Homo sapiens" id="d4323903632" />
                                <meta xsi:type="nex:LiteralMeta" property="genbank:mol_type" content="genomic DNA" id="d4323903696" />
                                <meta xsi:type="nex:LiteralMeta" property="genbank:db_xref" content="taxon:9606" id="d4323903760" />
                                <meta xsi:type="nex:LiteralMeta" property="genbank:chromosome" content="18" id="d4323903824" />
                                <meta xsi:type="nex:LiteralMeta" property="genbank:note" content="Ache" id="d4323903888" />
                            </meta>
                        </meta>
                        <meta xsi:type="nex:ResourceMeta" rel="genbank:INSDSeq_feature" id="d4323903568" >
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDFeature_key" content="misc_feature" id="d4323904016" />
                            <meta xsi:type="nex:LiteralMeta" property="genbank:INSDFeature_location" content="1..494" id="d4323904080" />
                            <meta xsi:type="nex:ResourceMeta" rel="genbank:INSDFeature_quals" id="d4323903952" >
                                <meta xsi:type="nex:LiteralMeta" property="genbank:note" content="non-coding region T864" id="d4323904208" />
                            </meta>
                        </meta>
                    </meta>
                </meta>
            </otu>
            <otu id="d4323884816" label="EU105475">
                <meta xsi:type="nex:ResourceMeta" rel="dcterms:source" href="http://www.ncbi.nlm.nih.gov/nucleotide/EU105475" id="d4324005904" >
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_locus" content="EU105475" id="d4324006032" />
                    <meta xsi:type="nex:LiteralMeta" property="genbank:INSDSeq_length" content="494" id="d4324006096" />
                    .
                    .
                    .
                    (etc.)
                </meta>
            </otu>
        </otus>
        .
        .
        .
        (etc.)
    </nex:nexml>


