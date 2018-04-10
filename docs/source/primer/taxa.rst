*************************************
Taxon Namespaces and Taxon Management
*************************************

Conceptual Background
=====================

Many elements of phylogenetic and, more generally, evo- or bioinformatic data
are associated with some element in the real world.
For example, a leaf :term:`node` on a :term:`tree` or a sequence in a character
matrix is typically associated with an individual biological organism, a
population or deme, a species or higher-level taxonomic group, and so on.
In classical phylogenetic literature, this referent is termed an "operational
taxonomic unit" or OTU.
In DendroPy, we use the term "taxon".
Regardless of whether the referent of the data represents an
individual organism (or even a less distinct subunit, e.g., a fragment from a
shotgun assay) or an actual taxonomic group, we apply the term "taxon".
We assign a (string) label to the concept of the entity (individual,
sub-individual, or group) represented by a "taxon", which allows us to relate
different elements of data to the same or different real-world referent.
These collections of labels representing taxon concepts are organized into
"taxon namespaces".

The concept of a "taxon namespace" is fundamental to managing data in DendroPy.
A "taxon namespace" represents a self-contained universe of *names* that
map to operational taxonomic unit *concepts*.
Operational taxonomic unit concepts are essentially names for groups of
organisms in the "real world". Operational taxonomic unit concepts are
organized into taxonomic namespaces. A taxonomic namespace is a self-contained
and functionally-complete collection of mutually-distinct operational taxonomic
unit concepts, and provide the semantic context in which operational taxonomic
units from across various data sources of different formats and provenances can
be related through correct interpretation of their taxon labels.

    * Operational taxonomic units are modeled by a |Taxon| object.

    * Taxonomic namespaces, in which operational taxonomic units are organized,
      are modeled by a |TaxonNamespace| object.

    * A |TaxonNamespace| manages a collection of |Taxon| objects, where each
      object represents a distinct operational taxonomic unit concept within
      the taxonomic namespace represented by that |TaxonNamespace| object.

    * Each |Taxon| object can belong to one and only one |TaxonNamespace|:
      |Taxon| objects are not shared across |TaxonNamespace| objects.

    * Each |Taxon| object has an attribute, ``label``, whose (string) value
      is the name of the operational taxon unit concept that it represents.

    * Different |Taxon| objects represent different operational taxonomic
      unit concepts, even if they have the same label value.

    * All client objects (`TaxonNamespaceAssociated` objects) that reference
      the same |TaxonNamespace| reference the same "universe" or domain of
      operational taxonomic unit concepts.

    * Operational taxonomic units from across different data sources are mapped
      to distinct |Taxon| objects within a particular |TaxonNamespace| based on
      matching the string values of labels of the |Taxon| object.

    * A particular taxonomic unit concept in one data source will only be
      correctly related to the same taxonomic unit concept (i.e, the same
      |Taxon| object) in another data source only if they have both
      been parsed with reference to the same taxonomic namespace (i.e., the
      same |TaxonNamespace| has been used).

    * A |TaxonNamespace| assigned an "accession index" to every |Taxon| object
      added to it. This is a stable and unique number within the context of any
      given |TaxonNamespace| object (though a |Taxon| object may have different
      accession indexes in different |TaxonNamespace| objects if it
      belongs to multiple namespaces). This number is will be used to
      calculate the "split bitmask" hash of the trivial split or external edge
      subtending the node to which this |Taxon| object is assigned on a tree.
      The concept of a "split bitmask" hash is fundamental to DendroPy's tree
      operations. The split bitmask is a hash that uniquely identifies every
      split on a tree.  It is calculated by OR'ing the split bitmask of all the
      child splits of the given split. Terminal edges, of course, do not have
      child edges, and their split bitmask is given by the accession index of
      the |Taxon| object at their head or target nodes.

Management of Shared Taxon Namespaces
=====================================

Operational taxonomic units in DendroPy are represented by |Taxon| objects, and distinct collections of operational taxonomic units are represented by |TaxonNamespace| objects.
Two distinct |Taxon| objects are considered distinct entities, *even if they share the same label*.
Understanding this is crucial to understanding management of data in DendroPy.
Many operations in DendroPy are based on the identity of the |Taxon| objects (e.g., counting of splits on trees).
Many errors by novices using DendroPy come from inadventently creating and using multiple |Taxon| objects to refer to the same taxon concept.

Every time a definition of taxa is encountered in a data source, for example, a "TAXA" block in a NEXUS file, a new |TaxonNamespace| object is created and populated with |Taxon| objects corresponding to the taxa defined in the data source.
Some data formats do not have explicit definition of taxa, e.g. a Newick tree source.
These nonetheless can be considered to have an implicit definition of a collection of operational taxonomic units given by the aggregate of all operational taxonomic units referenced in the data (i.e., the set of all distinct labels on trees in a Newick file).

Every time a reference to a taxon is encountered in a data source, such as a taxon label in a tree or matrix statement in a NEXUS file, the current |TaxonNamespace| object is searched for corresponding |Taxon| object with a matching label (see below for details on how the match is made).
If found, the |Taxon| object is used to represent the taxon.
If not, a new |Taxon| object is created, added to the |TaxonNamespace| object, and used to represent the taxon.

If multiple data sources are read, then with |TreeList| or |TreeArray| the |TaxonNamespace| instance associated with the collection through the ``taxon_namespace`` attribute will always be used to manage the |Taxon| objects, resulting in correct association of labels with |Taxon| objects across multiple reads. So, for example, the following:

.. literalinclude:: /examples/taxa_mgmt1.py

results in::

    ['A', 'B', 'C']
    ['A', 'B', 'C']

Note how the total number of taxa is three, and there is full correspondence between the taxa.
That is, the taxa referenced by "A", "B", and "C" in the second read operation were correctly mapped to the taxa from the second read operation.

With |DataSet| instances, however, each independent read operation will, by default, be managed under a *new* (i.e., independent and different) |TaxonNamespace|.

.. literalinclude:: /examples/taxa_mgmt2.py

So, if reading data from multiple data sources using a |DataSet| instance that should all be managed under the same taxon namespace, then the |TaxonNamespace| instance to use should be explicitly passed in using the "``taxon_namespace``" keyword argument:

.. literalinclude:: /examples/taxa_mgmt3.py


While each |TreeList| manages all its member |Tree| objects under the same |TaxonNamespace| reference, if two different |TreeList| instances have different |TaxonNamespace| references, then the |Taxon| objects read/managed by them *will* be necessarily different from each other, even if the labels are the same.

.. literalinclude:: /examples/taxa_mgmt4.py

Again, this can be addressed by ensuring that the |TaxonNamespace| reference is the same for |TreeList| instances that need to interact:

.. literalinclude:: /examples/taxa_mgmt5.py

The same obtains for |Tree| and |CharacterMatrix|-derived instances: if the associated |TaxonNamespace| references are different, then the associated |Taxon| objects will be different, even if the labels are the same. This will make comparison or any operation between them impossible:

.. literalinclude:: /examples/taxa_mgmt1a.py

So, if taxa are shared, then the |TaxonNamespace| to use should be passed in explicitly to ensure that each |Tree| or |CharacterMatrix|-derived instance also share the same |TaxonNamespace|:

.. literalinclude:: /examples/taxa_mgmt1b.py

Managing Taxon Name Mapping Within a Taxon Namespace
====================================================

DendroPy maps taxon definitions encountered in a data source to |Taxon| objects by the taxon label.
The labels have to match **exactly** for the taxa to be correctly mapped.
By default, this matching is case-insensitive, though case-sensitivity can be set by specifying "``case_sensitive_taxon_labels=True``".

Some quirks may arise due to some schema-specific idiosyncracies.
For example, the NEXUS standard dictates that an underscore ("_") should be substituted for a space in all labels.
Thus, when reading a NEXUS or Newick source, the taxon labels "Python_regius" and "Python regius" are exactly equivalent, and will be mapped to the same |Taxon| object.

However, this underscore-to-space mapping does **not** take place when reading, for example, a FASTA schema file.
Here, underscores are preserved, and thus "Python_regius" does not map to "Python regius".
This means that if you were to read a NEXUS file with the taxon label, "Python_regius", and later a read a FASTA file with the same taxon label, i.e., "Python_regius", these would map to different taxa!
This is illustrated by the following:

.. literalinclude:: /examples/taxon_labels1.py

Which produces the following, almost certainly incorrect, result::

    TaxonNamespace object at 0x43b4e0 (TaxonNamespace4437216): 4 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python regius'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python sebae'
        [2] Taxon object at 0x22867d0 (Taxon36202448): 'Python_regius'
        [3] Taxon object at 0x2286830 (Taxon36202544): 'Python_sebae'

Even more confusingly, if this file is written out in NEXUS schema, it would result in the space/underscore substitution taking place, resulting in two pairs of taxa with the same labels.

If you plan on mixing sources from different formats, it is important to keep in mind the space/underscore substitution that takes place by default with NEXUS/Newick formats, but does not take place with other formats.

You could simply avoid underscores and use only spaces instead:

.. literalinclude:: /examples/taxon_labels2.py

Which results in::

    TaxonNamespace object at 0x43b4e0 (TaxonNamespace4437216): 2 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python_regius'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python_sebae'

Or use underscores in the NEXUS-formatted data, but spaces in the non-NEXUS data:

.. literalinclude:: /examples/taxon_labels2b.py

Which results in the same as the preceding example::

    TaxonNamespace object at 0x43b4e0 (TaxonNamespace4437216): 2 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python regius'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python sebae'

You can also wrap the underscore-bearing labels in the NEXUS/Newick source in quotes, which preserves them from being substituted for spaces:

.. literalinclude:: /examples/taxon_labels3.py

Which will result in::

    TaxonNamespace object at 0x43c780 (TaxonNamespace4441984): 2 Taxa
        [0] Taxon object at 0x2386770 (Taxon37250928): 'Python_regius'
        [1] Taxon object at 0x2386790 (Taxon37250960): 'Python_sebae'

Finally, you can also override the default behavior of DendroPy's NEXUS/Newick parser by passing the keyword argument ``preserve_underscores=True`` to any "|read_from_methods|" or "|get_from_methods|" method. For example:

.. literalinclude:: /examples/taxon_labels4.py

will result in::

    TaxonNamespace object at 0x43c780 (TaxonNamespace4441984): 2 Taxa
        [0] Taxon object at 0x2386770 (Taxon37250928): 'Python_regius'
        [1] Taxon object at 0x2386790 (Taxon37250960): 'Python_sebae'

This may seem the simplest solution, in so far as it means that you need not maintain lexically-different taxon labels across files of different formats, but a gotcha here is that if writing to NEXUS/Newick schema, any label with underscores will be automatically quoted to preserve the underscores (again, as dictated by the NEXUS standard), which will mean that: (a) your output file will have quotes, and, as a result, (b) the underscores in the labels will be "hard" underscores if the file is read by PAUP* or DendroPy. So, for example, continuing from the previous example, the NEXUS-formatted output would look like::

    >>> print(d.as_string('nexus'))
    #NEXUS

    BEGIN TAXA;
        TITLE TaxonNamespace5736800;
        DIMENSIONS NTAX=2;
        TAXLABELS
            'Python_regius'
            'Python_sebae'
      ;
    END;

    BEGIN CHARACTERS;
        TITLE DnaCharacterMatrix37505040;
        LINK TAXA = TaxonNamespace5736800;
        DIMENSIONS  NCHAR=5;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    'Python_regius'    ACGTA
    'Python_sebae'      ACGTA
        ;
    END;

    BEGIN CHARACTERS;
        TITLE DnaCharacterMatrix37504848;
        LINK TAXA = TaxonNamespace5736800;
        DIMENSIONS  NCHAR=4;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    'Python_regius'    AAAA
    'Python_sebae'      ACGT
        ;
    END;

Note that the taxon labels have changed semantically between the input and the NEXUS output, as, according to the NEXUS standard, "Python_regius", while equivalent to "Python regius", is **not** equivalent to "'Python_regius'".
To control this, you can pass the keyword argument ``quote_underscores=False`` to any :meth:`write_to_*`, or :meth:`as_string()` method, which will omit the quotes even if the labels contain underscores::

    >>> print(d.as_string('nexus', quote_underscores=False))
    #NEXUS

    BEGIN TAXA;
        TITLE TaxonNamespace5736800;
        DIMENSIONS NTAX=2;
        TAXLABELS
            Python_regius
            Python_sebae
      ;
    END;

    BEGIN CHARACTERS;
        TITLE DnaCharacterMatrix37505040;
        LINK TAXA = TaxonNamespace5736800;
        DIMENSIONS  NCHAR=5;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    Python_regius    ACGTA
    Python_sebae      ACGTA
        ;
    END;

    BEGIN CHARACTERS;
        TITLE DnaCharacterMatrix37504848;
        LINK TAXA = TaxonNamespace5736800;
        DIMENSIONS  NCHAR=4;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    Python_regius    AAAA
    Python_sebae      ACGT
        ;
    END;
