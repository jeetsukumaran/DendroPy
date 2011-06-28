*************************
Taxa and Taxon Management
*************************

Operational taxonomic units in DendroPy are represented by |Taxon| objects, and distinct collections of operational taxonomic units are represented by |TaxonSet| objects.

Every time a definition of taxa is encountered in a data source, for example, a "TAXA" block in a NEXUS file, a new |TaxonSet| object is created and populated with |Taxon| objects corresponding to the taxa defined in the data source.
Some data formats do not have explicit definition of taxa, e.g. a Newick tree source.
These nonetheless can be considered to have an implicit definition of a collection of operational taxonomic units given by the aggregate of all operational taxonomic units referenced in the data (i.e., the set of all distinct labels on trees in a Newick file).

Every time a reference to a taxon is encountered in a data source, such as a taxon label in a tree or matrix statement in a NEXUS file, the current |TaxonSet| object is searched for corresponding |Taxon| object with a matching label (see below for details on how the match is made). If found, the |Taxon| object is used to represent the taxon. If not, a new |Taxon| object is created, added to the |TaxonSet| object, and used to represent the taxon.

DendroPy maps taxon definitions encountered in a data source to |Taxon| objects by the taxon label.
The labels have to match **exactly** for the taxa to be correctly mapped

Some special formats, such as NEXUS or Newick, treat the taxa labels as **case-insensitive**: "Python regius", "PYTHON REGIUS" and "python regius" will all be considered the same taxon (this can be turned off by specifying the appropriate |False| to ``case_insensitive_taxon_labels`` when reading the data in NEXUS or Newick formats).
Otherwise, in general, most other formats (e.g., PHYLIP, Fasta, NExML) treat the taxa labels as **case-sensitive**: "Python regius", "PYTHON REGIUS" and "python regius" will all be considered the different taxa.

Further quirks may arise due to some schema-specific idiosyncracies.
For example, the NEXUS standard dictates that an underscore ("_") should be substituted for a space in all labels.
Thus, when reading a NEXUS or Newick source, the taxon labels "Python_regius" and "Python regius" are exactly equivalent, and will be mapped to the same |Taxon| object.

However, this underscore-to-space mapping does **not** take place when reading, for example, a FASTA schema file.
Here, underscores are preserved, and thus "Python_regius" does not map to "Python regius".
This means that if you were to read a NEXUS file with the taxon label, "Python_regius", and later a read a FASTA file with the same taxon label, i.e., "Python_regius", these would map to different taxa!
This is illustrated by the following:

.. literalinclude:: /examples/taxon_labels1.py
    :linenos:

Which produces the following, almost certainly incorrect, result::

    TaxonSet object at 0x43b4e0 (TaxonSet4437216): 4 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python regius'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python sebae'
        [2] Taxon object at 0x22867d0 (Taxon36202448): 'Python_regius'
        [3] Taxon object at 0x2286830 (Taxon36202544): 'Python_sebae'

Even more confusingly, if this file is written out in NEXUS schema, it would result in the space/underscore substitution taking place, resulting in two pairs of taxa with the same labels.

If you plan on mixing sources from different formats, it is important to keep in mind the space/underscore substitution that takes place by default with NEXUS/Newick formats, but does not take place with other formats.

You could simply avoid underscores and use only spaces instead:

.. literalinclude:: /examples/taxon_labels2.py
    :linenos:

Which results in::

    TaxonSet object at 0x43b4e0 (TaxonSet4437216): 2 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python_regius'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python_sebae'

Or use underscores in the NEXUS-formatted data, but spaces in the non-NEXUS data:

.. literalinclude:: /examples/taxon_labels2b.py
    :linenos:

Which results in the same as the preceding example::

    TaxonSet object at 0x43b4e0 (TaxonSet4437216): 2 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python regius'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python sebae'

You can also wrap the underscore-bearing labels in the NEXUS/Newick source in quotes, which preserves them from being substituted for spaces:

.. literalinclude:: /examples/taxon_labels3.py
    :linenos:

Which will result in::

    TaxonSet object at 0x43c780 (TaxonSet4441984): 2 Taxa
        [0] Taxon object at 0x2386770 (Taxon37250928): 'Python_regius'
        [1] Taxon object at 0x2386790 (Taxon37250960): 'Python_sebae'

Finally, you can also override the default behavior of DendroPy's NEXUS/Newick parser by passing the keyword argument ``preserve_underscores=True`` to any :meth:`read_from_*`, :meth:`get_from_*` or stream-parsing constructor. For example:

.. literalinclude:: /examples/taxon_labels4.py
    :linenos:

will result in::

    TaxonSet object at 0x43c780 (TaxonSet4441984): 2 Taxa
        [0] Taxon object at 0x2386770 (Taxon37250928): 'Python_regius'
        [1] Taxon object at 0x2386790 (Taxon37250960): 'Python_sebae'

This may seem the simplest solution, in so far as it means that you need not maintain lexically-different taxon labels across files of different formats, but a gotcha here is that if writing to NEXUS/Newick schema, any label with underscores will be automatically quoted to preserve the underscores (again, as dictated by the NEXUS standard), which will mean that: (a) your output file will have quotes, and, as a result, (b) the underscores in the labels will be "hard" underscores if the file is read by PAUP* or DendroPy. So, for example, continuing from the previous example, the NEXUS-formatted output would look like::

    >>> print(d.as_string('nexus'))
    #NEXUS

    BEGIN TAXA;
        TITLE TaxonSet5736800;
        DIMENSIONS NTAX=2;
        TAXLABELS
            'Python_regius'
            'Python_sebae'
      ;
    END;

    BEGIN CHARACTERS;
        TITLE DnaCharacterMatrix37505040;
        LINK TAXA = TaxonSet5736800;
        DIMENSIONS  NCHAR=5;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    'Python_regius'    ACGTA
    'Python_sebae'      ACGTA
        ;
    END;

    BEGIN CHARACTERS;
        TITLE DnaCharacterMatrix37504848;
        LINK TAXA = TaxonSet5736800;
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
        TITLE TaxonSet5736800;
        DIMENSIONS NTAX=2;
        TAXLABELS
            Python_regius
            Python_sebae
      ;
    END;

    BEGIN CHARACTERS;
        TITLE DnaCharacterMatrix37505040;
        LINK TAXA = TaxonSet5736800;
        DIMENSIONS  NCHAR=5;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    Python_regius    ACGTA
    Python_sebae      ACGTA
        ;
    END;

    BEGIN CHARACTERS;
        TITLE DnaCharacterMatrix37504848;
        LINK TAXA = TaxonSet5736800;
        DIMENSIONS  NCHAR=4;
        FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;
        MATRIX
    Python_regius    AAAA
    Python_sebae      ACGT
        ;
    END;
