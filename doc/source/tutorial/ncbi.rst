**************************************************************
NCBI (National Center for Biotechnology Information) Databases
**************************************************************

The :mod:`~dendropy.interop.ncbi` module provides the :class:`~dendropy.interop.ncbi.Entrez` class, which wraps up some basic querying and retrieval of data from the NCBI (National Center for Biotechnology Information) life-sciences databases.

At the moment, the only functionality supported is for fetching data from the nucleotide database by accession numbers, which will be returned in the form of a :class:`~dendropy.dataobject.char.DnaCharacterMatrix` object. For example::

    >>> from dendropy.interop import ncbi
    >>> entrez = ncbi.Entrez()
    >>> data = entrez.fetch_nucleotide_accession_ids(['EU10574', 'EU10575'])
    >>> for t in data.taxon_set:
    ...     print(t)
    ...
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence
    gi|158930545|gb|EU105474.1| Homo sapiens Ache non-coding region T864 genomic sequence

Sometimes, it might be more convenient to specify the required accession numbers as a range::

    >>> from dendropy.interop import ncbi
    >>> entrez = ncbi.Entrez()
    >>> data = entrez.fetch_nucleotide_accession_range(10574, 10577, prefix="EU")
    >>> for t in data.taxon_set:
    ...     print(t)
    ...
    gi|158930547|gb|EU105476.1| Homo sapiens Bribri non-coding region T864 genomic sequence
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence
    gi|158930545|gb|EU105474.1| Homo sapiens Ache non-coding region T864 genomic sequence

Note that, consistent with Python's native ``range`` function, the last or end value is **not** part of the range. So, the expanded list of accession numbers specified by the example immediately above is ``['EU10574', 'EU10575', 'EU10576']``; ``EU10577`` is not included.

You can request the :class:`~dendropy.interop.ncbi.Entrez` object to generate labels that are little more compact and analysis friendly by passing ``generate_label=True`` to the constructor. This will generate a new taxon label for sequence based on the GenBank FASTA defline value. By default, it will compose a label in the form of:

    <GBNUM>_<Genus>_<species>_<other>

So, for example, a sequence with the defline:

    gi|158930547|gb|EU105476.1| Homo sapiens Bribri non-coding region T864 genomic sequence

will get the taxon label:

    EU105476_Homo_sapiens_Bribri

You can control details of the label construction by the following arguments to the constructor:

    - ``label_num_desc_components`` specifies the number of components from the defline to use. By default, this is 3, which usually corresponds (in a sensible defline) to the genus name, the species epithet, and either the sub-species or locality information.
    - ``label_separator`` specifies the string used in between different label components. By default, this is an underscore.
    - ``label_gbnum_in_front`` specifies whether the GenBank accession number should form the beginning
        (``True``; default) or tail (``False``) end of the label.

Furthermore, you can request that the data get sorted by label value by specifying ``sort_taxa_by_label=True``.

So, for example::

    >>> entrez = ncbi.Entrez(generate_labels=True, sort_taxa_by_label=True)
    >>> data = entrez.fetch_nucleotide_accession_ids(['EU105474', 'EU105475', 'EU105476'])
    >>> for t in data.taxon_set:
    ...     print(t)
    ...
    EU105474_Homo_sapiens_Ache
    EU105475_Homo_sapiens_Arara
    EU105476_Homo_sapiens_Bribri

Or::

    >>> entrez = ncbi.Entrez(generate_labels=True,
    ...         label_num_desc_components=2,
    ...         label_gbnum_in_front=False,
    ...         label_separator='.')
    >>> data = entrez.fetch_nucleotide_accession_ids(['EU105474', 'EU105475', 'EU105476'])
    >>> for t in data.taxon_set:
    ...     print(t)
    ...
    Homo.sapiens.EU105476
    Homo.sapiens.EU105475
    Homo.sapiens.EU105474

