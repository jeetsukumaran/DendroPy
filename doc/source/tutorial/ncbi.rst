**************************************************************
NCBI (National Center for Biotechnology Information) Databases
**************************************************************

The :mod:`~dendropy.interop.ncbi` module provides the :class:`~dendropy.interop.ncbi.Entrez` class, which wraps up some basic querying and retrieval of data from the NCBI (National Center for Biotechnology Information) life-sciences databases.

Retrieving Nucleotide Data
==========================

By Lists of Accession Identifiers
---------------------------------

The :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accession_ids` method takes a list of nucleotide accession numbers as arguments::

    >>> from dendropy.interop import ncbi
    >>> entrez = ncbi.Entrez()
    >>> data = entrez.fetch_nucleotide_accession_ids(['EU105474', 'EU105475'])
    >>> for t in data.taxon_set:
    ...     print(t)
    ...
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence
    gi|158930545|gb|EU105474.1| Homo sapiens Ache non-coding region T864 genomic sequence

As can be seen, If successful, it returns a :class:`~dendropy.dataobject.char.DnaCharacterMatrix` object that consists of sequences corresponding to the requested accessions. This can, of course, be written to a file or manipulated as needed::

    >>> data.write_to_path('eu105474-105475.nex', 'nexus')

By Ranges of Accession Identifiers
----------------------------------

Sometimes, it might be more convenient to specify the required accession numbers as a range. For example, a publication may list sequences accessioned as "EU10574-106045". The :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accession_range` takes three arguments: a numeric value indicating the start of the range, a numeric value indicating the end of the range, and string giving a prefix to be added to each number within the range to yield the full accession identifier.
Note that, unlike Python's native ``range`` function, the last or end value *is included* as part of the range.
So, to get the all the sequences given in a publication as "EU10574-106045"::

    >>> from dendropy.interop import ncbi
    >>> entrez = ncbi.Entrez()
    >>> data = entrez.fetch_nucleotide_accession_range(105474, 106045, prefix="EU")
    >>> for t in data.taxon_set:
    ...     print(t)
    ...
    gi|158930636|gb|EU105565.1| Homo sapiens Arara non-coding region T946 genomic sequence
    gi|158930635|gb|EU105564.1| Homo sapiens Ache non-coding region T946 genomic sequence
    gi|158930638|gb|EU105567.1| Homo sapiens Guatuso non-coding region T946 genomic sequenc
    .
    .
    .
    >>> data.write_to_path('data2.fas', 'fasta')

Error Handling
--------------

By default, if you were to give non-existing accession numbers, an exception will be thrown::

    >>> entrez = ncbi.Entrez()
    >>> data = entrez.fetch_nucleotide_accession_ids(['zzz0', 'zzz1'])
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "dendropy/interop/ncbi.py", line 232, in fetch_nucleotide_accession_ids
        raise Entrez.AccessionFetchError(missing_ids)
    dendropy.interop.ncbi.AccessionFetchError: Failed to retrieve accessions: zzz0, zzz1

By default, an exception will be thrown even if some of the specified accessions are valid:::

    >>> data = entrez.fetch_nucleotide_accession_ids(['zzz0', 'zzz1', 'EU105475'])
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "dendropy/interop/ncbi.py", line 232, in fetch_nucleotide_accession_ids
        raise Entrez.AccessionFetchError(missing_ids)
    dendropy.interop.ncbi.AccessionFetchError: Failed to retrieve accessions: zzz0, zzz1

By passing in ``verify=False`` to :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accession_ids` or :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accession_range`, you can request that data retrieval failures can be ignore, and only existing accessions be returned::

    >>> data = entrez.fetch_nucleotide_accession_ids(['zzz0', 'zzz1', 'EU105475'], verify=False)
    >>> len(data)
    1
    >>> for t in data.taxon_set:
    ...     print(t.label)
    ...
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence


Specifying ``verify=False`` means that you might end up with empty :class:`~dendropy.dataobject.char.DnaCharacterMatrix`  objects::

    >>> data = entrez.fetch_nucleotide_accession_ids(['zzz0', 'zzz1'], verify=False)
    >>> len(data)
    0

(Auto-)Generating Analysis-Friendly Sequence Labels
===================================================

When fetching nucleotides, you can request the :class:`~dendropy.interop.ncbi.Entrez` object to generate labels that are little more compact and analysis friendly by passing ``generate_label=True`` to the constructor. This will generate a new taxon label for sequence based on the GenBank FASTA defline value. By default, it will compose a label in the form of:

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
    >>> data.write_to_path('gb2.nex', 'nexus')

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
    >>> data.write_to_path('seqs.dat', 'phylip', strict=False)

