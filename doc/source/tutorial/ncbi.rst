**************************************************************
NCBI (National Center for Biotechnology Information) Databases
**************************************************************

The :mod:`~dendropy.interop.ncbi` module provides the :class:`~dendropy.interop.ncbi.Entrez` class, which wraps up some basic querying and retrieval of data from the NCBI (National Center for Biotechnology Information) life-sciences databases.

Retrieving Nucleotide Data
==========================

By Lists of Accession Identifiers
---------------------------------

The :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accessions` method takes a list of nucleotide accession numbers as arguments::

    >>> from dendropy.interop import ncbi
    >>> entrez = ncbi.Entrez()
    >>> data = entrez.fetch_nucleotide_accessions(['EU105474', 'EU105475'])
    >>> for t in data.taxon_set:
    ...     print(t)
    ...
    gi|158930545|gb|EU105474.1| Homo sapiens Ache non-coding region T864 genomic sequence
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence

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
    gi|158930545|gb|EU105474.1| Homo sapiens Ache non-coding region T864 genomic sequence
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence
    gi|158930547|gb|EU105476.1| Homo sapiens Bribri non-coding region T864 genomic sequence
    gi|158930548|gb|EU105477.1| Homo sapiens Guatuso non-coding region T864 genomic sequence
    .
    .
    .
    >>> data.write_to_path('data2.fas', 'fasta')

RNA vs. DNA
-----------

As noted, it is assumed that the data is DNA, and thus the query will result in a :class:`~dendropy.dataobject.char.DnaCharacterMatrix` object. By specifying ``matrix_type=RnaCharacterMatrix`` to :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accessions` or :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accession_range`, you can retrieve RNA data as well.

Error Handling
--------------

By default, if you were to give non-existing accession numbers, an exception will be thrown::

    >>> entrez = ncbi.Entrez()
    >>> data = entrez.fetch_nucleotide_accessions(['zzz0', 'zzz1'])
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "dendropy/interop/ncbi.py", line 476, in fetch_genbank_records
          raise Entrez.AccessionFetchError(missing)
      dendropy.interop.ncbi.AccessionFetchError: Failed to retrieve accessions: zzz0, zzz1

An exception will be thrown even if some of the specified accessions are valid:::

    >>> data = entrez.fetch_nucleotide_accessions(['zzz0', 'zzz1', 'EU105475'])
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "dendropy/interop/ncbi.py", line 476, in fetch_genbank_records
        raise Entrez.AccessionFetchError(missing)
    dendropy.interop.ncbi.AccessionFetchError: Failed to retrieve accessions: zzz0, zzz1

By passing in ``verify=False`` to :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accessions` or :meth:`~dendropy.interop.ncbi.Entrez.fetch_nucleotide_accession_range`, you can request that data retrieval failures can be ignore, and only existing accessions be returned::

    >>> data = entrez.fetch_nucleotide_accessions(['zzz0', 'zzz1', 'EU105475'], verify=False)
    >>> len(data)
    1
    >>> for t in data.taxon_set:
    ...     print(t.label)
    ...
    gi|158930546|gb|EU105475.1| Homo sapiens Arara non-coding region T864 genomic sequence


Note that specifying ``verify=False`` means that you might end up with empty :class:`~dendropy.dataobject.char.DnaCharacterMatrix`  objects::

    >>> data = entrez.fetch_nucleotide_accessions(['zzz0', 'zzz1'], verify=False)
    >>> len(data)
    0

Also, perhaps more of a concern, turning off verification may lead to **wrong** sequences being retrieved.
For example, when trying to download a range of accessions, but inadvertently omitting to specify a ``prefix`` value to be pre-pended to identifiers might result in matching the wrong sequences, based on GI values::

    >>> data = entrez.fetch_nucleotide_accession_range(1000, 1001, verify=False)
    >>> print(len(data))
    2
    >>> for t in data.taxon_set:
    ...     print t.label
    ...
    gi|1000|gb|Z18639.1| D.leucas gene for large subunit rRNA
    gi|1001|gb|Z18638.1| D.leucas gene for small subunit rRNA

Here, the sequences were retrieved based on matching GI numbers (1000, 1001) rather than the accession ids (e.g., "AY1000", "AY1001").

Switching off verification can also lead to some confusing errors. For example::

    >>> data = entrez.fetch_nucleotide_accession_range(1000, 1003, verify=False)
    ---
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File ":dendropy/interop/ncbi.py", line 282, in fetch_nucleotide_accession_range

      File "dendropy/interop/ncbi.py", line 245, in fetch_nucleotide_accessions
        sys.stderr.write("---\nNCBI Entrez Query returned:\n%s\n---\n" % results_str)
      File "dendropy/utility/iosys.py", line 199, in get_from_string
        readable.read_from_string(src, schema, **kwargs)
      File "dendropy/utility/iosys.py", line 260, in read_from_string
        return self.read(stream=s, schema=schema, **kwargs)
      File "dendropy/dataobject/char.py", line 653, in read
        d = DataSet(stream=stream, schema=schema, **kwargs)
      File "dendropy/dataobject/dataset.py", line 90, in __init__
        self.process_source_kwargs(**kwargs)
      File "dendropy/utility/iosys.py", line 221, in process_source_kwargs
        self.read(stream=stream, schema=schema, **kwargs)
      File "dendropy/dataobject/dataset.py", line 172, in read
        raise x
    dendropy.utility.error.DataParseError: Error parsing data source on line 8 at column 3: Unrecognized sequence symbol "P"

Here, the sequence with GI number of "1003" was a protein sequence, so it included characters not part of the DNA alphabet, resulting in the :class:`~dendropy.utility.error.DataParseError` exception being raised.

(Auto-)Generating Analysis-Friendly Sequence Labels
===================================================

When fetching nucleotides, you can request new taxon labels are generated by specifying ``generate_new_labels=True``.
By default, this results in taxon labels consisting of the GenBank record accession and the organism name, separated by a space::

    >>> entrez = ncbi.Entrez()
    >>> data = entrez.fetch_nucleotide_accessions(['EU105474', 'EU105475', 'EU105476'],
    ...         generate_new_labels=True)
    ...
    >>> for taxon in data:
    ...     print taxon.label
    ...
    >>> EU105474 Homo sapiens
    >>> EU105475 Homo sapiens
    >>> EU105476 Homo sapiens


You can customize the label composition formula by passing in a list of field names of the GenBank record to be used using the ``label_components`` argument::

    >>> data = entrez.fetch_nucleotide_accessions(['EU105474', 'EU105475', 'EU105476'],
    ...         generate_new_labels=True,
    ...         label_components=["organism", "accession"])
    ...
    >>> for taxon in data:
    ...     print taxon.label
    ...
    Homo sapiens EU105474
    Homo sapiens EU105475
    Homo sapiens EU105476

The token separating the label components, which defaults to a space, can be
changed using the ``label_component_separator`` argument::

    >>> data = entrez.fetch_nucleotide_accessions(['EU105474', 'EU105475', 'EU105476'],
    ...         generate_new_labels=True,
    ...         label_components=["organism", "moltype", "gi" ],
    ...         label_component_separator=".")
    ...
    >>> for taxon in data:
    ...     print taxon.label
    ...
    Homo.sapiens.DNA.158930545
    Homo.sapiens.DNA.158930546
    Homo.sapiens.DNA.158930547

