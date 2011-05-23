*************************
Writing Phylogenetic Data
*************************

Writing to Streams, Filepaths, or Strings
=========================================

The |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes all support the following instance methods for writing data:

    :meth:`write_to_stream(dest, schema, **kwargs)`
        Takes a file or file-like object opened for writing the data as the first argument, and a string specifying the schema as the second.

    :meth:`write_to_path(dest, schema, **kwargs)`
        Takes a string specifying the path to the file as the first argument, and a string specifying the schema as the second.

    :meth:`as_string(schema, **kwargs)`
        Takes a string specifying the schema as the first argument, and returns a string containing the formatted-representation of the data.

.. _Specifying_the_Data_Writing_Format:

Specifying the Data Writing Format
==================================

The schema specification string can be one of the following:

    "``nexus``"
        To write |Tree|, |TreeList|, |CharacterMatrix|, or |DataSet| objects in NEXUS format.

    "``newick``"
        To write |Tree|, |TreeList|, or |DataSet| objects in Newick format. With |DataSet| objects, only tree data will be written.

    "``fasta``"
        To write |CharacterMatrix| or |DataSet| objects in FASTA format. With |DataSet| objects, only character data will be written.

    "``phylip``"
        To write |CharacterMatrix| or |DataSet| objects in PHYLIP format. With |DataSet| objects, only character data will be written.

.. _Customizing_the_Data_Writing_Format:

Customizing the Data Writing Format
===================================

The writing of data can be controlled or fine-tuned using keyword arguments. As with reading, some of these arguments apply generally, while others are only available or make sense for a particular format.

.. _Customizing_Writing_All_Formats:

All Formats
^^^^^^^^^^^

    ``taxon_set``
        When writing a |DataSet| object, if passed a specific |TaxonSet|, then **only** |TreeList| and |CharacterMatrix| objects associated with this |TaxonSet| will be written. By default, this is :keyword:`None`, meaning that all data in the |DataSet| object will be written.

    ``exclude_trees``
        When writing a |DataSet| object, if :keyword:`True`, then **no** tree data will be written (i.e., all |TreeList| objects in the |DataSet| will be skipped in the output). By default, this is :keyword:`False`, meaning that all tree data will be written.

    ``exclude_chars``
        When writing a |DataSet| object, if :keyword:`True`, then **no** characer data will be written (i.e., all |CharacterMatrix| objects in the |DataSet| will be skipped in the output). By default, this is :keyword:`False`, meaning that all character data will be written.

.. _Customizing_Writing_NEXUS_and_Newick:

NEXUS/Newick
^^^^^^^^^^^^

    ``simple``
        When writing NEXUS-formatted data, if :keyword:`True`, then character data will be represented as a single "``DATA``" block, instead of separate "``TAXA``" and "``CHARACTERS``" blocks. By default this is :keyword:`False`.

    ``write_rooting``
        If :keyword:`False`, then tree rooting statements (e.g., "``[&R]``" or "``[&U]``") will not be prefixed to the tree statements. By default, this is :keyword:`True`, i.e., rooting statements will be written.

    ``edge_lengths``
        If :keyword:`False`, then edge or branch lengths will not be written as part of the tree statements. By default, this is :keyword:`True`, i.e., edge lengths will be written.

    ``internal_labels``
        If :keyword:`False`, then labels for internal nodes (if given) will not be written as part of the tree statements. By default, this is :keyword:`True`, i.e., internal node labels will be written.

    ``write_item_comments``
        If :keyword:`True`, then comments associated with nodes on trees will be written. When writing NEXUS formats, this defaults to :keyword:`True` for NEXUS formats, *unless* ``simple=True`` (see above) is specified, in which case it defaults to :keyword:`False` unless explicitly overridden by calling code. When writing NEWICK formats, this defaults to :keyword:`False`.

    ``block_titles``
        When writing NEXUS-formatted data, if :keyword:`False`, then title statements will not be added to the various NEXUS blocks (i.e., "``TAXA``", "``CHARACTERS``", and "``TREES``"). By default, this is :keyword:`True`, i.e., block titles will be written.

    ``preserve_spaces``
        If :keyword:`True`, then no attempt will be made to produce unquoted labels by substituting spaces for underscores. By default, this is :keyword:`False`, i.e., any label that includes spaces but no other special punctuation character or underscores will have all spaces replaced by underscores so as to allow the label to be represented without quotes.

    ``quote_underscores``
        If :keyword:`False`, then labels will not be wrapped in quotes even if they contain underscores (meaning that the underscores will be interpreted as spaces according to the NEXUS standard). By default, this is :keyword:`True`, meaning that any label that contains underscores will be wrapped in quotes. Note that if a label has any other characters requiring quote protection as specified by the NEXUS standard, then the label will be quoted regardless of the value of this keyword argument.

    ``comment``
        When writing NEXUS-formatted data, then the contents of this variable will be added as NEXUS comment to the output. By default, this is :keyword:`None`.

.. _Customizing_Writing_PHYLIP:

PHYLIP
^^^^^^

    ``strict``
        Write in "strict" PHYLIP format, i.e., with taxon labels truncated to 10-characters, and sequence characters beginning on column 11.

    ``spaces_to_underscores``
        Replace all spaces in taxon labels with underscores; useful if writing in relaxed mode, where spaces are used to delimit the beginning of sequence characters.
