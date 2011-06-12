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

NEXUS (only)
^^^^^^^^^^^^

    ``simple``
        When writing NEXUS-formatted data, if :keyword:`True`, then character data will be represented as a single "``DATA``" block, instead of separate "``TAXA``" and "``CHARACTERS``" blocks. By default this is :keyword:`False`.
    ``block_titles``
        When writing NEXUS-formatted data, if :keyword:`False`, then title statements will not be added to the various NEXUS blocks (i.e., "``TAXA``", "``CHARACTERS``", and "``TREES``"). By default, this is :keyword:`True`, i.e., block titles will be written.
    ``comment``
        When writing NEXUS-formatted data, then the contents of this variable will be added as NEXUS comment to the output. By default, this is :keyword:`None`.

.. _Customizing_Writing_NEXUS_and_Newick_Trees:

NEXUS/Newick Trees
^^^^^^^^^^^^^^^^^^

    ``suppress_rooting``
        If :keyword:`True`, will not write rooting statement. Default is :keyword:`False`.  NOTE: this keyword argument replaces the ``write_rooting`` argument which has now been deprecated.
    ``suppress_edge_lengths``
        If :keyword:`True`, will not write edge lengths. Default is :keyword:`False`.  NOTE: this keyword argument replaces the ``edge_lengths`` argument which has now been deprecated.
    ``suppress_internal_labels``
        If :keyword:`True`, internal labels will not be written. Default is :keyword:`False`.  NOTE: this keyword argument replaces the ``internal_labels`` argument which has now been deprecated.
    ``unquoted_underscores``
        If :keyword:`True`, labels with underscores will not be quoted, which will mean that they will be interpreted as spaces if read again ("soft" underscores).  If :keyword:`False`, then labels with underscores will be quoted, resulting in "hard" underscores.  Default is :keyword:`False`.  NOTE: this keyword argument replaces the ``quote_underscores`` argument which has now been deprecated.
    ``preserve_spaces``
        If :keyword:`True`, spaces not mapped to underscores in labels. Default is :keyword:`False`.
    ``store_tree_weights``
        If :keyword:`True`, tree weights are written. Default is :keyword:`False`.
    ``annotations_as_comments``
        If :keyword:`True`, will write annotations as comments. Default is :keyword:`False`.
    ``annotations_as_nhx``
        If :keyword:`True`, will write annotation as NHX statements. Default is
        :keyword:`False`.
    ``write_item_comments``
        If :keyword:`True`, will write any additional comments. Default is :keyword:`False`.

.. _Customizing_Writing_PHYLIP:

PHYLIP
^^^^^^

    ``strict``
        Write in "strict" PHYLIP format, i.e., with taxon labels truncated to 10-characters, and sequence characters beginning on column 11.

    ``spaces_to_underscores``
        Replace all spaces in taxon labels with underscores; useful if writing in relaxed mode, where spaces are used to delimit the beginning of sequence characters.
