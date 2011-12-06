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
        When writing a |DataSet| object, if passed a specific |TaxonSet|, then **only** |TreeList| and |CharacterMatrix| objects associated with this |TaxonSet| will be written. By default, this is |None|, meaning that all data in the |DataSet| object will be written.

    ``exclude_trees``
        When writing a |DataSet| object, if |True|, then **no** tree data will be written (i.e., all |TreeList| objects in the |DataSet| will be skipped in the output). By default, this is |False|, meaning that all tree data will be written.

    ``exclude_chars``
        When writing a |DataSet| object, if |True|, then **no** characer data will be written (i.e., all |CharacterMatrix| objects in the |DataSet| will be skipped in the output). By default, this is |False|, meaning that all character data will be written.

.. _Customizing_Writing_NEXUS_and_Newick:

NEXUS and Newick
^^^^^^^^^^^^^^^^

The following code fragment shows a typical invocation of a NEXUS-format write operation using all supported keywords with their defaults:

.. literalinclude:: /examples/skel_dataio_nexus_writing.py

The following code fragment shows a typical invocation of a Newick-format write operation using all supported keyword arguments with their default values:

.. literalinclude:: /examples/skel_dataio_newick_writing.py

.. NEXUS and Newick share mostly the same format for writing tree statements. As such, in DendroPy the same set of keyword arguments can be used to control and customize both NEXUS and Newick output (though the defaults for a few of these keywords vary between formats). In addtion, because it is more extensive than Newick, several other keyword arguments are supported when writing in NEXUS format.

The special keywords supported for writing NEXUS-formatted output include:

    ``simple``
        When writing NEXUS-formatted data, if |True|, then character data will be represented as a single "``DATA``" block, instead of separate "``TAXA``" and "``CHARACTERS``" blocks. By default this is |False|.
    ``block_titles``
        When writing NEXUS-formatted data, if |False|, then title statements will not be added to the various NEXUS blocks (i.e., "``TAXA``", "``CHARACTERS``", and "``TREES``"). By default, this is |True|, i.e., block titles will be written.
    ``suppress_taxa_block``
        If |True|, do not write a "TAXA" block. Default is |False|.
    ``exclude_trees``
        When writing NEXUS-formatted data, if |True|, then **no** tree data will be written (i.e., all |TreeList| objects in the |DataSet| will be skipped in the output). By default, this is |False|, meaning that all tree data will be written.
    ``exclude_chars``
        When writing NEXUS-formatted data, if |True|, then **no** characer data will be written (i.e., all |CharacterMatrix| objects in the |DataSet| will be skipped in the output). By default, this is |False|, meaning that all character data will be written.
    ``preamble_blocks``
        When writing NEXUS-formatted data, a list of other blocks (or strings) to be written at the beginning of the file.
    ``supplemental_blocks``
        When writing NEXUS-formatted data, a list of other blocks (or strings) to be written at the end of the file.
    ``file_comments``
        When writing NEXUS-formatted data, then the contents of this variable (a string or a list of strings) will be added as a NEXUS comment to the file (at the top). By default, this is |None|.

The special keywords supported for writing both NEXUS- or Newick-formatted trees include:

    ``suppress_leaf_taxon_labels``
        If |True|, then taxon labels will not be printed for leaves.  Default is |False|.
    ``suppress_leaf_node_labels``
        If |False|, then node labels (if available) will be printed for leaves. Defaults to |True|. Note that DendroPy distinguishes between *taxon* labels and *node* labels. In a typical NEWICK string, taxon labels are printed for leaf nodes, while leaf node labels are ignored (hence the default '|True|' setting, to ignore leaf *node* labels).
    ``suppress_internal_taxon_labels``
        If |True|, then taxon labels will not be printed for internal nodes.  Default is |False|.  NOTE: this replaces the ``internal_labels`` argument which has been deprecated.
    ``suppress_internal_node_labels``
        If |True|, internal node labels will not be written. Default is |False|.  NOTE: this replaces the ``internal_labels`` argument which has been deprecated.
    ``suppress_rooting``
        If |True|, will not write rooting statement. Default is |False|.  NOTE: this keyword argument replaces the ``write_rooting`` argument which has now been deprecated.
    ``suppress_edge_lengths``
        If |True|, will not write edge lengths. Default is |False|.  NOTE: this keyword argument replaces the ``edge_lengths`` argument which has now been deprecated.
    ``unquoted_underscores``
        If |True|, labels with underscores will not be quoted, which will mean that they will be interpreted as spaces if read again ("soft" underscores).  If |False|, then labels with underscores will be quoted, resulting in "hard" underscores.  Default is |False|.  NOTE: this keyword argument replaces the ``quote_underscores`` argument which has now been deprecated.
    ``preserve_spaces``
        If |True|, spaces not mapped to underscores in labels. Default is |False|.
    ``store_tree_weights``
        If |True|, tree weights are written. Default is |False|.
    ``suppress_annotations``
        If |True|, will **not** write annotated attributes as comments. Default is |False| if writing in NEXUS format *and* ``simple`` is |False|; otherwise, if writing in NEWICK format or NEXUS format with ``simple`` set to |True|, then defaults to |True|.
    ``annotations_as_nhx``
        If |True| and ``suppress_annotations`` is |True|, then annotations will be written in NHX format ('[&&field=value:field=value]'), as opposed to a more generic 'hot comment' format with only one leading ampersand ('[&field=value,field=value,field={value,value}]'). Defaults to |False|.
    ``suppress_item_comments``
        If |True|, will **not** write any additional comments associated with (tree) items. Default is |False| if writing in NEXUS format *and* ``simple`` is |False|; otherwise, if writing in NEWICK format or NEXUS format with ``simple`` set to |True|, then defaults to |True|.
    ``node_label_element_separator``
        If both ``suppress_leaf_taxon_labels`` and ``suppress_leaf_node_labels`` are |False|, then this will be the string used to join them. Defaults to ' '.
    ``node_label_compose_func``
        If not None, should be a function that takes a |Node| object as an argument and returns the string to be used to represent the node in the tree statement. The return value from this function is used unconditionally to print a node representation in a tree statement, by-passing the default labelling function (and thus ignoring ``suppress_leaf_taxon_labels``, ``suppress_leaf_node_labels=True``, ``suppress_internal_taxon_labels``, ``suppress_internal_node_labels``, etc.). Defaults to |None|.
    ``edge_label_compose_func``
        If not None, should be a function that takes an |Edge| object as
        an argument, and returns the string to be used to represent the
        edge length in the tree statement.

.. _Customizing_Writing_FASTA:

FASTA
^^^^^

The following code fragment shows a typical invocation of a FASTA-format write operation using all supported keywords with their defaults:

.. literalinclude:: /examples/skel_dataio_fasta_writing.py

The special keywords supported for writing FASTA-formatted data include:

    ``wrap``
        If |True|, then sequences will be wrapped at ``wrap_width`` characters. Defaults to |False|. Output is prettier, but writing operations are considerably slower.

    ``wrap_width``
        If ``wrap`` is |True|, then sequences will be wrapped at these many characters. Defaults to 70.


.. _Customizing_Writing_PHYLIP:

PHYLIP
^^^^^^

The following code fragment shows a typical invocation of a PHYLIP-format write operation using all supported keywords with their defaults:

.. literalinclude:: /examples/skel_dataio_phylip_writing.py

The special keywords supported for writing PHYLIP-formatted data include:

    ``strict``
        If |True|, write in "strict" PHYLIP format, i.e., with taxon labels truncated to 10-characters, and sequence characters beginning on column 11. Defaults to |False|: writes in "relaxed" format (taxon labels not truncated, and separated from sequence characters by more two consecutive spaces).

    ``spaces_to_underscores``
        If |True|, replace all spaces in taxon labels with underscores; useful if writing in relaxed mode, where spaces are used to delimit the beginning of sequence characters. Defaults to |False|: labels not changed.

    ``force_unique_taxon_labels``
        If |True|, then identical taxon labels (or labels that are identical due to truncation) will be disambiguated through the appending of indexes.
