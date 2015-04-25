*************************
Writing Phylogenetic Data
*************************

Writing to Files, Filepaths, or Strings
=======================================

The |Tree|, |TreeList|, |CharacterMatrix|-derived (i.e., |DnaCharacterMatrix|,
|ProteinCharacterMatrix|, |StandardCharacterMatrix|, etc.), and |DataSet|
classes all support a "|write|" instance method for serialization of data to an
external data source.
This method takes two mandatory keyword arguments:

    -   One and exactly one of the following to specify the *destination*:
        -   a path to a file (specified using the keyword argument "``path``")
        -   a file or a file-like object opened for writing (specified using the keyword argument ``"file"``)

    -   A ":ref:`schema specification string <Specifying_the_Data_Source_Format>`" given by the keyword argument "``schema``", to identify the schema or format for the output.

Alternatively, the |Tree|, |TreeList|, |CharacterMatrix|-derived, or |DnaCharacterMatrix| objects may also be represented as a string by calling the "``as_string()``" method, which requires at least one single mandatory argument, "``schema``", giving the ":ref:`schema specification string <Specifying_the_Data_Source_Format>`" to identify the format of the output.

In either case, the ":ref:`schema specification string <Specifying_the_Data_Source_Format>`" can be one of: ":doc:`fasta </schemas/fasta>`", ":doc:`newick </schemas/newick>`", ":doc:`nexus </schemas/nexus>`", ":doc:`nexml </schemas/nexml>`", or ":doc:`phylip </schemas/phylip>`".

For example:

.. code-block:: python


    tree.write(path="output.tre", schema="newick")
    dest = open("output.xml", "w")
    tree_list.write(file=dest, schema="nexml")
    print(dna_character_matrix.as_string(schema="fasta"))


As with the "|get|" and "|read|" methods, further keyword arguments can be specified to control behavior.
These are covered in detail in the ":doc:`/schemas/index`" section.

.. note::

    The |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes also support a "|write_to_methods|" family of instance methods that can be seen as specializations of the "|write|" method for various types of destinations:

        :meth:`write_to_stream(dest, schema, \*\*kwargs)`
            Takes a file or file-like object opened for writing the data as the first argument, and a string specifying the schema as the second.

        :meth:`write_to_path(dest, schema, \*\*kwargs)`
            Takes a string specifying the path to the file as the first argument, and a string specifying the schema as the second.

        :meth:`as_string(schema, \*\*kwargs)`
            Takes a string specifying the schema as the first argument, and returns a string containing the formatted-representation of the data.

