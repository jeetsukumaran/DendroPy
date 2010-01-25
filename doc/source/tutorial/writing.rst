*************************
Writing Phylogenetic Data
*************************

The |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes all support the following instance methods for writing data:

    - :meth:`write_to_stream(dest, schema, **kwargs)`
        Takes a file or file-like object opened for writing the data as the first argument, and a string specifying the schema as the second.

    - :meth:`write_to_path(dest, schema, **kwargs)`
        Takes a string specifying the path to the file as the first argument, and a string specifying the schema as the second.

    - :meth:`as_string(schema, **kwargs)`
        Takes a string specifying the schema as the first argument, and returns a string containing the formatted-representation of the data.

As above, the schema specification can be any supported and type-apppropriate schema, such as "nexus", "newick", "nexml", "dnafasta", "rnafasta", "proteinfasta" etc.

