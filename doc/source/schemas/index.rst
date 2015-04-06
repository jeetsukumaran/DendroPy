####################################################################
DendroPy Schemas: Phylogenetic and Evolutionary Biology Data Formats
####################################################################

All the "|get_from_methods|", "|read_from_methods|", and "|write_to_methods|" methods require a "``schema``" argument, which takes a :term:`schema` specification string value that describes the format (":term:`schema`) of the data.
Furthermore, these methods take other keywords arguments that vary depending on the schema to customize the way the data is parsed when read or formatted when written.
These schemas and their associated schema-specific keyword arguments are summarized in this section of the documentation.

.. toctree::
    :maxdepth: 3

    fasta
    newick
    nexml
    nexus
    phylip
