####################################################################
DendroPy Schemas: Phylogenetic and Evolutionary Biology Data Formats
####################################################################

.. _Specifying_the_Data_Source_Format:
.. _Specifying_the_Data_Writing_Format:

In DendroPy, the format of the data is called its "schema".
All the data import and export methods require specification of the data format through a "``schema``" keyword argument, which takes a *schema specification string* as a value.
This is a string identifer that uniquely maps to a particular format, and should be one of the following values:

    - ":doc:`fasta </schemas/fasta>`"
    - ":doc:`newick </schemas/newick>`"
    - ":doc:`nexus </schemas/nexus>`"
    - ":doc:`nexml </schemas/nexml>`"
    - ":doc:`phylip </schemas/phylip>`"


.. _Schema_Specific_Keyword_Arguments:
.. _Customizing_the_Data_Writing_Format:

Furthermore, the various data reading and writing methods take other keywords arguments that vary depending on the schema to customize the way the data is parsed when read or formatted when written.
These *schema-specific keyword arguments* are detailed in this section under each schema description, along with information and code templates/examples.


Contents
========

.. toctree::
    :maxdepth: 3

    fasta
    newick
    nexml
    nexus
    phylip
