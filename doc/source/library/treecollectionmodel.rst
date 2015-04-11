
*********************************************************************
Collections of Trees -- :mod:`dendropy.datamodel.treecollectionmodel`
*********************************************************************

.. module:: dendropy.datamodel.treecollectionmodel

.. toctree::
    :maxdepth: 2

The |TreeList| Class
====================
.. autoclass:: dendropy.datamodel.treecollectionmodel.TreeList
    :members:
    :exclude-members: get,put

.. classmethod:: TreeList.get(\*\*kwargs)

    Instantiate and return a *new* |TreeList| object from a data source providing one or more collections of trees.

    :Mandatory Source-Specification Keyword Arguments (one and exactly one of the following required):

        - **file** (*file*) -- File or file-like object with data opened for reading.
        - **path** (*str*) -- Path to file with data.
        - **url** (*str*) -- URL providing data.
        - **value** (*str*) -- Data given directly.

    :Mandatory Schema-Specification Keyword Argument:

        - **schema** (*str*) -- Identifier of format of data given by the "``file``", "``path``", "``value``", or "``url``" argument specified above: ":doc:`newick </schemas/newick>`", ":doc:`nexus </schemas/nexus>`", or ":doc:`nexml </schemas/nexml>`".

    :Optional General Keyword Arguments:

        - **label** (*str*) -- Name or identifier to be assigned to the new |TreeList|; if not given, will be assigned the one specified in the data source, or `None` otherwise.
        - **taxon_namespace** (|TaxonNamespace|) -- The |TaxonNamespace| instance to use to :doc:`manage the taxon names </primer/taxa>`. If not specified, a new one will be created.
        - **collection_offset** (*int*) -- 0-based index of tree block or collection in source to be parsed. If not specified then the first collection (offset = 0) is assumed.
        - **tree_offset** (*int*) -- 0-based index of first tree within the collection specified by ``collection_offset`` to be parsed (i.e., skipping the first ``tree_offset`` trees). If not specified, then the first tree (offset = 0) is assumed.

    :Optional Schema-Specific Keyword Arguments:

        -  These provide control over how the data is interpreted and processed, and supported argument names and values depend on the schema as specified by the value passed as the "``schema``" argument:
            -   :ref:`"newick" (Newick) <schema_specific_keyword_arguments_reading_newick>`
            -   :ref:`"nexus" (Nexus) <schema_specific_keyword_arguments_reading_nexus>`
            -   :ref:`"nexml" (NeXML) <schema_specific_keyword_arguments_reading_nexml>`

.. method:: TreeList.read(\*\*kwargs)

    Add |Tree| objects to existing |TreeList| from data source providing one or more collections of trees.

    :Mandatory Source-Specification Keyword Arguments (one and exactly one of the following required):

        - **file** (*file*) -- File or file-like object with data opened for reading.
        - **path** (*str*) -- Path to file with data.
        - **url** (*str*) -- URL providing data.
        - **value** (*str*) -- Data given directly.

    :Mandatory Schema-Specification Keyword Argument:

        - **schema** (*str*) -- Identifier of format of data given by the "``file``", "``path``", "``value``", or "``url``" argument specified above: ":doc:`newick </schemas/newick>`", ":doc:`nexus </schemas/nexus>`", or ":doc:`nexml </schemas/nexml>`".

    :Optional General Keyword Arguments:

        - **label** (*str*) -- Name or identifier to be assigned to the new |TreeList|; if not given, will be assigned the one specified in the data source, or `None` otherwise.
        - **taxon_namespace** (|TaxonNamespace|) -- The |TaxonNamespace| instance to use to :doc:`manage the taxon names </primer/taxa>`. If not specified, a new one will be created.
        - **collection_offset** (*int*) -- 0-based index of tree block or collection in source to be parsed. If not specified then the first collection (offset = 0) is assumed.
        - **tree_offset** (*int*) -- 0-based index of first tree within the collection specified by ``collection_offset`` to be parsed (i.e., skipping the first ``tree_offset`` trees). If not specified, then the first tree (offset = 0) is assumed.

    :Optional Schema-Specific Keyword Arguments:

        -  These provide control over how the data is interpreted and processed, and supported argument names and values depend on the schema as specified by the value passed as the "``schema``" argument:
            -   :ref:`"newick" (Newick) <schema_specific_keyword_arguments_reading_newick>`
            -   :ref:`"nexus" (Nexus) <schema_specific_keyword_arguments_reading_nexus>`
            -   :ref:`"nexml" (NeXML) <schema_specific_keyword_arguments_reading_nexml>`

.. method:: TreeList.put(\*\*kwargs)

    Write out collection of trees to file.

    :Mandatory Destimation-Specification Keyword Arguments (one and exactly one of the following required):

        - **file** (*file*) -- File or file-like object opened for writing.
        - **path** (*str*) -- Path to file to which to write.

    :Mandatory Schema-Specification Keyword Argument:

        - **schema** (*str*) -- Identifier of format of data given by the "``file``", "``path``", "``value``", or "``url``" argument specified above: ":doc:`newick </schemas/newick>`", ":doc:`nexus </schemas/nexus>`", or ":doc:`nexml </schemas/nexml>`".

    :Optional Schema-Specific Keyword Arguments:

        -  These provide control over how the data is formatted, and supported argument names and values depend on the schema as specified by the value passed as the "``schema``" argument:
            -   :ref:`"newick" (Newick) <schema_specific_keyword_arguments_writing_newick>`
            -   :ref:`"nexus" (Nexus) <schema_specific_keyword_arguments_writing_nexus>`
            -   :ref:`"nexml" (NeXML) <schema_specific_keyword_arguments_writing_nexml>`


The |TreeArray| Class
=====================
.. autoclass:: dendropy.datamodel.treecollectionmodel.TreeArray
    :members:

The |SplitDistribution| Class
=============================
.. autoclass:: dendropy.datamodel.treecollectionmodel.SplitDistribution
    :members:

The |SplitDistributionSummarizer| Class
=======================================
.. autoclass:: dendropy.datamodel.treecollectionmodel.SplitDistributionSummarizer
    :members:

