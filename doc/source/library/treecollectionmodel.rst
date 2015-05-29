*******************************************************************
:mod:`dendropy.datamodel.treecollectionmodel`: Collections of Trees
*******************************************************************

.. module:: dendropy.datamodel.treecollectionmodel

.. toctree::
    :maxdepth: 2

The |TreeList| Class
====================
.. autoclass:: dendropy.datamodel.treecollectionmodel.TreeList
    :members:
    :inherited-members:


.. method:: TreeList.read(\*\*kwargs)

.. method:: TreeList.put(\*\*kwargs)

    Write out collection of trees to file.

    :Mandatory Destimation-Specification Keyword Arguments (one and exactly one of the following required):

        - **file** (*file*) -- File or file-like object opened for writing.
        - **path** (*str*) -- Path to file to which to write.

    :Mandatory Schema-Specification Keyword Argument:

        - **schema** (*str*) -- Identifier of format of data given by the "``file``", "``path``", "``data``", or "``url``" argument specified above: ":doc:`newick </schemas/newick>`", ":doc:`nexus </schemas/nexus>`", or ":doc:`nexml </schemas/nexml>`".

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

