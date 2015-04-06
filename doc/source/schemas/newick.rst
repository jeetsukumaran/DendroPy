******
NEWICK
******

Reading
=======

|Tree|
------

Create and Return a New |Tree| Instance
.......................................

Supported Methods
^^^^^^^^^^^^^^^^^

    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_string`


Representative Examples
^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: /schemas/interfaces/newick_tree_get_from_path.py

Add Data to an Existing |Tree| Instance
.......................................

Not supported.

|TreeList|
----------

Create and Return a New |TreeList| Instance from a |Newick|-formatted Data Source
.................................................................................

    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_string`

Populate an Existing |TreeList| Instance with Data from a |Newick|-Formatted Source
...................................................................................

    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_string`

|TreeArray|
-----------

Create and Return a New |TreeArray| Instance from a |Newick|-formatted Data Source
..................................................................................

Not supported.


Populate an Existing |TreeArray| Instance with Data from a |Newick|-Formatted Source
....................................................................................

    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_string`

.. _schema_specific_keyword_arguments_reading_newick:

Schema-Specific Keyword Arguments
=================================

.. autokeywordargumentsonly:: dendropy.dataio.newickreader.NewickReader.__init__
