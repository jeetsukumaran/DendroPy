******
NEWICK
******

.. contents::
    :local:
    :backlinks: none

Reading
=======

Examples
--------

.. literalinclude:: /schemas/interfaces/newick_tree_get_from_path.py

.. schema_specific_keyword_arguments_reading_newick:

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.newickreader.NewickReader.__init__

Supported Methods
-----------------

|Tree|
......

    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_string`

|TreeList|
..........

    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_string`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_string`

|TreeArray|
...........

    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_string`

|DataSet|
.........

    -   :meth:`~dendropy.datamodel.treemodel.DataSet.read_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.DataSet.read_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.DataSet.read_from_string`


