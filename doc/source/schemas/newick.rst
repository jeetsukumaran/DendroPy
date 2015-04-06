******
NEWICK
******

.. contents::
    :local:
    :backlinks: none

Reading
=======

.. schema_specific_keyword_arguments_reading_newick:

Supported Methods
-----------------

|Tree|
......

:meth:`~dendropy.datamodel.treemodel.Tree.get`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: /schemas/interfaces/newick_tree_get_from_stream.py


|TreeList|
..........

    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_stream`
            .. literalinclude:: /schemas/interfaces/newick_treelist_get_from_path.py
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_path`
            .. literalinclude:: /schemas/interfaces/newick_treelist_get_from_path.py
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.get_from_string`
            .. literalinclude:: /schemas/interfaces/newick_treelist_get_from_path.py
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_stream`
            .. literalinclude:: /schemas/interfaces/newick_treelist_get_from_path.py
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_path`
            .. literalinclude:: /schemas/interfaces/newick_treelist_get_from_path.py
    -   :meth:`~dendropy.datamodel.treemodel.TreeList.read_from_string`
            .. literalinclude:: /schemas/interfaces/newick_treelist_get_from_path.py

|TreeArray|
...........

    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.TreeArray.read_from_string`

.. literalinclude:: /schemas/interfaces/newick_treearray_read_from_path.py

|DataSet|
.........

    -   :meth:`~dendropy.datamodel.treemodel.DataSet.read_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.DataSet.read_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.DataSet.read_from_string`

.. literalinclude:: /schemas/interfaces/newick_dataset_get_from_path.py

Schema-Specific Keyword Arguments
---------------------------------

.. autokeywordargumentsonly:: dendropy.dataio.newickreader.NewickReader.__init__
