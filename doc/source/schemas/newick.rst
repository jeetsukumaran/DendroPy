******
NEWICK
******
Reading
=======

|Tree|
------

Create and return a new |Tree| instance from a |Newick|-formatted data source
.............................................................................

    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_stream`
    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_path`
    -   :meth:`~dendropy.datamodel.treemodel.Tree.get_from_string`

Populate an existing |Tree| instance with data from a |Newick|-formatted source
...............................................................................

Not supported.

Usage Examples
..............

.. literalinclude:: /schemas/examples/newick_read_tree.py


Details
=======

.. autodocstringonly:: dendropy.dataio.newickreader.NewickReader.__init__
