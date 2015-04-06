********************************************
Trees -- :mod:`dendropy.datamodel.treemodel`
********************************************

.. module:: dendropy.datamodel.treemodel

.. toctree::
    :maxdepth: 2

The :class:`Tree` Class
=======================
.. autoclass:: dendropy.datamodel.treemodel.Tree
    :members:
    :inherited-members:
    :exclude-members: get_from_stream

.. classmethod:: Tree.get_from_stream(src, schema, label=None, taxon_namespace=None, collection_offset=None, tree_offset=None, \*\*kwargs,)

    Instantiate and return a *single* |Tree| object from a file opened for
    reading given by ``src``.

    :Parameters:    - **src** (*file* or *str*) -- File-like object (``get_from_stream``), path to file (``get_from_path``), or string (``get_from_string``) with :term:`tree` data.
                    - **schema** (*str*) -- Identifier of format of data given in ``src``.
                    - **label** (*str*) -- Name or identifier to be assigned to the new |Tree|; if not given, will be assigned the one specified in the data source, or `None` otherwise.
                    - **taxon_namespace** (|TaxonNamespace|) -- The |TaxonNamespace| instance to use to :doc:`manage the taxon names </primer/taxa>`. If not specified, a new one will be created.
                    - **collection_offset** (*int*) -- 0-based index of tree block or collection in source to be parsed. If ``tree_offset`` is specified, then ``collection_offset`` must be also be specified. If neither is specified, then the first tree given in the source will be selected.
                    - **tree_offset** (*int*) -- 0-based index of tree within the collection specified by ``collection_offset`` to be parsed to be parsed. If ``tree_offset`` is specified, then ``collection_offset`` must be also be specified. If neither is specified, then the first tree given in the source will be selected.
                    - **\*\*kwargs** -- schema-specific keyword arguments. These provide control over how the data is interpreted and processed, depending on the schemas used.

The :class:`Node` Class
=======================
.. autoclass:: dendropy.datamodel.treemodel.Node
    :members:
    :inherited-members:

The :class:`Edge` Class
=======================
.. autoclass:: dendropy.datamodel.treemodel.Edge
    :members:
    :inherited-members:

The :class:`Bipartition` Class
==============================
.. autoclass:: dendropy.datamodel.treemodel.Bipartition
    :members:
    :inherited-members:

