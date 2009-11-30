*******************************************
Working with Trees and Collections of Trees
*******************************************

Trees in are represented by the class :class:`~dendropy.dataobject.tree.Tree`.
Every :class:`~dendropy.dataobject.tree.Tree` object has a :attr:`~dendropy.dataobject.tree.Tree.seed_node` attribute. If the tree is rooted, then this is the root node. If the tree is not rooted, however, then this is an artificial node that serves as the "starting point" for the tree.
The :attr:`~dendropy.dataobject.tree.Tree.seed_node`, like every other node on the tree, is a :class:`~dendropy.dataobject.tree.Node` object.
Every :class:`~dendropy.dataobject.tree.Node` object maintains a list of its immediate child :class:`~dendropy.dataobject.tree.Node` objects as well as a reference to its parent :class:`~dendropy.dataobject.tree.Node` object.
You can request a shallow-copy :func:`~list` of child :class:`~dendropy.dataobject.tree.Node` objects using the :meth:`~dendropy.dataobject.tree.Node.child_nodes()` method, and you can access the parent :class:`~dendropy.dataobject.tree.Node` object directly through the :attr:`~dendropy.dataobject.tree.Node.parent_node` attribute.
By definition, the :attr:`~dendropy.dataobject.tree.Tree.seed_node` has no parent node, leaf nodes have no child nodes, and internal nodes have both parent nodes and child nodes.

Customizing Tree Creation and Reading
=====================================

Interpreting Rootings
---------------------

The rooting state of a :class:`~dendropy.dataobject.tree.Tree` object is set by the :attr:`~dendropy.dataobject.tree.Tree.is_rooted` property.
When parsing NEXUS- and Newick-formatted data, the rooting states of the resulting :class:`~dendropy.dataobject.tree.Tree` objects are given by ``[&R]`` (for rooted) or ``[&U]`` (for unrooted) comment tags preceding the tree definition in the data source.
If these tags are not present, then the trees are assumed to be unrooted.
This behavior can be changed by specifying keyword arguments to the :meth:`get_from_*()`,  or :meth:`read_from_*()` methods of both the :class:`~dendropy.dataobject.tree.Tree` and :class:`~dendropy.dataobject.tree.TreeList` classes, or the constructors of these classes when specifying a data source from which to construct the tree:

The ``as_rooted`` keyword argument, if :keyword:`True`, forces all trees to be interpreted as rooted, regardless of whether or not the ``[&R]``/``[&U]`` comment tags are given.
Conversely, if :keyword:`False`, all trees will be interpreted as unrooted.
For semantic clarity, you can also specify ``as_unrooted`` to be :keyword:`True` to force all trees to be unrooted.

.. literalinclude:: /examples/tree_rootings1.py
    :linenos:

In addition, you can specify a ``default_as_rooted`` keyword argument, which, if :keyword:`True`, forces all trees to be interpreted as rooted, *if* the ``[&R]``/``[&U]`` comment tags are *not* given.
Otherwise the rooting will follow the ``[&R]``/``[&U]`` commands.
Conversely, if ``default_as_rooted`` is :keyword:`False`, all trees will be interpreted as unrooted if the ``[&R]``/``[&U]`` comment tags are not given.
Again, for semantic clarity, you can also specify ``default_as_unrooted`` to be :keyword:`True` to assume all trees are unrooted if not explicitly specified, though, as this is default behavior, this should not be neccessary.


Tree Traversal
==============

Iterating Over Nodes
--------------------

The following example shows how you might evolve a continuous character on a tree by recursively visting each node, and setting the value of the character to one drawn from a normal distribution centered on the value of the character of the node's ancestor and standard deviation given by the length of the edge subtending the node:

.. literalinclude:: /examples/tree_evolve_char1.py
    :linenos:

While the previous example works, it is probably clearer and more efficient to use one of the pre-defined node iterator methods:

    :meth:`~dendropy.dataobject.tree.Tree.preorder_node_iter()`
        Iterates over nodes in a :class:`~dendropy.dataobject.tree.Tree` object in a `depth-first <http://en.wikipedia.org/wiki/Depth-first_traversal>`_ search pattern, i.e., "visiting" a node before visiting the children of the node. This is the same traversal order as the previous example. This traversal order is useful if you require ancestral nodes to be processed before descendent nodes, as, for example, when evolving sequences over a tree.

    :meth:`~dendropy.dataobject.tree.Tree.postorder_node_iter()`
        Iterates over nodes in a :class:`~dendropy.dataobject.tree.Tree` object in a `postorder <http://en.wikipedia.org/wiki/Tree_traversal>`_ search pattern, i.e., visiting the children of the node before visiting the node itself. This traversal order is useful if you require descendent nodes to be processed before ancestor nodes, as, for example, when calculating ages of nodes.

    :meth:`~dendropy.dataobject.tree.Tree.level_order_node_iter()`
        Iterates over nodes in a :class:`~dendropy.dataobject.tree.Tree` object in a  `breadth-first <http://en.wikipedia.org/wiki/Breadth-first_traversal>`_  search pattern, i.e., every node at a particular level is visited before proceeding to the next level.

    :meth:`~dendropy.dataobject.tree.Tree.leaf_iter()`
        Iterates over the leaf or tip nodes of a :class:`~dendropy.dataobject.tree.Tree` object.

The previous example would thus be better implemented as follows:

.. literalinclude:: /examples/tree_evolve_char2.py
    :linenos:

The nodes returned by each of these iterators can be filtered if a filter function is passed as a second argument to the iterator.
This filter function should take a :class:`~dendropy.dataobject.tree.Node` object as an argument, and return :keyword:`True` if the node is to be returned or :keyword:`False` if it is not. For example, the following iterates over all nodes that have more than two children:

.. literalinclude:: /examples/preorder_filtered_node_iteration.py
    :linenos:

Iterating Over Edges
--------------------

The :class:`~dendropy.dataobject.tree.Edge` objects associated with each :class:`~dendropy.dataobject.tree.Node` can be accessed through the :attr:`~dendropy.dataobject.tree.Node.edge` attribute of the :class:`~dendropy.dataobject.tree.Node` object.
So it is possible to iterate over every edge on a tree by iterating over the nodes and referencing the :attr:`~dendropy.dataobject.tree.Node.edge` attribute of the node when processing the node.
But it is clearer and probably more convenient to use one of the :class:`~dendropy.dataobject.tree.Edge` iterators:

    :meth:`~dendropy.dataobject.tree.Tree.preorder_edge_iter()`
        Iterates over edges in a :class:`~dendropy.dataobject.tree.Tree` object in a `depth-first <http://en.wikipedia.org/wiki/Depth-first_traversal>`_ search pattern, i.e., "visiting" an edge before visiting the edges descending from that edge. This is the same traversal order as the previous example. This traversal order is useful if you require ancestral edges to be processed before descendent edges, as, for example, when calculating the sum of edge lengths from the root.

    :meth:`~dendropy.dataobject.tree.Tree.postorder_edge_iter()`
        Iterates over edges in a :class:`~dendropy.dataobject.tree.Tree` object in a `postorder <http://en.wikipedia.org/wiki/Tree_traversal>`_ search pattern, i.e., visiting the descendents of the edge before visiting the edge itself. This traversal order is useful if you require descendent edges to be processed before ancestral edges, as, for example, when calculating the sum of edge lengths from the tip

    :meth:`~dendropy.dataobject.tree.Tree.level_order_edge_iter()`
        Iterates over edges in a :class:`~dendropy.dataobject.tree.Tree` object in a  `breadth-first <http://en.wikipedia.org/wiki/Breadth-first_traversal>`_  search pattern, i.e., every edge at a particular level is visited before proceeding to the next level.

The following example sets the edge lengths of a tree to the proportions of the total tree length that they represent:

.. literalinclude:: /examples/rescale_tree_length.py
    :linenos:

While this one removes the edge lengths entirely:

.. literalinclude:: /examples/remove_branch_lengths.py
    :linenos:

Like the node iterators, the edge iterators also optionally take a filter function as a second argument, except here the filter function should take an :class:`~dendropy.dataobject.tree.Edge` object as an argument.
The following example shows how you might iterate over all edges with lengths less than some value:

.. literalinclude:: /examples/preorder_filtered_edge_iteration.py
    :linenos:

Finding Nodes on Trees
======================

Nodes with Taxa
---------------

To retrieve a node associated with a particular taxon, we can use the :meth:`~dendropy.dataobject.tree.Tree.find_taxon_node()` method, which takes a filter function as an argument.
The filter function should take a :class:`~dendropy.dataobject.taxon.Taxon` object as an argument and return :keyword:`True` if the taxon is to be returned.
For example:

.. literalinclude:: /examples/find_taxon_node1.py
    :linenos:

Because we might find it easier to refer to :class:`~dendropy.dataobject.taxon.Taxon` objects by their labels, a convenience method that wraps the retrieval of nodes associated with :class:`~dendropy.dataobject.taxon.Taxon` objects of particular label is provided:

.. literalinclude:: /examples/find_taxon_node2.py
    :linenos:

Most Recent Common Ancestors
----------------------------

The MRCA (most recent common ancestor) of taxa or nodes can be retrieved by the instance method :meth:`~dendropy.dataobject.tree.Tree.mrca()`.
This method takes a list of :class:`~dendropy.dataobject.taxon.Taxon` objects given by the ``taxa`` keyword argument, or a list of taxon labels given by the ``taxon_labels`` keyword argument, and returns a :class:`~dendropy.dataobject.tree.Node` object that corresponds to the MRCA of the specified taxa.
For example:

.. literalinclude:: /examples/mrca.py
    :linenos:

Note that this method is inefficient when you need to resolve MRCA's for multiple sets or pairs of taxa.
In this context, the :class:`~dendropy.treecalc.PatristicDistanceMatrix` offers a more efficient approach, and should be preferred for applications such as calculating the patristic distances between all pairs of taxa.

Tree Metrics
============

Tree Length
-----------

The :meth:`~dendropy.dataobject.tree.Tree.length()` method returns the sum of edge lengths of a :class:`~dendropy.dataobject.tree.Tree` object, with edges that do not have any length assigned being treated as edges with length 0.
The following example shows how to identify the "critical" value for an `Archie-Faith-Cranston or PTP test <http://hymenoptera.tamu.edu/courses/ento606/Suggested%20Readings/Slowinksi_Crother_1998.pdf>`_ from a sample of :class:`~dendropy.dataobject.tree.Tree` objects, i.e. a tree length equal to or greater than 95% of the trees in the sample:

.. literalinclude:: /examples/tree_length_crit.py
    :linenos:

Node Ages
---------

Pybus-Harvey Gamma
-------------------

Patristic Distances
-------------------


.. SCRATCH
    Each :class:`~dendropy.dataobject.tree.Tree` object has an attribute, :attr:`~dendropy.dataobject.tree.Tree.taxon_set`, which is a ``TaxaBlock`` object, and manages all the :class:`~dendropy.dataobject.taxon.Taxon` objects associated with the tree.
    The ``TaxaBlock`` object referenced by a :class:`~dendropy.dataobject.tree.Tree` object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set` might be shared by many other elements of the dataset, including other :class:`~dendropy.dataobject.tree.Tree` objects and :class:`~dendropy.dataobject.char.CharacterArray` objects, so any modification of elements of a :class:`~dendropy.dataobject.tree.Tree` object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set` will probably have dataset-wide effects.
    That is, if you were to change the label of a :class:`~dendropy.dataobject.taxon.Taxon` object maintained by a particular :class:`~dendropy.dataobject.tree.Tree` object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set`, all other :class:`~dendropy.dataobject.tree.Tree` objects in the dataset referencing the same ``TaxaBlock`` will be effected.



    Every :class:`~dendropy.dataobject.tree.Node` object also has an :attr:`~dendropy.dataobject.tree.Node.edge` attribute, which points to an :class:`~dendropy.dataobject.tree.Edge` object representing the branch subtending the node. :class:`~dendropy.dataobject.tree.Edge` objects have a :attr:`~dendropy.dataobject.tree.Edge.length` attribute, which is typically either a ``float`` or ``int`` value, representing the weight or length of the branch.
    If branch lengths have not been specified, then the value of :attr:`~dendropy.dataobject.tree.Edge.length` is :keyword:`~None`.
    Even if the source tree has had branch lengths specified, if the tree is unrooted, then the edge of the :attr:`~dendropy.dataobject.tree.Tree.seed_node` is usually :keyword:`~None`.

    :class:`~dendropy.dataobject.tree.Node` objects also have a :attr:`~dendropy.dataobject.tree.Node.label` and :attr:`~dendropy.dataobject.tree.Node.taxon` attribute. Leaf nodes usually have their :attr:`~dendropy.dataobject.tree.Node.taxon` attribute set, pointing to :class:`~dendropy.dataobject.taxon.Taxon` object associated with that tip of the tree. The :attr:`~dendropy.dataobject.tree.Node.label` attribute will be set if the source tree has internal node labels, though, of course, you can also assign a value to this programmatically.
