*****************************
Tree Traversal and Navigation
*****************************

Tree Traversal
==============

Iterating Over Nodes
--------------------

The following example shows how you might evolve a continuous character on a tree by recursively visting each node, and setting the value of the character to one drawn from a normal distribution centered on the value of the character of the node's ancestor and standard deviation given by the length of the edge subtending the node:

.. literalinclude:: /examples/tree_evolve_char1.py
    :linenos:

While the previous example works, it is probably clearer and more efficient to use one of the pre-defined node iterator methods:

    :meth:`~dendropy.dataobject.tree.Tree.preorder_node_iter()`
        Iterates over nodes in a |Tree| object in a `depth-first <http://en.wikipedia.org/wiki/Depth-first_traversal>`_ search pattern, i.e., "visiting" a node before visiting the children of the node. This is the same traversal order as the previous example. This traversal order is useful if you require ancestral nodes to be processed before descendent nodes, as, for example, when evolving sequences over a tree.

    :meth:`~dendropy.dataobject.tree.Tree.postorder_node_iter()`
        Iterates over nodes in a |Tree| object in a `postorder <http://en.wikipedia.org/wiki/Tree_traversal>`_ search pattern, i.e., visiting the children of the node before visiting the node itself. This traversal order is useful if you require descendent nodes to be processed before ancestor nodes, as, for example, when calculating ages of nodes.

    :meth:`~dendropy.dataobject.tree.Tree.level_order_node_iter()`
        Iterates over nodes in a |Tree| object in a  `breadth-first <http://en.wikipedia.org/wiki/Breadth-first_traversal>`_  search pattern, i.e., every node at a particular level is visited before proceeding to the next level.

    :meth:`~dendropy.dataobject.tree.Tree.leaf_iter()`
        Iterates over the leaf or tip nodes of a |Tree| object.

The previous example would thus be better implemented as follows:

.. literalinclude:: /examples/tree_evolve_char2.py
    :linenos:

The nodes returned by each of these iterators can be filtered if a filter function is passed as a second argument to the iterator.
This filter function should take a |Node| object as an argument, and return |True| if the node is to be returned or |False| if it is not. For example, the following iterates over all nodes that have more than two children:

.. literalinclude:: /examples/preorder_filtered_node_iteration.py
    :linenos:

Iterating Over Edges
--------------------

The |Edge| objects associated with each |Node| can be accessed through the :attr:`~dendropy.dataobject.tree.Node.edge` attribute of the |Node| object.
So it is possible to iterate over every edge on a tree by iterating over the nodes and referencing the :attr:`~dendropy.dataobject.tree.Node.edge` attribute of the node when processing the node.
But it is clearer and probably more convenient to use one of the |Edge| iterators:

    :meth:`~dendropy.dataobject.tree.Tree.preorder_edge_iter()`
        Iterates over edges in a |Tree| object in a `depth-first <http://en.wikipedia.org/wiki/Depth-first_traversal>`_ search pattern, i.e., "visiting" an edge before visiting the edges descending from that edge. This is the same traversal order as the previous example. This traversal order is useful if you require ancestral edges to be processed before descendent edges, as, for example, when calculating the sum of edge lengths from the root.

    :meth:`~dendropy.dataobject.tree.Tree.postorder_edge_iter()`
        Iterates over edges in a |Tree| object in a `postorder <http://en.wikipedia.org/wiki/Tree_traversal>`_ search pattern, i.e., visiting the descendents of the edge before visiting the edge itself. This traversal order is useful if you require descendent edges to be processed before ancestral edges, as, for example, when calculating the sum of edge lengths from the tip

    :meth:`~dendropy.dataobject.tree.Tree.level_order_edge_iter()`
        Iterates over edges in a |Tree| object in a  `breadth-first <http://en.wikipedia.org/wiki/Breadth-first_traversal>`_  search pattern, i.e., every edge at a particular level is visited before proceeding to the next level.

The following example sets the edge lengths of a tree to the proportions of the total tree length that they represent:

.. literalinclude:: /examples/rescale_tree_length.py
    :linenos:

While this one removes the edge lengths entirely:

.. literalinclude:: /examples/remove_branch_lengths.py
    :linenos:

Like the node iterators, the edge iterators also optionally take a filter function as a second argument, except here the filter function should take an |Edge| object as an argument.
The following example shows how you might iterate over all edges with lengths less than some value:

.. literalinclude:: /examples/preorder_filtered_edge_iteration.py
    :linenos:

Finding Nodes on Trees
======================

Nodes with Particular Taxa
--------------------------

To retrieve a node associated with a particular taxon, we can use the :meth:`~dendropy.dataobject.tree.Tree.find_taxon_node()` method, which takes a filter function as an argument.
The filter function should take a |Taxon| object as an argument and return |True| if the taxon is to be returned.
For example:

.. literalinclude:: /examples/find_taxon_node1.py
    :linenos:

Because we might find it easier to refer to |Taxon| objects by their labels, a convenience method that wraps the retrieval of nodes associated with |Taxon| objects of particular label is provided:

.. literalinclude:: /examples/find_taxon_node2.py
    :linenos:

Most Recent Common Ancestors
----------------------------

The MRCA (most recent common ancestor) of taxa or nodes can be retrieved by the instance method :meth:`~dendropy.dataobject.tree.Tree.mrca()`.
This method takes a list of |Taxon| objects given by the ``taxa`` keyword argument, or a list of taxon labels given by the ``taxon_labels`` keyword argument, and returns a |Node| object that corresponds to the MRCA of the specified taxa.
For example:

.. literalinclude:: /examples/mrca.py
    :linenos:

Note that this method is inefficient when you need to resolve MRCA's for multiple sets or pairs of taxa.
In this context, the :class:`~dendropy.treecalc.PatristicDistanceMatrix` offers a more efficient approach, and should be preferred for applications such as calculating the patristic distances between all pairs of taxa.
