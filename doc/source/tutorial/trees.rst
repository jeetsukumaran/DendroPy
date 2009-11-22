*****************
Trees in DendroPy
*****************

Trees in are represented by the class :class:`~dendropy.dataobject.tree.Tree`.
Every :class:`~dendropy.dataobject.tree.Tree` object has a :attr:`~dendropy.dataobject.tree.Tree.seed_node` attribute. If the tree is rooted, then this is the root node. If the tree is not rooted, however, then this is an artificial node that serves as the "starting point" for the tree.
The :attr:`~dendropy.dataobject.tree.Tree.seed_node`, like every other node on the tree, is a :class:`~dendropy.dataobject.tree.Node` object.
Every :class:`~dendropy.dataobject.tree.Node` object maintains a list of its immediate child :class:`~dendropy.dataobject.tree.Node` objects as well as a reference to its parent :class:`~dendropy.dataobject.tree.Node` object.
You can request a shallow-copy :func:`~list` of child :class:`~dendropy.dataobject.tree.Node` objects using the :meth:`~dendropy.dataobject.tree.Node.child_nodes()` method, and you can access the parent :class:`~dendropy.dataobject.tree.Node` object directly through the :attr:`~dendropy.dataobject.tree.Node.parent_node` attribute.
By definition, the :attr:`~dendropy.dataobject.tree.Tree.seed_node` has no parent node, leaf nodes have no child nodes, and internal nodes have both parent nodes and child nodes.

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
        Iterates over the leaf or tip nodes :class:`~dendropy.dataobject.tree.Tree` object.

The previous example would thus be better implemented as follows:

.. literalinclude:: /examples/tree_evolve_char2.py
    :linenos:

The nodes returned by each of these iterators can be filtered if a filter function is passed as a second argument to the iterator.
This filter function should take a :class:`~dendropy.dataobject.tree.Node` object as an argument, and return :keyword:`True` if the node is to be returned or :keyword:`False` if it is not. For example, the following iterates over all nodes that have more than two children:

.. literalinclude:: /examples/preorder_filtered_node_iteration.py

Iterating Over Edges
--------------------

The :class:`~dendropy.dataobject.tree.Edge` objects associated with each :class:`~dendropy.dataobject.tree.Node` can be accessed through the :attr:`~dendropy.dataobject.tree.Node.edge` attribute of the :class:`~dendropy.dataobject.tree.Node` object.
So it is possible to iterate over every edge on a tree by iterating over the nodes and referencing the :attr:`~dendropy.dataobject.tree.Node.edge` attribute of the node when processing the node.
But it clearer and probably more convenient to use one of the :class:`~dendropy.dataobject.tree.Edge` iterators:

    :meth:`~dendropy.dataobject.tree.Tree.preorder_edge_iter()`
        Iterates over edges in a :class:`~dendropy.dataobject.tree.Tree` object in a `depth-first <http://en.wikipedia.org/wiki/Depth-first_traversal>`_ search pattern, i.e., "visiting" an edge before visiting the edges descending from that edge. This is the same traversal order as the previous example. This traversal order is useful if you require ancestral edges to be processed before descendent edges, as, for example, when calculating the sum of edge lengths from the root.

    :meth:`~dendropy.dataobject.tree.Tree.postorder_edge_iter()`
        Iterates over edges in a :class:`~dendropy.dataobject.tree.Tree` object in a `postorder <http://en.wikipedia.org/wiki/Tree_traversal>`_ search pattern, i.e., visiting the descendents of the edge before visiting the edge itself. This traversal order is useful if you require descendent edges to be processed before ancestral edges, as, for example, when calculating the sum of edge lengths from the tip

    :meth:`~dendropy.dataobject.tree.Tree.level_order_edge_iter()`
        Iterates over edges in a :class:`~dendropy.dataobject.tree.Tree` object in a  `breadth-first <http://en.wikipedia.org/wiki/Breadth-first_traversal>`_  search pattern, i.e., every edge at a particular level is visited before proceeding to the next level.


.. SCRATCH
    Each :class:`~dendropy.dataobject.tree.Tree` object has an attribute, :attr:`~dendropy.dataobject.tree.Tree.taxon_set`, which is a ``TaxaBlock`` object, and manages all the :class:`~dendropy.dataobject.taxon.Taxon` objects associated with the tree.
    The ``TaxaBlock`` object referenced by a :class:`~dendropy.dataobject.tree.Tree` object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set` might be shared by many other elements of the dataset, including other :class:`~dendropy.dataobject.tree.Tree` objects and :class:`~dendropy.dataobject.char.CharacterArray` objects, so any modification of elements of a :class:`~dendropy.dataobject.tree.Tree` object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set` will probably have dataset-wide effects.
    That is, if you were to change the label of a :class:`~dendropy.dataobject.taxon.Taxon` object maintained by a particular :class:`~dendropy.dataobject.tree.Tree` object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set`, all other :class:`~dendropy.dataobject.tree.Tree` objects in the dataset referencing the same ``TaxaBlock`` will be effected.



    Every :class:`~dendropy.dataobject.tree.Node` object also has an :attr:`~dendropy.dataobject.tree.Node.edge` attribute, which points to an :class:`~dendropy.dataobject.tree.Edge` object representing the branch subtending the node. :class:`~dendropy.dataobject.tree.Edge` objects have a :attr:`~dendropy.dataobject.tree.Edge.length` attribute, which is typically either a ``float`` or ``int`` value, representing the weight or length of the branch.
    If branch lengths have not been specified, then the value of :attr:`~dendropy.dataobject.tree.Edge.length` is :keyword:`~None`.
    Even if the source tree has had branch lengths specified, if the tree is unrooted, then the edge of the :attr:`~dendropy.dataobject.tree.Tree.seed_node` is usually :keyword:`~None`.

    :class:`~dendropy.dataobject.tree.Node` objects also have a :attr:`~dendropy.dataobject.tree.Node.label` and :attr:`~dendropy.dataobject.tree.Node.taxon` attribute. Leaf nodes usually have their :attr:`~dendropy.dataobject.tree.Node.taxon` attribute set, pointing to :class:`~dendropy.dataobject.taxon.Taxon` object associated with that tip of the tree. The :attr:`~dendropy.dataobject.tree.Node.label` attribute will be set if the source tree has internal node labels, though, of course, you can also assign a value to this programmatically.
