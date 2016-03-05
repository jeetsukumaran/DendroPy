*****
Trees
*****

The |Tree| Class
================

:term:`Trees <tree>` in DendroPy are represented by objects of the class |Tree|.
:term:`Trees <tree>` consist of a collection of :term:`nodes <node>`, represented by objects of the class |Node|, connected to each other in parent-child or ancestor-descendent relationships by objects of the class |Edge|.
The first or initial :term:`node` of a |Tree| is the head of the data structure, and is represented by the :attr:`seed_node` attribute of a |Tree| object.
If the tree is :term:`rooted <rooted tree>`, then this is the :term:`root node`.
If the tree is :term:`unrooted <unrooted tree>`, however, then this is an artificial node that *only* serves as the initial node when iterating over a tree in :term:`preorder <preorder traversal>` or the final node when iterating over the tree in :term:`postorder <postorder traversal>`, but does not have any informational significance by itself: it is an algorithmic artifact.

The :attr:`~dendropy.datamodel.treemodel.Tree.seed_node` attribute of a |Tree| object, like every other node on the tree, is an object of the |Node| class.
Every |Node| object maintains a list of its immediate child |Node| objects as well as a reference to its parent |Node| object.
You can iterate over the :term:`child <child node>` of a particular |Node| object using the :meth:`~dendropy.datamodel.treemodel.Node.child_node_iter()` method, get a shallow-copy list of child |Node| objects using the :meth:`~dendropy.datamodel.treemodel.Node.child_nodes()` method, or access the :term:`parent <parent node>` |Node| object directly through the :attr:`~dendropy.datamodel.treemodel.Node.parent_node` attribute of the |Node|.
By definition, the :attr:`~dendropy.datamodel.treemodel.Tree.seed_node` has no :term:`parent node`, :term:`leaf nodes <leaf node>` have no :term:`child nodes <child node>`, and term:`internal nodes <internal node>` have both :term:`parent nodes <parent node>` and :term:`child nodes <child node>`.

Every |Node| object has an attribute, :attr:`~dendropy.datamodel.treemodel.Node.edge`, which is an |Edge| object representing the :term:`edge` that is :term:`incident to or subtends <incident edge>` the :term:`node` represented by that |Node| object.
Each |Edge|, in turn, has an attribute, :attr:`~dendropy.datamodel.treemodel.Edge.head_node`, which is the |Node| object representing the :term:`node` that the edge subtends.

The |Tree|, |Node|, and |Edge| classes all have "``annotations``" as an attribute, which is a :class:`~dendropy.datamodel.basemodel.AnnotationSet` object, i.e. a collection of :class:`~dendropy.datamodel.basemodel.Annotation` instances tracking metadata.
More information on working with metadata can be found in the ":doc:`/primer/working_with_metadata_annotations`" section.

Reading and Writing |Tree| Instances
====================================

The |Tree| class supports the ":meth:`~dendropy.datamodel.treemodel.Tree.get`" factory class method for simultaneously instantiating and populating a |Tree| instance, taking a data source as the first argument and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` ("``nexus``", "``newick``", "``nexml``", "``fasta``", or "``phylip``", etc.) as the second::

    import dendropy
    tree = dendropy.Tree.get(
        path='pythonidae.mcmc.nex',
        schema='nexus')

A |Tree| object can be written to an external resource using the ":meth:`~dendropy.datamodel.treemodel.Tree.write`" method::

    import dendropy
    tree = dendropy.Tree.get(
        path="trees1.nex",
        schema="nexus",
        tree_offset=2,
        )
    tree.write(
        path="trees1.newick",
        schema="newick",
        )

It can also be represented as a string using the ":meth:`~dendropy.datamodel.treemodel.Tree.as_string`" method::

    import dendropy
    tree = dendropy.Tree.get(
        path="trees1.nex",
        schema="nexus",
        )
    print(tree.as_string(schema="newick",)

More information on reading operations is available in the :doc:`/primer/reading_and_writing` section.

Cloning/Copying a |Tree|
========================

You can make a "taxon namespace-scoped" copy of a |Tree| instance, i.e., where all |Node| and associated |Edge| instances of a |Tree| are cloned, but references to |Taxon| objects are preserved, you can call :meth:`dendropy.datamodel.treemodel.Tree.clone` with a "``depth``" argument value of 1 or by copy construction:

.. literalinclude:: /examples/tree_copy1.py

For a true and complete deep-copy, where even the |Taxon| and |TaxonNamespace| references are copied, call :func:`copy.deepcopy`:

.. literalinclude:: /examples/tree_copy2.py

Alternatively, many times you want a "light" or "thin" copy, where just the tree structure (node and edge relationships) and basic information are retained (e.g., edge lengths, taxon associations, and node and edge labels), but not, e.g. the rich annotations. Then the :meth:`~dendropy.datamodel.treemode.Tree.extract_tree` method is what you are looking for:

.. literalinclude:: /examples/tree_extract1.py

Tree Traversal
==============

Iterating Over Nodes
--------------------

The following example shows how you might evolve a continuous character on a tree by recursively visting each node, and setting the value of the character to one drawn from a normal distribution centered on the value of the character of the node's ancestor and standard deviation given by the length of the edge subtending the node:

.. literalinclude:: /examples/tree_evolve_char1.py
    :linenos:

While the previous example works, it is probably clearer and more efficient to use one of the pre-defined node iterator methods:

    :meth:`~dendropy.datamodel.treemodel.Tree.preorder_node_iter()`
        Iterates over nodes in a |Tree| object in a `depth-first <http://en.wikipedia.org/wiki/Depth-first_traversal>`_ search pattern, i.e., "visiting" a node before visiting the children of the node. This is the same traversal order as the previous example. This traversal order is useful if you require ancestral nodes to be processed before descendent nodes, as, for example, when evolving sequences over a tree.

    :meth:`~dendropy.datamodel.treemodel.Tree.postorder_node_iter()`
        Iterates over nodes in a |Tree| object in a `postorder <http://en.wikipedia.org/wiki/Tree_traversal>`_ search pattern, i.e., visiting the children of the node before visiting the node itself. This traversal order is useful if you require descendent nodes to be processed before ancestor nodes, as, for example, when calculating ages of nodes.

    :meth:`~dendropy.datamodel.treemodel.Tree.level_order_node_iter()`
        Iterates over nodes in a |Tree| object in a  `breadth-first <http://en.wikipedia.org/wiki/Breadth-first_traversal>`_  search pattern, i.e., every node at a particular level is visited before proceeding to the next level.

    :meth:`~dendropy.datamodel.treemodel.Tree.leaf_node_iter()`
        Iterates over the leaf or tip nodes of a |Tree| object.

The previous example would thus be better implemented as follows:

.. literalinclude:: /examples/tree_evolve_char2.py

The nodes returned by each of these iterators can be filtered if a filter function is passed as a second argument to the iterator.
This filter function should take a |Node| object as an argument, and return |True| if the node is to be returned or |False| if it is not. For example, the following iterates over all nodes that have more than two children:

.. literalinclude:: /examples/preorder_filtered_node_iteration.py
    :linenos:

Iterating Over Edges
--------------------

The |Edge| objects associated with each |Node| can be accessed through the :attr:`~dendropy.datamodel.treemodel.Node.edge` attribute of the |Node| object.
So it is possible to iterate over every edge on a tree by iterating over the nodes and referencing the :attr:`~dendropy.datamodel.treemodel.Node.edge` attribute of the node when processing the node.
But it is clearer and probably more convenient to use one of the |Edge| iterators:

    :meth:`~dendropy.datamodel.treemodel.Tree.preorder_edge_iter()`
        Iterates over edges in a |Tree| object in a `depth-first <http://en.wikipedia.org/wiki/Depth-first_traversal>`_ search pattern, i.e., "visiting" an edge before visiting the edges descending from that edge. This is the same traversal order as the previous example. This traversal order is useful if you require ancestral edges to be processed before descendent edges, as, for example, when calculating the sum of edge lengths from the root.

    :meth:`~dendropy.datamodel.treemodel.Tree.postorder_edge_iter()`
        Iterates over edges in a |Tree| object in a `postorder <http://en.wikipedia.org/wiki/Tree_traversal>`_ search pattern, i.e., visiting the descendents of the edge before visiting the edge itself. This traversal order is useful if you require descendent edges to be processed before ancestral edges, as, for example, when calculating the sum of edge lengths from the tip

    :meth:`~dendropy.datamodel.treemodel.Tree.level_order_edge_iter()`
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

To retrieve a node associated with a particular taxon, we can use the :meth:`~dendropy.datamodel.treemodel.Tree.find_taxon_node()` method, which takes a filter function as an argument.
The filter function should take a |Taxon| object as an argument and return |True| if the taxon is to be returned.
For example:

.. literalinclude:: /examples/find_taxon_node1.py

Because we might find it easier to refer to |Taxon| objects by their labels, a convenience method that wraps the retrieval of nodes associated with |Taxon| objects of particular label is provided:

.. literalinclude:: /examples/find_taxon_node2.py

Most Recent Common Ancestors
----------------------------

The MRCA (most recent common ancestor) of taxa or nodes can be retrieved by the instance method :meth:`~dendropy.datamodel.treemodel.Tree.mrca()`.
This method takes a list of |Taxon| objects given by the ``taxa`` keyword argument, or a list of taxon labels given by the ``taxon_labels`` keyword argument, and returns a |Node| object that corresponds to the MRCA of the specified taxa.
For example:

.. literalinclude:: /examples/mrca.py

Note that this method is inefficient when you need to resolve MRCA's for multiple sets or pairs of taxa.
In this context, the :class:`~dendropy.calculate.treemeasure.PhylogeneticDistanceMatrix` offers a more efficient approach, and should be preferred for applications such as calculating the patristic distances between all pairs of taxa. An instance of this class will be returned when you call :meth:`~dendropy.datamodel.treemodel.Tree.phylogenetic_distance_matrix()`:

.. literalinclude:: /examples/mrca2.py

Note that the |PhylogeneticDistanceMatrix| object does not automatically update if the original |Tree| changes: it is essentially a snapshot of |Tree| at the point in which it is instantiated.
If the original |Tree| changes, you should create a new instance of the corresponding |PhylogeneticDistanceMatrix| object.

Viewing and Displaying Trees
============================

Sometimes it is useful to get a visual representation of a |Tree|.

For quick inspection, the :meth:`~dendropy.datamodel.treemodel.Tree.print_plot()` will write an ASCII text plot to the standard output stream::

    >>> t = dendropy.Tree.get_from_string("(A,(B,(C,D)))", "newick")
    >>> t.print_plot()
    /----------------------------------------------- A
    +
    |                /------------------------------ B
    \----------------+
                     |          /------------------- C
                     \----------+
                                \------------------- D

If you need to store this representation as a string instead, you can use :meth:`~dendropy.datamodel.treemodel.Tree.as_ascii_plot()`::

    >>> s = t.as_ascii_plot()
    >>> print(s)
    /----------------------------------------------- A
    +
    |                /------------------------------ B
    \----------------+
                     |          /------------------- C
                     \----------+
                                \------------------- D

You can also, as mentioned above, using the :meth:`~dendropy.datamodel.treemodel.Tree.as_string` method to represent a |Tree| as string in any format::

    t = dendropy.Tree.get_from_string("(A,(B,(C,D)))", "newick")
    print(t.as_string(schema="nexus"))
    print(t.as_string(schema="newick"))


Building a Tree Programmatically
================================

For example:

.. literalinclude:: /examples/build_tree_programmatically.py

produces the following::

    ((A:1,B:2):1,(C:1,D:2):1);

                                       /---------------------------------- A
    /----------------------------------+
    |                                  \---------------------------------- B
    +
    |                                  /---------------------------------- C
    \----------------------------------+
                                       \---------------------------------- D
