*******************************************
Working with Trees and Collections of Trees
*******************************************

Trees
-----

Trees in are represented by the class |Tree|.
Every |Tree| object has a :attr:`~dendropy.dataobject.tree.Tree.seed_node` attribute. If the tree is rooted, then this is the root node. If the tree is not rooted, however, then this is an artificial node that serves as the "starting point" for the tree.
The :attr:`~dendropy.dataobject.tree.Tree.seed_node`, like every other node on the tree, is a |Node| object.
Every |Node| object maintains a list of its immediate child |Node| objects as well as a reference to its parent |Node| object.
You can request a shallow-copy :func:`~list` of child |Node| objects using the :meth:`~dendropy.dataobject.tree.Node.child_nodes()` method, and you can access the parent |Node| object directly through the :attr:`~dendropy.dataobject.tree.Node.parent_node` attribute.
By definition, the :attr:`~dendropy.dataobject.tree.Tree.seed_node` has no parent node, leaf nodes have no child nodes, and internal nodes have both parent nodes and child nodes.

Tree Lists
----------

|TreeList| objects are lists of |Tree| objects constrained to sharing the same |TaxonSet|.
Any |Tree| object added to a |TreeList| will have its :attr:`~dendropy.dataobject.tree.Tree.taxon_set` attribute assigned to the |TaxonSet| object of the |TreeList|, and all referenced |Taxon| objects will be mapped to the same or corresponding |Taxon| objects of this new |TaxonSet|, with new |Taxon| objects created if no suitable match is found.

.. _Customizing_Tree_Creation_and_Reading:

Customizing |Tree| and |TreeList| Creation and Reading
======================================================

Both the |Tree| and |TreeList| classes support the :meth:`get_from_*()` factory class methods for simultaneously instantiating and populating objects, as well as the  :meth:`read_from_*()` instance methods for populating existing objects.
In the case of |Tree| objects, calling :meth:`read_from_*()` **repopulates** (i.e., redefines) the |Tree| with data from the data source, while in the case of |TreeList| objects, calling :meth:`read_from_*()` **appends** the tree definitions in the data source to the |TreeList|.

Using a Specific |TaxonSet|
---------------------------
Passing a |TaxonSet| object using the ``taxon_set`` argument when instantiating a |Tree| or |TreeList| object (using, for example, the meth:`get_from_*()` or :meth:`read_from_*()` methods) results in the |Tree| or |TreeList| object being bound to the specified |TaxonSet| object.


Selecting Specific Trees or Subsets of Trees
--------------------------------------------

The ``tree_offset`` and ``collection_offset`` keywords allow you to control which tree defintions are parsed from the data source:

    ``tree_offset``
        A non-negative integer specifying the 0-based index of a tree within a collection in the data source.
        The default is 0, which means that the first tree definition is used.
        If passed to :meth:`get_from_*()`, :meth:`read_from_*()` or a constructor of |Tree|, this selects a specific tree definition in the source (i.e, ``tree_offset=2`` will create or populate the |Tree| object based on the 3rd tree definition). If passed to  :meth:`get_from_*()`, :meth:`read_from_*()` or a constructor of |TreeList| or |DataSet| object, this effectively skips all the tree definitions preceding the specified index from being created (i.e, ``tree_offset=200`` will populate the |TreeList| object starting with the 201st tree definition).


        For example, the following creates a |Tree| object from the second tree definition in the data source::

            >>> import dendropy
            >>> t = dendropy.Tree.get_from_path('pythonidae.best-trees.tre', \
                        'nexus', tree_offset=1)

        While this effectively skips over the first 200 trees as burn-in from an MCMC sample of trees::

            >>> import dendropy
            >>> pp_trees = dendropy.TreeList.get_from_path('pythonidae_mcmc.tre', \
                    'nexus', tree_offset=200)

    ``collection_offset``
        A non-negative integer specifying the 0-based index of a collection (e.g., a NEXUS "TREES" block) of trees in the data source.
        A negative value means that a union of all the tree collections in the data source will be used.
        The default is -1, i.e., all the collections will be aggregated.
        For example, the following selects the third tree collection to populate a |TreeList| object::

            >>> import dendropy
            >>> trees = dendropy.Tree.get_from_path('pythonidae.nex', 'nexus', \
                    collection_offset=4)

        While this reads all the trees from all "TREES" block in the data source::

            >>> import dendropy
            >>> trees = dendropy.TreeList.get_from_path('pythonidae.nex', 'nexus', \
                    collection_offset=-1)

        The following selects the second tree from the third "TREES" block in the data source::

            >>> import dendropy
            >>> trees = dendropy.Tree.get_from_path('pythonidae.nex', 'nexus', \
                    collection_offset=2, tree_offset=1)

        The following selects the 30th tree defined in the data source across all tree collections, with the first tree in the first collection treated as having index 0::

            >>> import dendropy
            >>> tree_31 = dendropy.Tree.get_from_path('pythonidae.nex', 'nexus', \
                    collection_offset=-1, tree_offset=29)

Interpreting Rootings
---------------------

The rooting state of a |Tree| object is set by the :attr:`~dendropy.dataobject.tree.Tree.is_rooted` property.
When parsing NEXUS- and Newick-formatted data, the rooting states of the resulting |Tree| objects are given by ``[&R]`` (for rooted) or ``[&U]`` (for unrooted) comment tags preceding the tree definition in the data source.
If these tags are not present, then the trees are assumed to be unrooted.
This behavior can be changed by specifying keyword arguments to the :meth:`get_from_*()`,  or :meth:`read_from_*()` methods of both the |Tree| and |TreeList| classes, or the constructors of these classes when specifying a data source from which to construct the tree:

The ``as_rooted`` keyword argument, if :keyword:`True`, forces all trees to be interpreted as rooted, regardless of whether or not the ``[&R]``/``[&U]`` comment tags are given.
Conversely, if :keyword:`False`, all trees will be interpreted as unrooted.
For semantic clarity, you can also specify ``as_unrooted`` to be :keyword:`True` to force all trees to be unrooted.

.. literalinclude:: /examples/tree_rootings1.py
    :linenos:

In addition, you can specify a ``default_as_rooted`` keyword argument, which, if :keyword:`True`, forces all trees to be interpreted as rooted, *if* the ``[&R]``/``[&U]`` comment tags are *not* given.
Otherwise the rooting will follow the ``[&R]``/``[&U]`` commands.
Conversely, if ``default_as_rooted`` is :keyword:`False`, all trees will be interpreted as unrooted if the ``[&R]``/``[&U]`` comment tags are not given.
Again, for semantic clarity, you can also specify ``default_as_unrooted`` to be :keyword:`True` to assume all trees are unrooted if not explicitly specified, though, as this is default behavior, this should not be neccessary.

Custom Handling of Underscores
------------------------------
With NEXUS and NEWICK data sources, you can also specify ``preserve_underscores=True``.
The NEXUS standard dictates that underscores are equivalent to spaces, and thus any underscore found in any unquoted label in a NEXUS/NEWICK data source will be substituted for spaces.
Specifying ``preserve_underscores=True`` will force DendroPy to keep the underscores. More details on using this keyword to manage taxon references and mapping can be found in here: :ref:`Taxon_Label_Mapping`.


Efficiently Iterating Over Trees in a File
==========================================

If you need to process a collection of trees defined in a file source, you can read the trees into a |TreeList| object and iterate over the resulting collection::

    >>> import dendropy
    >>> trees = dendropy.TreeList.get_from_path('pythonidae.beast-mcmc.trees', 'nexus')
    >>> for tree in trees:
    ...     print(tree.as_string('newick'))

In the above, the entire data source is parsed and stored in the ``trees`` object before being processed in the subsequent lines.
In some cases, you might not need to maintain all the trees in memory at the same time.
For example, you might be interested in calculating the distribution of a statistic over a collection of trees, but have no need to refer to any of the trees after the statistic has been calculated.
In this case, it might be more efficient to use the :func:`~dendropy.dataio.tree_source_iter()` function.
This takes a file-like object as its first argument and a schema specification as the second and returns an iterator over the trees in the file.
Additional keyword arguments to customize the parsing are the same as that for the general :meth:`get_from_*()` and :meth:`read_from_*()` methods.
For example, the following script reads a model tree from a file, and then iterates over a collection of MCMC trees in another file, calculating a storing the symmetric distance between the model tree and each of the MCMC trees one at time:

.. literalinclude:: /examples/tree_iter1.py
    :linenos:

Note how a |TaxonSet| object is created and passed to both the :meth:`~dendropy.dataobject.Tree.get_from_path()` and the :func:`~dendropy.dataio.tree_source_iter()` functions using the ``taxon_set`` keyword argument.
This is to ensure that the corresponding taxa in both sources get mapped to the same |Taxon| objects in DendroPy object space, so as to enable comparisons of the trees.
If this was not done, then each tree would have its own distinct |TaxonSet| object (and associated |Taxon| objects), making comparisons impossible.

Also note how the ``tree_offset`` keyword is used to skip over the burn-in trees from the MCMC sample.

If you want to iterate over trees in multiple sources, you can use the :func:`~dendropy.dataio.multi_tree_source_iter()`.
This takes a list of file-like objects *or* a list of filepath strings as its first argument, and a schema-specification string as its second argument.
Again, other keyword arguments supported by the general :meth:`get_from_*()` and :meth:`read_from_*()` methods are also available.

For example:

.. literalinclude:: /examples/tree_iter2.py
    :linenos:

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
This filter function should take a |Node| object as an argument, and return :keyword:`True` if the node is to be returned or :keyword:`False` if it is not. For example, the following iterates over all nodes that have more than two children:

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
The filter function should take a |Taxon| object as an argument and return :keyword:`True` if the taxon is to be returned.
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

Tree Metrics
============

Tree Length
-----------

The :meth:`~dendropy.dataobject.tree.Tree.length()` method returns the sum of edge lengths of a |Tree| object, with edges that do not have any length assigned being treated as edges with length 0.
The following example shows how to identify the "critical" value for an `Archie-Faith-Cranston or PTP test <http://hymenoptera.tamu.edu/courses/ento606/Suggested%20Readings/Slowinksi_Crother_1998.pdf>`_ from a sample of |Tree| objects, i.e. a tree length equal to or greater than 95% of the trees in the sample:

.. literalinclude:: /examples/tree_length_crit.py
    :linenos:

Node Ages
---------

The :meth:`~dendropy.dataobject.tree.Tree.add_ages_to_nodes()` method calculates the age of a node (i.e., the sum of edge lengths from the node to a tip) and assigns it to a new attribute of the node: :attr:`~dendropy.dataobject.tree.Node.age`. The following example iterates through the post-burn-in of an MCMC sample of ultrametric trees, calculating the age of the MRCA of two taxa, and reports the mean age of the node.

.. literalinclude:: /examples/node_ages1.py
    :linenos:

Pybus-Harvey Gamma
------------------

The Pybus-Harvey Gamma statistic is given by the :meth:`~dendropy.dataobject.tree.Tree.pybus_harvey_gamma()` instance method. The following example iterates through the post-burn-in of an MCMC sample of trees, reporting the mean Pybus-Harvey Gamma statistic:

.. literalinclude:: /examples/pbhg.py
    :linenos:

Patristic Distances
-------------------

The :class:`~dendropy.treecalc.PatristicDistanceMatrix` is the most efficient way to calculate the patristic distances between any pair of taxa on a tree.
Its constructor takes a |Tree| object as an argument, and the object return is callable, taking two |Taxon| objects as arguments and returning the sum of edge lengths between the two. The following example reports the pairwise distances between all taxa on the input tree:

.. literalinclude:: /examples/pdm.py
    :linenos:


.. SCRATCH
    Each |Tree| object has an attribute, :attr:`~dendropy.dataobject.tree.Tree.taxon_set`, which is a ``TaxaBlock`` object, and manages all the |Taxon| objects associated with the tree.
    The ``TaxaBlock`` object referenced by a |Tree| object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set` might be shared by many other elements of the dataset, including other |Tree| objects and |CharacterMatrix| objects, so any modification of elements of a |Tree| object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set` will probably have dataset-wide effects.
    That is, if you were to change the label of a |Taxon| object maintained by a particular |Tree| object's :attr:`~dendropy.dataobject.tree.Tree.taxon_set`, all other |Tree| objects in the dataset referencing the same ``TaxaBlock`` will be effected.



    Every |Node| object also has an :attr:`~dendropy.dataobject.tree.Node.edge` attribute, which points to an |Edge| object representing the branch subtending the node. |Edge| objects have a :attr:`~dendropy.dataobject.tree.Edge.length` attribute, which is typically either a ``float`` or ``int`` value, representing the weight or length of the branch.
    If branch lengths have not been specified, then the value of :attr:`~dendropy.dataobject.tree.Edge.length` is :keyword:`~None`.
    Even if the source tree has had branch lengths specified, if the tree is unrooted, then the edge of the :attr:`~dendropy.dataobject.tree.Tree.seed_node` is usually :keyword:`~None`.

    |Node| objects also have a :attr:`~dendropy.dataobject.tree.Node.label` and :attr:`~dendropy.dataobject.tree.Node.taxon` attribute. Leaf nodes usually have their :attr:`~dendropy.dataobject.tree.Node.taxon` attribute set, pointing to |Taxon| object associated with that tip of the tree. The :attr:`~dendropy.dataobject.tree.Node.label` attribute will be set if the source tree has internal node labels, though, of course, you can also assign a value to this programmatically.
