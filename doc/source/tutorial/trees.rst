*******************************************
Working with Trees and Collections of Trees
*******************************************

Trees and Tree Lists
====================

Trees
-----

Trees in DendroPy are represented by the class |Tree|.
Every |Tree| object has a :attr:`~dendropy.dataobject.tree.Tree.seed_node` attribute. If the tree is rooted, then this is the root node. If the tree is not rooted, however, then this is an artificial node that serves as the "starting point" for the tree.
The :attr:`~dendropy.dataobject.tree.Tree.seed_node`, like every other node on the tree, is a |Node| object.
Every |Node| object maintains a list of its immediate child |Node| objects as well as a reference to its parent |Node| object.
You can request a shallow-copy :func:`~list` of child |Node| objects using the :meth:`~dendropy.dataobject.tree.Node.child_nodes()` method, and you can access the parent |Node| object directly through the :attr:`~dendropy.dataobject.tree.Node.parent_node` attribute.
By definition, the :attr:`~dendropy.dataobject.tree.Tree.seed_node` has no parent node, leaf nodes have no child nodes, and internal nodes have both parent nodes and child nodes.

Tree Lists
----------

|TreeList| objects are lists of |Tree| objects constrained to sharing the same |TaxonSet|.
Any |Tree| object added to a |TreeList| will have its :attr:`~dendropy.dataobject.tree.Tree.taxon_set` attribute assigned to the |TaxonSet| object of the |TreeList|, and all referenced |Taxon| objects will be mapped to the same or corresponding |Taxon| objects of this new |TaxonSet|, with new |Taxon| objects created if no suitable match is found.

|Tree| and |TreeList| Creation and Reading
==========================================

Creating a New |Tree| or |TreeList| from a Data Source
-------------------------------------------------------

Both the |Tree| and |TreeList| classes support the :meth:`get_from_stream()`, :meth:`get_from_path()`, and :meth:`get_from_string()` factory class methods for simultaneously instantiating and populating objects, taking a data source as the first argument and a data format or schema specification as the second:

    >>> import dendropy
    >>> tree = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
    >>> treelist = dendropy.TreeList.get_from_path('pythonidae.mcmc.nex', 'nexus')

Valid schema specification strings include: "``nexus``", "``newick``", "``nexml``", "``dnafasta``", "``rnafasta``", "``proteinfasta``" etc.

Reading into an Existing |Tree| or |TreeList| from a Data Source
----------------------------------------------------------------

The :meth:`read_from_stream()`, :meth:`read_from_path()`, and :meth:`read_from_string()` instance methods for populating existing objects are also supported, taking the same arguments:

    >>> import dendropy
    >>> tree = dendropy.Tree()
    >>> tree.read_from_path('pythonidae.mle.nex', 'nexus')
    >>> treelist = dendropy.TreeList()
    >>> treelist.read_from_path('pythonidae_cytb.mb.run1.t', 'nexus')
    >>> treelist.read_from_path('pythonidae_cytb.mb.run2.t', 'nexus')
    >>> treelist.read_from_path('pythonidae_cytb.mb.run3.t', 'nexus')
    >>> treelist.read_from_path('pythonidae_cytb.mb.run4.t', 'nexus')

In the case of |Tree| objects, calling :meth:`read_from_*()` *repopulates* (i.e., redefines) the |Tree| with data from the data source, while in the case of |TreeList| objects, calling :meth:`read_from_*()` *appends* the tree definitions in the data source to the |TreeList|.

Cloning an Existing |Tree| or |TreeList|
----------------------------------------

You can also clone existing |Tree| and |TreeList| objects by passing them as arguments to their respective constructors.

For example, to create a clone of a |Tree| object:

    >>> import dendropy
    >>> tree1 = dendropy.Tree.get_from_path('pythonidae.mle.tree', 'nexus')
    >>> tree2 = dendropy.Tree(tree1)

With this, ``tree2`` will be an exact clone of ``tree1``, and can be independentally manipulated (e.g., derooted, branches pruned, splits collapsed, etc.) without effecting ``tree1``.
Note, however, that the |Taxon| objects remain linked: changing the label, for example, of a |Taxon| object on ``tree2`` will result in the label of the corresponding |Taxon| object in ``tree1`` being similarly affected.

To create a clone of a |TreeList| object:

    >>> import dendropy
    >>> treelist1 = dendropy.TreeList.get_from_path('pythonidae.mcmc.nex', 'nexus')
    >>> treelist2 = dendropy.TreeList(treelist1)

Here, ``treelist2`` will be a *deep-copy* of ``treelist1``, i.e., with each |Tree| object in ``treelist2`` being a clone of the corresponding |Tree| object in ``treelist1``.
The same constraint regarding |Taxon| object applies: i.e., the cloning does not extend to |Taxon| objects, and these are shared across all |Tree| objects in both ``treelist1`` and ``treelist2``, as well as the |TreeList| objects themselves.

.. _Customizing_Tree_Creation_and_Reading:

Customizing |Tree| and |TreeList| Creation and Reading
------------------------------------------------------

Using a Specific |TaxonSet|
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Passing a |TaxonSet| object using the ``taxon_set`` argument when instantiating a |Tree| or |TreeList| object (using, for example, the :meth:`get_from_*()` or :meth:`read_from_*()` methods) results in the |Tree| or |TreeList| object being bound to the specified |TaxonSet| object.

Selecting Specific Trees or Subsets of Trees
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With NEXUS and NEWICK data sources, you can also specify ``preserve_underscores=True``.
The NEXUS standard dictates that underscores are equivalent to spaces, and thus any underscore found in any unquoted label in a NEXUS/NEWICK data source will be substituted for spaces.
Specifying ``preserve_underscores=True`` will force DendroPy to keep the underscores.

|Tree| and |TreeList| Saving and Writing
========================================

Writing to Files
----------------

The :meth:`write_to_stream()`, and :meth:`write_to_path()` instance methods allow you to write the data of |Tree| and |TreeList| objects to a file-like object or a file path respectively.
These methods take a file-like object (in the case of :meth:`write_to_stream()`) or a string specifying a filepath (in the case of :meth:`write_to_path()`) as the first argument, and a format or schema specification string as the second argument.

The following example aggregates the post-burn in MCMC samples from a series of NEXUS-formatted files, and saves the collection as a NEWICK-formatted file:

    >>> import dendropy
    >>> treelist = dendropy.TreeList()
    >>> treelist.read_from_path('pythonidae_cytb.mb.run1.t', 'nexus', tree_offset=200)
    >>> treelist.read_from_path('pythonidae_cytb.mb.run2.t', 'nexus', tree_offset=200)
    >>> treelist.read_from_path('pythonidae_cytb.mb.run3.t', 'nexus', tree_offset=200)
    >>> treelist.read_from_path('pythonidae_cytb.mb.run4.t', 'nexus', tree_offset=200)
    >>> treelist.write_to_path('pythonidae_cytb.mcmc-postburnin.tre', 'newick')

Composing a String
------------------

If you do not want to actually write to a file, but instead simply need a string representing the data in a particular format, you can call the instance method :meth:`as_string()`, passing a schema or format specification string as the first argument:

    >>> import dendropy
    >>> tree = dendropy.Tree()
    >>> tree.read_from_path('pythonidae.mle.nex', 'nexus')
    >>> s = tree.as_string('newick')
    >>> print(s)
    >>> (Python_molurus:0.0779719244,((Python_sebae:0.1414715009,(((((Morelia_tracyae:0.0435011998,(Morelia_amethistina:0.0305993564,((Morelia_nauta:0.0092774432,Morelia_kinghorni:0.0093145395):0.005595,Morelia_clastolepis:0.0052046980):0.023435):0.012223):0.025359,Morelia_boeleni:0.0863199106):0.019894,((Python_reticulatus:0.0828549023,Python_timoriensis:0.0963051344):0.072003,Morelia_oenpelliensis:0.0820543043):0.002785):0.002740,((((Morelia_viridis:0.0925974416,(Morelia_carinata:0.0943697342,(Morelia_spilota:0.0237557178,Morelia_bredli:0.0357358071):0.041377):0.005225):0.004424,(Antaresia_maculosa:0.1141193265,((Antaresia_childreni:0.0363195704,Antaresia_stimsoni:0.0188535952):0.043287,Antaresia_perthensis:0.0947695442):0.019148):0.007921):0.022413,(Leiopython_albertisii:0.0698883547,Bothrochilus_boa:0.0811607602):0.020941):0.007439,((Liasis_olivaceus:0.0449896545,(Liasis_mackloti:0.0331564496,Liasis_fuscus:0.0230286886):0.058253):0.016766,Apodora_papuana:0.0847328612):0.008417):0.006539):0.011557,(Aspidites_ramsayi:0.0349772256,Aspidites_melanocephalus:0.0577536309):0.042499):0.036177):0.016859,Python_brongersmai:0.1147218285):0.001271,Python_regius:0.1800489093):0.000000;


Customizing |Tree| and |TreeList| Saving and Writing
-----------------------------------------------------

The following keyword arguments, when passed to :meth:`write_to_stream()`, :meth:`write_to_path()`, or :meth:`as_string()`, allow you to control the formatting of the output:

    ``exclude_taxa``
        When writing NEXUS-formatted data, if :keyword:`True`, then a "``TAXA``" block will not be written. By default, this is :keyword:`False`, i.e., "``TAXA``" blocks will be written.

    ``write_rooting``
        When writing NEXUS-formatted or NEWICK-formatted data, if :keyword:`False`, then tree rooting statements (e.g., "``[&R]``" or "``[&U]``") will not be prefixed to the tree statements. By default, this is :keyword:`True`, i.e., rooting statements will be written.

    ``edge_lengths``
        When writing NEXUS-formatted or NEWICK-formatted data, if :keyword:`False`, then edge or branch lengths will not be written as part of the tree statements. By default, this is :keyword:`True`, i.e., edge lengths will be written.

    ``internal_labels``
        When writing NEXUS-formatted or NEWICK-formatted data, if :keyword:`False`, then labels for internal nodes (if given) will not be written as part of the tree statements. By default, this is :keyword:`True`, i.e., internal node labels will be written.

    ``block_titles``
        When writing NEXUS-formatted data, if :keyword:`False`, then title statements will not be added to the various NEXUS blocks. By default, this is :keyword:`True`, i.e., block titles will be written.

    ``preserve_spaces``
        When writing NEXUS-formatted or NEWICK-formatted data, if :keyword:`True`, then no attempt will be made to produce unquoted labels by substituting spaces for underscores. By default, this is :keyword:`False`, i.e., any label that includes spaces but no other special punctuation character or underscores will have all spaces replaced by underscores so as to allow the label to be represented without quotes.

    ``quote_underscores``
        When writing NEXUS-formatted or NEWICK-formatted data, if :keyword:`False`, then labels will not be wrapped in quotes even if they contain underscores (meaning that the underscores will be interpreted as spaces according to the NEXUS standard). By default, this is :keyword:`True`, meaning that any label that contains underscores will be wrapped in quotes. Note that if a label has any other characters requiring quote protection as specified by the NEXUS standard, then the label will be quoted regardless of the value of this keyword argument.

    ``comment``
        When writing NEXUS-formatted data, then the contents of this variable will be added as NEXUS comment to the output. By default, this is :keyword:`None`.


Taxon Management with Trees and Tree Lists
==========================================

Taxon Management with Trees
---------------------------

It is important to recognize that, by default, DendroPy will create new |TaxonSet| object every time a data source is parsed (and, if the data source has multiple taxon objects, there may be more than one |TaxonSet| created).

Consider the following example::

    >>> import dendropy
    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick')
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick')
    >>> print(t1.description(2))
    Tree object at 0x64b130 (Tree6599856): 7 Nodes, 7 Edges
        [Taxon Set]
            TaxonSet object at 0x64c270 (TaxonSet6603376): 4 Taxa
        [Tree]
            ((A,B),(C,D))
    >>> print(t2.description(2))
    Tree object at 0x64b190 (Tree6600560): 7 Nodes, 7 Edges
        [Taxon Set]
            TaxonSet object at 0x64c1e0 (TaxonSet6603232): 4 Taxa
        [Tree]
            ((A,B),(C,D))

We now have two distinct |Tree| objects, each associated with a distinct |TaxonSet| objects, each with its own set of |Taxon| objects that, while having the same labels, are distinct from one another::

    >>> t1.leaf_nodes()[0].taxon == t2.leaf_nodes()[0].taxon
    False
    >>> t1.leaf_nodes()[0].taxon.label == t2.leaf_nodes()[0].taxon.label
    True

This means that even though the tree shape and structure is identical between the two trees, they exist in different universes as far as DendroPy is concerned, and many operations that involving comparing trees will fail::

    >>> from dendropy import treecalc
    >>> treecalc.robinson_foulds_distance(t1, t2)
    ------------------------------------------------------------
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>

      File "/Users/jeet/Documents/Projects/dendropy/dendropy/treecalc.py", line 263, in robinson_foulds_distance
        value_type=float)

      File "/Users/jeet/Documents/Projects/dendropy/dendropy/treecalc.py", line 200, in splits_distance
        % (hex(id(tree1.taxon_set)), hex(id(tree2.taxon_set))))

    TypeError: Trees have different TaxonSet objects: 0x101f630 vs. 0x103bf30

The solution is to explicitly specify the same ``taxon_set`` when creating the trees. In DendroPy all phylogenetic data classes that are associated with |TaxonSet| objects have constructors, factory methods, and ``read_from_*`` methods take a specific :class:`TaxonSet` object as an argument using the ``taxon_set`` a keyword. For example::

    >>> taxa = dendropy.TaxonSet()
    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick', taxon_set=taxa)
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick', taxon_set=taxa)
    >>> treecalc.robinson_foulds_distance(t1, t2)
    0.0

Taxon Management with Tree Lists
--------------------------------

The |TreeList| class is designed to manage collections of |Tree| objects that share the same |TaxonSet|.
As |Tree| objects are appended to a |TreeList| object, the |TreeList| object will automatically take care of remapping the |TaxonSet| and associated |Taxon| objects::

    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick')
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick')
    >>> print(repr(t1.taxon_set))
    <TaxonSet object at 0x1243a20>
    >>> repr(t1.taxon_set)
    '<TaxonSet object at 0x1243a20>'
    >>> repr(t2.taxon_set)
    '<TaxonSet object at 0x12439f0>'
    >>> trees = dendropy.TreeList()
    >>> trees.append(t1)
    >>> trees.append(t2)
    >>> repr(t1.taxon_set)
    '<TaxonSet object at 0x1243870>'
    >>> repr(t2.taxon_set)
    '<TaxonSet object at 0x1243870>'
    >>> treecalc.robinson_foulds_distance(t1, t2)
    0.0

The same applies when using the :meth:`read_from_*` method of a |TreeList| object: all trees read from the data source will be assigned the same |TaxonSet| object, and the taxa referenced in the tree definition will be mapped to corresponding |Taxon| objects, identified by label, in the |TaxonSet|, with new |Taxon| objects created if no suitable match is found.

While |TreeList| objects ensure that all |Tree| objects created, read or added using them all have the same |TaxonSet| object reference, if two |TreeList| objects are independentally created, they will each have their own, distinct, |TaxonSet| object reference.
For example, if you want to read in two collections of trees and compare trees between the two collections, the following will **not** work:


    >>> import dendropy
    >>> mcmc1 = dendropy.TreeList.get_from_path('pythonidae.mcmc1.nex', 'nexus')
    >>> mcmc2 = dendropy.TreeList.get_from_path('pythonidae.mcmc2.nex', 'nexus')

Of course, reading both data sources into the same  |TreeList| object *will* work insofar as ensuring all the |Tree| objects have the same |TaxonSet|  reference, but then you will lose the distinction between the two sources, unless you keep track of the indexes of where one source begins and the other ends, which error-prone and tedious.
A better approach would be simply to create a |TaxonSet| object, and pass it to the factory methods of both  |TreeList| objects::

    >>> import dendropy
    >>> taxa = dendropy.TaxonSet()
    >>> mcmc1 = dendropy.TreeList.get_from_path('pythonidae.mcmc1.nex', 'nexus', taxon_set=taxa)
    >>> mcmc2 = dendropy.TreeList.get_from_path('pythonidae.mcmc2.nex', 'nexus', taxon_set=taxa)

Now both ``mcmc1`` and ``mcmc2`` share the same |TaxonSet|, and thus so do the |Tree| objects created within them, which means the |Tree| objects can be compared both within and between the collections.

You can also pass the |TaxonSet| to the constructor of |TreeList|.
So, for example, the following is logically identical to the previous::

    >>> import dendropy
    >>> taxa = dendropy.TaxonSet()
    >>> mcmc1 = dendropy.TreeList(taxon_set=taxa)
    >>> mcmc1.read_from_path('pythonidae.mcmc1.nex', 'nexus')
    >>> mcmc2 = dendropy.TreeList(taxon_set=taxa)
    >>> mcmc2.read_from_path('pythonidae.mcmc2.nex', 'nexus')


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

Tree Statistics and Metrics
===========================

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

Probability Under the Coalescent and Counting of Deep Coalescences
------------------------------------------------------------------

The :mod:`~dendropy.coalescent` module provides a range of methods for simulations and calculations under Kingman's coalescent framework and related models:

    :func:`~dendropy.coalescent.log_probability_of_coalescent_tree`
        Given a |Tree| object as the first argument, and the haploid population size as the second, returns the log probability of the |Tree| under the neutral coalescent.

    :func:`~dendropy.coalescent.num_deep_coalescences_with_fitted_tree`
        Given two |Tree| objects, a gene tree and a species tree, sharing the same leaf-set, this returns the number of deep coalescences resulting from fitting the gene tree to the species tree.

    :func:`~dendropy.coalescent.num_deep_coalescences_with_grouping`
        Given a |Tree| object as the first argument, and a list of lists of
        |Taxon| objects representing the expected monophyletic partitioning of the |TaxonSet| of the |Tree| as the second argument, this returns the number of deep coalescences found in the relationships implied by the |Tree| object, conditional on the taxon groupings given by the second argument.

Majority-Rule Consensus Tree from a Collection of Trees
-------------------------------------------------------

To get the majority-rule consensus tree of a |TreeList| object, you can call the :meth:`~dendropy.dataobject.tree.TreeList.consensus()` instance method.
You can specify the frequency threshold for the consensus tree by the ``min_freq`` argument, which default to 0.5 (i.e., a 50% majority rule tree).
The following example aggregates the post-burn-in trees from four MCMC samples into a single |TreeList| object, and prints the 95% majority-rule consensus as a NEWICK string:

.. literalinclude:: /examples/majrule.py
    :linenos:

Frequency of a Split in a Collection of Trees
---------------------------------------------

The :meth:`~dendropy.dataobject.tree.TreeList.frequency_of_split()` method of a |TreeList| object returns the frequency of occurrence of a single split across all the |Tree| objects in the |TreeList|.
The split can be specified by passing a split bitmask directly using the ``split_bitmask`` keyword argument, as a list of |Taxon| objects using the ``taxa`` keyword argument, or as a list of taxon labels using the ``labels`` keyword argument.
The following example shows how to calculate the frequency of a split defined by two taxa, "Morelia amethistina" and "Morelia tracyae", from the post-burn-in trees aggregated across four MCMC samples:

.. literalinclude:: /examples/splitfreq.py
    :linenos:

Tree Distances
--------------

The :mod:`~dendropy.treecalc` module provides a number of functions to calculate the distance between two trees passed as arguments:

    :func:`~dendropy.treecalc.symmetric_distance`
        This function returns the symmetric distance between two trees. The symmetric distance between two trees is the sum of splits found in one of the trees but not the other. It is common to see this statistic called the "Robinson-Foulds distance", but in DendroPy we reserve this term to apply to the Robinson-Foulds distance in the strict sense, i.e., the weighted symmetric distance (see below).

    :func:`~dendropy.treecalc.euclidean_distance`
        This function returns the "branch length distance" of Felsenstein (2004), i.e. the sum of absolute differences in branch lengths for equivalent splits between two trees, with the branch length for a missing split taken to be 0.0.

    :func:`~dendropy.treecalc.robinson_foulds_distance`
        This function returns the Robinsons-Foulds distance between two trees, i.e., the sum of the square of differences in branch lengths for equivalent splits between two trees, with the branch length for a missing split taken to be 0.0.

