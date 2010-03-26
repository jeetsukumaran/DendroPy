******************************
Trees and Collections of Trees
******************************

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

Both the |Tree| and |TreeList| classes support the :meth:`get_from_stream()`, :meth:`get_from_path()`, and :meth:`get_from_string()` factory class methods for simultaneously instantiating and populating objects, taking a data source as the first argument and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` ("``nexus``", "``newick``", "``nexml``", "``fasta``", or "``phylip``", etc.) as the second:

    >>> import dendropy
    >>> tree = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
    >>> treelist = dendropy.TreeList.get_from_path('pythonidae.mcmc.nex', 'nexus')

In addition, fine-grained control over the parsing of the data source is available through various :ref:`keyword arguments <Customizing_Data_Creation_and_Reading>`.


Reading into an Existing |Tree| or |TreeList| from a Data Source
----------------------------------------------------------------

The :meth:`read_from_stream()`, :meth:`read_from_path()`, and :meth:`read_from_string()` instance methods for populating existing objects are also supported, taking the same arguments (i.e., a data source, a :ref:`schema specification string <Specifying_the_Data_Source_Format>`, as well as optional :ref:`keyword arguments <Customizing_Data_Creation_and_Reading>` to customize the parse behavior):

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


|Tree| and |TreeList| Saving and Writing
========================================

Writing to Files
----------------

The :meth:`write_to_stream()`, and :meth:`write_to_path()` instance methods allow you to write the data of |Tree| and |TreeList| objects to a file-like object or a file path respectively.
These methods take a file-like object (in the case of :meth:`write_to_stream()`) or a string specifying a filepath (in the case of :meth:`write_to_path()`) as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Writing_Format>` as the second argument.

The following example aggregates the post-burn in MCMC samples from a series of NEXUS-formatted files, and saves the collection as a Newick-formatted file:

    >>> import dendropy
    >>> treelist = dendropy.TreeList()
    >>> treelist.read_from_path('pythonidae_cytb.mb.run1.t', 'nexus', tree_offset=200)
    >>> treelist.read_from_path('pythonidae_cytb.mb.run2.t', 'nexus', tree_offset=200)
    >>> treelist.read_from_path('pythonidae_cytb.mb.run3.t', 'nexus', tree_offset=200)
    >>> treelist.read_from_path('pythonidae_cytb.mb.run4.t', 'nexus', tree_offset=200)
    >>> treelist.write_to_path('pythonidae_cytb.mcmc-postburnin.tre', 'newick')

Fine-grained control over the output format can be specified using :ref:`keyword arguments <Customizing_the_Data_Writing_Format>`.

Composing a String
------------------

If you do not want to actually write to a file, but instead simply need a string representing the data in a particular format, you can call the instance method :meth:`as_string()`, passing a :ref:`schema specification string <Specifying_the_Data_Writing_Format>` as the first argument:

    >>> import dendropy
    >>> tree = dendropy.Tree()
    >>> tree.read_from_path('pythonidae.mle.nex', 'nexus')
    >>> s = tree.as_string('newick')
    >>> print(s)
    >>> (Python_molurus:0.0779719244,((Python_sebae:0.1414715009,(((((Morelia_tracyae:0.0435011998,(Morelia_amethistina:0.0305993564,((Morelia_nauta:0.0092774432,Morelia_kinghorni:0.0093145395):0.005595,Morelia_clastolepis:0.0052046980):0.023435):0.012223):0.025359,Morelia_boeleni:0.0863199106):0.019894,((Python_reticulatus:0.0828549023,Python_timoriensis:0.0963051344):0.072003,Morelia_oenpelliensis:0.0820543043):0.002785):0.002740,((((Morelia_viridis:0.0925974416,(Morelia_carinata:0.0943697342,(Morelia_spilota:0.0237557178,Morelia_bredli:0.0357358071):0.041377):0.005225):0.004424,(Antaresia_maculosa:0.1141193265,((Antaresia_childreni:0.0363195704,Antaresia_stimsoni:0.0188535952):0.043287,Antaresia_perthensis:0.0947695442):0.019148):0.007921):0.022413,(Leiopython_albertisii:0.0698883547,Bothrochilus_boa:0.0811607602):0.020941):0.007439,((Liasis_olivaceus:0.0449896545,(Liasis_mackloti:0.0331564496,Liasis_fuscus:0.0230286886):0.058253):0.016766,Apodora_papuana:0.0847328612):0.008417):0.006539):0.011557,(Aspidites_ramsayi:0.0349772256,Aspidites_melanocephalus:0.0577536309):0.042499):0.036177):0.016859,Python_brongersmai:0.1147218285):0.001271,Python_regius:0.1800489093):0.000000;

As above, fine-grained control over the output format can be specified using :ref:`keyword arguments <Customizing_the_Data_Writing_Format>`.

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

Viewing and Displaying Trees
============================

Sometimes it is useful to get a visual representation of a |Tree|.

For quick inspection, the :meth:`~dendropy.dataobject.tree.Tree.print_plot()` will write an ASCII text plot to the standard output stream::

    >>> t = dendropy.Tree.get_from_string("(A,(B,(C,D)))", "newick")
    >>> t.print_plot()
    /----------------------------------------------- A
    +
    |                /------------------------------ B
    \----------------+
                     |          /------------------- C
                     \----------+
                                \------------------- D

If you need to store this representation as a string instead, you can use :meth:`~dendropy.dataobject.tree.Tree.as_ascii_plot()`::

    >>> s = t.as_ascii_plot()
    >>> print(s)
    /----------------------------------------------- A
    +
    |                /------------------------------ B
    \----------------+
                     |          /------------------- C
                     \----------+
                                \------------------- D

While the :meth:`~dendropy.dataobject.tree.Tree.write_to_path()`, :meth:`~dendropy.dataobject.tree.Tree.write_to_stream()` and :meth:`~dendropy.dataobject.tree.Tree.as_string()` methods provide for a rich and flexible way to write representations of a |Tree| in various formats to various destinations, the :meth:`~dendropy.dataobject.tree.Tree.print_newick()` provides a quick-and-dirty way to get a snapshot NEWICK string of the tree::

    >>> t.print_newick()
    (A,(B,(C,D)))





