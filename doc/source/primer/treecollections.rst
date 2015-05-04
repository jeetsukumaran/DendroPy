********************
Collections of Trees
********************

Collections of Trees: The |TreeList| Class
==========================================

|TreeList| objects are collections of |Tree| objects constrained to sharing the same |TaxonNamespace|.
Any |Tree| object added to a |TreeList| will have its :attr:`~dendropy.datamodel.treemodel.Tree.taxon_namespace` attribute assigned to the |TaxonNamespace| object of the |TreeList|, and all referenced |Taxon| objects will be mapped to the same or corresponding |Taxon| objects of this new |TaxonNamespace|, with new |Taxon| objects created if no suitable match is found.

Reading and Writing |TreeList| Instances
========================================

The |TreeList| class supports the ":meth:`~dendropy.datamodel.treecollectionmodel.TreeList.get`" factory class method for simultaneously instantiating and populating |TreeList| instances, taking a data source as the first argument and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` ("``nexus``", "``newick``", "``nexml``", "``fasta``", or "``phylip``", etc.) as the second::

    import dendropy
    treelist = dendropy.TreeList.get(path='pythonidae.mcmc.nex', schema='nexus')

The ":meth:`~dendropy.datamodel.treecollectionmodel.TreeList.read`" instance method can be used to add trees from a data source to an existing |TreeList| instance:

.. literalinclude:: /examples/tree_list_add1.py

A |TreeList| object can be written to an external resource using the ":meth:`~dendropy.datamodel.treecollectionmodel.TreeList.write`" method::

    import dendropy
    treelist = dendropy.TreeList.get(
        path="trees1.nex",
        schema="nexus",
        )
    treelist.write(
        path="trees1.newick",
        schema="newick",
        )


It can also be represented as a string using the ":meth:`~dendropy.datamodel.treecollectionmodel.TreeList.as_string`" method::

    import dendropy
    treelist = dendropy.TreeList.get(
        path="trees1.nex",
        schema="nexus",
        )
    print(treelist.as_string(schema="newick",)

More information on reading operations is available in the :doc:`/primer/reading_and_writing` section.

Using and Managing the Collections of Trees
===========================================

A |TreeList| behaves very much like a list, supporting iteration, indexing, slices, removal, indexing, sorting, etc.:

.. literalinclude:: /examples/tree_list_ops1.py


The |TreeList| class supports the native Python ``list`` interface methods of adding individual |Tree| instances through
        :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.append`,
        :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.extend`,
        :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.insert`,
        and other methods, but with the added aspect of :doc:`taxon namespace </primer/taxa>` migration:

.. literalinclude:: /examples/tree_list_ops2.py


Cloning/Copying a |TreeList|
============================

You can make a *shallow*-copy of a |TreeList| calling :meth:`dendropy.datamodel.treecollectionmodel.TreeList.clone` with a "``depth``" argument value of 0 or by slicing:

.. literalinclude:: /examples/tree_list_copy1.py

With a shallow-copy, the actual |Tree| instances are shared between lists (as is the |TaxonNamespace|).

For a taxon namespace-scoped *deep*-copy, on the other hand, i.e., where the |Tree| instances are also cloned but the |Taxon| and |TaxonNamespace| references are preserved, you can call :meth:`dendropy.datamodel.treecollectionmodel.TreeList.clone` with a "``depth``" argument value of 1 or by copy construction:

.. literalinclude:: /examples/tree_list_copy2.py

Finally, for a true and complete deep-copy, where even the |Taxon| and |TaxonNamespace| references are copied, call :func:`copy.deepcopy`:

.. literalinclude:: /examples/tree_list_copy3.py

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

Fine-grained control over the output format can be specified using :doc:`keyword arguments </primer/reading_and_writing>`.

Composing a String
------------------

If you do not want to actually write to a file, but instead simply need a string representing the data in a particular format, you can call the instance method :meth:`as_string()`, passing a :ref:`schema specification string <Specifying_the_Data_Writing_Format>` as the first argument:

    >>> import dendropy
    >>> tree = dendropy.Tree()
    >>> tree.read_from_path('pythonidae.mle.nex', 'nexus')
    >>> s = tree.as_string('newick')
    >>> print(s)
    >>> (Python_molurus:0.0779719244,((Python_sebae:0.1414715009,(((((Morelia_tracyae:0.0435011998,(Morelia_amethistina:0.0305993564,((Morelia_nauta:0.0092774432,Morelia_kinghorni:0.0093145395):0.005595,Morelia_clastolepis:0.0052046980):0.023435):0.012223):0.025359,Morelia_boeleni:0.0863199106):0.019894,((Python_reticulatus:0.0828549023,Python_timoriensis:0.0963051344):0.072003,Morelia_oenpelliensis:0.0820543043):0.002785):0.002740,((((Morelia_viridis:0.0925974416,(Morelia_carinata:0.0943697342,(Morelia_spilota:0.0237557178,Morelia_bredli:0.0357358071):0.041377):0.005225):0.004424,(Antaresia_maculosa:0.1141193265,((Antaresia_childreni:0.0363195704,Antaresia_stimsoni:0.0188535952):0.043287,Antaresia_perthensis:0.0947695442):0.019148):0.007921):0.022413,(Leiopython_albertisii:0.0698883547,Bothrochilus_boa:0.0811607602):0.020941):0.007439,((Liasis_olivaceus:0.0449896545,(Liasis_mackloti:0.0331564496,Liasis_fuscus:0.0230286886):0.058253):0.016766,Apodora_papuana:0.0847328612):0.008417):0.006539):0.011557,(Aspidites_ramsayi:0.0349772256,Aspidites_melanocephalus:0.0577536309):0.042499):0.036177):0.016859,Python_brongersmai:0.1147218285):0.001271,Python_regius:0.1800489093):0.000000;

As above, fine-grained control over the output format can be specified using :doc:`keyword arguments </primer/reading_and_writing>`.

Taxon Management with Trees and Tree Lists
==========================================

Taxon Management with Trees
---------------------------

It is important to recognize that, by default, DendroPy will create new |TaxonNamespace| object every time a data source is parsed (and, if the data source has multiple taxon objects, there may be more than one |TaxonNamespace| created).

Consider the following example::

    >>> import dendropy
    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick')
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick')
    >>> print(t1.description(2))
    Tree object at 0x64b130 (Tree6599856): 7 Nodes, 7 Edges
        [Taxon Set]
            TaxonNamespace object at 0x64c270 (TaxonNamespace6603376): 4 Taxa
        [Tree]
            ((A,B),(C,D))
    >>> print(t2.description(2))
    Tree object at 0x64b190 (Tree6600560): 7 Nodes, 7 Edges
        [Taxon Set]
            TaxonNamespace object at 0x64c1e0 (TaxonNamespace6603232): 4 Taxa
        [Tree]
            ((A,B),(C,D))

We now have two distinct |Tree| objects, each associated with a distinct |TaxonNamespace| objects, each with its own set of |Taxon| objects that, while having the same labels, are distinct from one another::

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
        % (hex(id(tree1.taxon_namespace)), hex(id(tree2.taxon_namespace))))

    TypeError: Trees have different TaxonNamespace objects: 0x101f630 vs. 0x103bf30

The solution is to explicitly specify the same ``taxon_namespace`` when creating the trees. In DendroPy all phylogenetic data classes that are associated with |TaxonNamespace| objects have constructors, factory methods, and "|read_from_methods|" methods take a specific :class:`TaxonNamespace` object as an argument using the ``taxon_namespace`` a keyword. For example::

    >>> taxa = dendropy.TaxonNamespace()
    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick', taxon_namespace=taxa)
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick', taxon_namespace=taxa)
    >>> treecalc.robinson_foulds_distance(t1, t2)
    0.0

Taxon Management with Tree Lists
--------------------------------

The |TreeList| class is designed to manage collections of |Tree| objects that share the same |TaxonNamespace|.
As |Tree| objects are appended to a |TreeList| object, the |TreeList| object will automatically take care of remapping the |TaxonNamespace| and associated |Taxon| objects::

    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick')
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', schema='newick')
    >>> print(repr(t1.taxon_namespace))
    <TaxonNamespace object at 0x1243a20>
    >>> repr(t1.taxon_namespace)
    '<TaxonNamespace object at 0x1243a20>'
    >>> repr(t2.taxon_namespace)
    '<TaxonNamespace object at 0x12439f0>'
    >>> trees = dendropy.TreeList()
    >>> trees.append(t1)
    >>> trees.append(t2)
    >>> repr(t1.taxon_namespace)
    '<TaxonNamespace object at 0x1243870>'
    >>> repr(t2.taxon_namespace)
    '<TaxonNamespace object at 0x1243870>'
    >>> treecalc.robinson_foulds_distance(t1, t2)
    0.0

The same applies when using the "|read_from_methods|" method of a |TreeList| object: all trees read from the data source will be assigned the same |TaxonNamespace| object, and the taxa referenced in the tree definition will be mapped to corresponding |Taxon| objects, identified by label, in the |TaxonNamespace|, with new |Taxon| objects created if no suitable match is found.

While |TreeList| objects ensure that all |Tree| objects created, read or added using them all have the same |TaxonNamespace| object reference, if two |TreeList| objects are independentally created, they will each have their own, distinct, |TaxonNamespace| object reference.
For example, if you want to read in two collections of trees and compare trees between the two collections, the following will **not** work:


    >>> import dendropy
    >>> mcmc1 = dendropy.TreeList.get_from_path('pythonidae.mcmc1.nex', 'nexus')
    >>> mcmc2 = dendropy.TreeList.get_from_path('pythonidae.mcmc2.nex', 'nexus')

Of course, reading both data sources into the same  |TreeList| object *will* work insofar as ensuring all the |Tree| objects have the same |TaxonNamespace|  reference, but then you will lose the distinction between the two sources, unless you keep track of the indexes of where one source begins and the other ends, which error-prone and tedious.
A better approach would be simply to create a |TaxonNamespace| object, and pass it to the factory methods of both  |TreeList| objects::

    >>> import dendropy
    >>> taxa = dendropy.TaxonNamespace()
    >>> mcmc1 = dendropy.TreeList.get_from_path('pythonidae.mcmc1.nex', 'nexus', taxon_namespace=taxa)
    >>> mcmc2 = dendropy.TreeList.get_from_path('pythonidae.mcmc2.nex', 'nexus', taxon_namespace=taxa)

Now both ``mcmc1`` and ``mcmc2`` share the same |TaxonNamespace|, and thus so do the |Tree| objects created within them, which means the |Tree| objects can be compared both within and between the collections.

You can also pass the |TaxonNamespace| to the constructor of |TreeList|.
So, for example, the following is logically identical to the previous::

    >>> import dendropy
    >>> taxa = dendropy.TaxonNamespace()
    >>> mcmc1 = dendropy.TreeList(taxon_namespace=taxa)
    >>> mcmc1.read_from_path('pythonidae.mcmc1.nex', 'nexus')
    >>> mcmc2 = dendropy.TreeList(taxon_namespace=taxa)
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
Additional keyword arguments to customize the parsing are the same as that for the general "|get_from_methods|" and "|read_from_methods|" methods.
For example, the following script reads a model tree from a file, and then iterates over a collection of MCMC trees in another file, calculating a storing the symmetric distance between the model tree and each of the MCMC trees one at time:

.. literalinclude:: /examples/tree_iter1.py
    :linenos:

Note how a |TaxonNamespace| object is created and passed to both the :meth:`~dendropy.datamodel.treemodel.Tree.get_from_path()` and the :func:`~dendropy.dataio.tree_source_iter()` functions using the ``taxon_namespace`` keyword argument.
This is to ensure that the corresponding taxa in both sources get mapped to the same |Taxon| objects in DendroPy object space, so as to enable comparisons of the trees.
If this was not done, then each tree would have its own distinct |TaxonNamespace| object (and associated |Taxon| objects), making comparisons impossible.

Also note how the ``tree_offset`` keyword is used to skip over the burn-in trees from the MCMC sample.

If you want to iterate over trees in multiple sources, you can use the :func:`~dendropy.dataio.multi_tree_source_iter()`.
This takes a list of file-like objects *or* a list of filepath strings as its first argument, and a schema-specification string as its second argument.
Again, other keyword arguments supported by the general "|get_from_methods|" and "|read_from_methods|" methods are also available.

For example:

.. literalinclude:: /examples/tree_iter2.py
    :linenos:

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

While the :meth:`~dendropy.datamodel.treemodel.Tree.write_to_path()`, :meth:`~dendropy.datamodel.treemodel.Tree.write_to_stream()` and :meth:`~dendropy.datamodel.treemodel.Tree.as_string()` methods provide for a rich and flexible way to write representations of a |Tree| in various formats to various destinations, the :meth:`~dendropy.datamodel.treemodel.Tree.print_newick()` provides a quick-and-dirty way to get a snapshot NEWICK string of the tree::

    >>> t.print_newick()
    (A,(B,(C,D)))





