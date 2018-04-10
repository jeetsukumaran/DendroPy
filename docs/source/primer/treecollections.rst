********************
Collections of Trees
********************

Collections of Trees: The |TreeList| Class
==========================================

|TreeList| objects are collections of |Tree| objects constrained to sharing the same |TaxonNamespace|.
Any |Tree| object added to a |TreeList| will have its :attr:`~dendropy.datamodel.treemodel.Tree.taxon_namespace` attribute assigned to the |TaxonNamespace| object of the |TreeList|, and all referenced |Taxon| objects will be mapped to the same or corresponding |Taxon| objects of this new |TaxonNamespace|, with new |Taxon| objects created if no suitable match is found.
Objects of the |TreeList| class have an "``annotations``" attribute, which is a :class:`~dendropy.datamodel.basemodel.AnnotationSet` object, i.e. a collection of :class:`~dendropy.datamodel.basemodel.Annotation` instances tracking metadata.
More information on working with metadata can be found in the ":doc:`/primer/working_with_metadata_annotations`" section.

Reading and Writing |TreeList| Instances
----------------------------------------

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
-------------------------------------------

A |TreeList| behaves very much like a list, supporting iteration, indexing, slices, removal, indexing, sorting, etc.:

.. literalinclude:: /examples/tree_list_ops1.py


The |TreeList| class supports the native Python ``list`` interface methods of adding individual |Tree| instances through
        :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.append`,
        :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.extend`,
        :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.insert`,
        and other methods, but with the added aspect of :doc:`taxon namespace </primer/taxa>` migration:

.. literalinclude:: /examples/tree_list_ops2.py


Cloning/Copying a |TreeList|
----------------------------

You can make a *shallow*-copy of a |TreeList| calling :meth:`dendropy.datamodel.treecollectionmodel.TreeList.clone` with a "``depth``" argument value of 0 or by slicing:

.. literalinclude:: /examples/tree_list_copy1.py

With a shallow-copy, the actual |Tree| instances are shared between lists (as is the |TaxonNamespace|).

For a taxon namespace-scoped *deep*-copy, on the other hand, i.e., where the |Tree| instances are also cloned but the |Taxon| and |TaxonNamespace| references are preserved, you can call :meth:`dendropy.datamodel.treecollectionmodel.TreeList.clone` with a "``depth``" argument value of 1 or by copy construction:

.. literalinclude:: /examples/tree_list_copy2.py

Finally, for a true and complete deep-copy, where even the |Taxon| and |TaxonNamespace| references are copied, call :func:`copy.deepcopy`:

.. literalinclude:: /examples/tree_list_copy3.py

Efficiently Iterating Over Trees in a File
==========================================

If you need to process a collection of trees defined in a file source, you can, of course, read the trees into a |TreeList| object and iterate over the resulting collection::

    import dendropy
    trees = dendropy.TreeList.get(path='pythonidae.beast-mcmc.trees', schema='nexus')
    for tree in trees:
        print(tree.as_string('newick'))

In the above, the entire data source is parsed and stored in the ``trees`` object before being processed in the subsequent lines.
In some cases, you might not need to maintain all the trees in memory at the same time.
For example, you might be interested in calculating the distribution of a statistic over a collection of trees, but have no need to refer to any of the trees after the statistic has been calculated.
In this case, it will be more efficient to use the :meth:`~dendropy.datamodel.treemodel.Tree.yield_from_files` function.
This takes a *list* or any other iterable of file-like objects or strings (giving filepaths) as the first argument ("``files``") and a mandatory :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second argument ("``schema``).
Additional keyword arguments to customize the parsing are the same as that for the general "|get|" and "|read|" methods.
For example, the following script reads a model tree from a file, and then iterates over a collection of MCMC trees in a set of files, calculating and storing the symmetric distance between the model tree and each of the MCMC trees one at time:

.. literalinclude:: /examples/tree_iter1.py

Note how a |TaxonNamespace| object is created and passed to both the :meth:`~dendropy.datamodel.treemodel.Tree.get` and the :meth:`~dendropy.datamodel.treemodel.Tree.yield_from_files` functions using the ``taxon_namespace`` keyword argument.
This is to ensure that the corresponding taxa in both sources get mapped to the same |Taxon| objects in DendroPy object space, so as to enable comparisons of the trees.
If this was not done, then each tree would have its own distinct |TaxonNamespace| object (and associated |Taxon| objects), making comparisons impossible.

When the number of trees are large or the trees themselves are large or both, iterating over trees in files using :meth:`~dendropy.datamodel.treemodel.Tree.yield_from_files` is almost always going to give the best performance, sometimes orders of magnitude faster.
This is due to avoiding the Python virtual machine itself from slowing down due to memory usage.
