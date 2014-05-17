*********
Data Sets
*********

The |DataSet| class provides for objects that allow you to manage multiple types of phylogenetic data.

It has three primary attributes:

    :attr:`~dendropy.dataobject.dataset.DataSet.taxon_namespaces`
        A list of all |TaxonNamespace| objects in the |DataSet|, in the order that they were added or read, include |TaxonNamespace| objects added implicitly through being associated with added |TreeList| or |CharacterMatrix| objects.

    :attr:`~dendropy.dataobject.dataset.DataSet.tree_lists`
        A list of all |TreeList| objects in the |DataSet|, in the order that they were added or read.

    :attr:`~dendropy.dataobject.dataset.DataSet.char_matrices`
        A list of all |CharacterMatrix| objects in the |DataSet|, in the order that they were added or read.

|DataSet| Creation and Reading
===============================

Creating a new |DataSet| from a Data Source
--------------------------------------------

You can use the :meth:`get_from_stream()`, :meth:`get_from_path()`, and :meth:`get_from_string()` factory class methods for simultaneously instantiating and populating an object, taking a data source as the first argument and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` ("``nexus``", "``newick``", "``nexml``", "``fasta``", "``phylip``", etc.) as the second:

    >>> import dendropy
    >>> ds = dendropy.DataSet.get_from_path('pythonidae.nex', 'nexus')

In addition, fine-grained control over the parsing of the data source is available through various :ref:`keyword arguments <Customizing_Data_Creation_and_Reading>`.
Reading into an Existing |DataSet| from a Data Source
-----------------------------------------------------

The :meth:`read_from_stream()`, :meth:`read_from_path()`, and :meth:`read_from_string()` instance methods for populating existing objects are also supported, taking the same arguments (i.e., a data source, a :ref:`schema specification string <Specifying_the_Data_Source_Format>`, as well as optional :ref:`keyword arguments <Customizing_Data_Creation_and_Reading>` to customize the parse behavior)

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.attach_taxon_namespace()
    >>> ds = dendropy.DataSet.read_from_path('pythonidae.cytb.fasta', 'dnafasta')
    >>> ds = dendropy.DataSet.read_from_path('pythonidae.mle.nex', 'nexus')

Note how the :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_namespace()` method is called before invoking any :meth:`read_from_*()` statements, to ensure that all the taxon references in the data sources get mapped to the same |TaxonNamespace| instance.

Cloning an Existing |DataSet|
-----------------------------

You can also clone an existing |DataSet| object by passing it as an argument to the |DataSet| constructor:

    >>> import dendropy
    >>> ds1 = dendropy.DataSet.get_from_path('pythonidae.cytb.fasta', 'dnafasta')
    >>> ds2 = dendropy.DataSet(ds1)

Following this, ``ds2`` will be a *full* deep-copy clone of ``ds1``, with distinct and independent, but identical, |Taxon|, |TaxonNamespace|, |TreeList|, |Tree| and |CharacterMatrix| objects.
Note that, in distinction to the similar cloning methods of |Tree| and |TreeList|, even the |Taxon| and |TaxonNamespace| objects are cloned, meaning that you manipulate the |Taxon| and |TaxonNamespace| objects of ``ds2`` without in any way effecting those of ``ds1``.

Creating a New |DataSet| from Existing |TreeList| and |CharacterMatrix| Objects
-------------------------------------------------------------------------------

You can add independentally created or parsed data objects to a |DataSet| by passing them as unnamed arguments to the constructor:

    >>> import dendropy
    >>> treelist1 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run1.t', 'nexus')
    >>> treelist2 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run2.t', 'nexus')
    >>> treelist3 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run3.t', 'nexus')
    >>> treelist4 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run4.t', 'nexus')
    >>> cytb = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.fasta', 'dnafasta')
    >>> ds = dendropy.DataSet(cytb, treelist1, treelist2, treelist3, treelist4)
    >>> ds.unify_taxa()

Note how we call the instance method :meth:`~dendropy.dataobject.dataset.DataSet.unify_taxa()` after the creation of the |DataSet| object.
This method will remove all existing |TaxonNamespace| objects from the |DataSet|, create and add a new one, and then map all taxon references in all contained |TreeList| and |CharacterMatrix| objects to this new, unified |TaxonNamespace|.

Adding Data to an Exisiting |DataSet|
-------------------------------------

You can add independentally created or parsed data objects to a |DataSet| using the :meth:`~dendropy.dataobject.dataset.DataSet.add()` method:

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> treelist1 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run1.t', 'nexus')
    >>> treelist2 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run2.t', 'nexus')
    >>> treelist3 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run3.t', 'nexus')
    >>> treelist4 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run4.t', 'nexus')
    >>> cytb = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.fasta', 'dnafasta')
    >>> ds.add(treelist1)
    >>> ds.add(treelist2)
    >>> ds.add(treelist3)
    >>> ds.add(treelist4)
    >>> ds.add(cytb)
    >>> ds.unify_taxa()

Here, again, we call the :meth:`~dendropy.dataobject.dataset.DataSet.unify_taxa()` to map all taxon references to the same, common, unified |TaxonNamespace|.

|DataSet| Saving and Writing
=============================

Writing to Files
----------------

The :meth:`write_to_stream()`, and :meth:`write_to_path()` instance methods allow you to write the data of a |DataSet| object to a file-like object or a file path respectively.
These methods take a file-like object (in the case of :meth:`write_to_stream()`) or a string specifying a filepath (in the case of :meth:`write_to_path()`) as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Writing_Format>` as the second argument.

The following example aggregates the post-burn in MCMC samples from a series of NEXUS-formatted tree files into a single |TreeList|, then, adds the |TreeList| as well as the original character data into a single |DataSet| object, which is then written out as NEXUS-formatted file:

.. literalinclude:: /examples/dsrw1.py
    :linenos:

Fine-grained control over the output format can be specified using :ref:`keyword arguments <Customizing_the_Data_Writing_Format>`.

Composing a String
------------------

If you do not want to actually write to a file, but instead simply need a string representing the data in a particular format, you can call the instance method :meth:`as_string()`, passing a :ref:`schema specification string <Specifying_the_Data_Writing_Format>` as the first argument:

    >>> import dendropy
    >>> ds = dendropy.DataSet(attached_taxon_namespace=True)
    >>> ds.read_from_path('pythonidae.cytb.fasta', 'dnafasta')
    >>> s = ds.as_string('nexus')

As above, fine-grained control over the output format can be specified using :ref:`keyword arguments <Customizing_the_Data_Writing_Format>`.

Taxon Management with Data Sets
===============================

The |DataSet| object, representing a meta-collection of phylogenetic data, differs in one important way from all the other phylogenetic data objects discussed so far with respect to taxon management, in that it is not associated with any particular |TaxonNamespace| object.
Rather, it maintains a list (in the property :attr:`~dendropy.dataobject.char.DataSet.taxon_namespaces`) of *all* the |TaxonNamespace| objects referenced by its contained |TreeList| objects (in the property :attr:`~dendropy.dataobject.char.DataSet.tree_lists`) and |CharacterMatrix| objects (in the property :attr:`~dendropy.dataobject.char.DataSet.char_matrices`).

With respect to taxon management, |DataSet| objects operate in one of two modes: "detached taxon set" mode and "attached taxon set" mode.

Detached (Multiple) Taxon Set Mode
----------------------------------

In the "detached taxon set" mode, which is the default, |DataSet| object tracks all |TaxonNamespace| references of their other data members in the property :attr:`~dendropy.dataobject.char.DataSet.taxon_namespaces`, but no effort is made at taxon management as such.
Thus, every time a data source is read with a "detached taxon set" mode |DataSet| object, by deault, a new  |TaxonNamespace| object will be created and associated with the |Tree|, |TreeList|, or |CharacterMatrix| objects created from each data source, resulting in multiple |TaxonNamespace| independent references.
As such, "detached taxon set" mode |DataSet| objects are suitable for handling data with multiple distinct sets of taxa.

For example::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read_from_path("primates.nex", "nexus")
    >>> ds.read_from_path("snakes.nex", "nexus")

The dataset, ``ds``, will now contain two distinct sets of |TaxonNamespace| objects, one for the taxa defined in "primates.nex", and the other for the taxa defined for "snakes.nex".
In this case, this behavior is correct, as the two files do indeed refer to different sets of taxa.

However, consider the following::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")

Here, even though all the data files refer to the same set of taxa, the resulting  |DataSet| object will actually have 4 distinct  |TaxonNamespace| objects, one for each of the independent reads, and a taxon with a particular label in the first file (e.g., "Python regius" of "pythonidae_cytb.fasta") will map to a completely distinct |Taxon| object than a taxon with the same label in the second file (e.g., "Python regius" of "pythonidae_aa.nex").
This is incorrect behavior, and to achieve the correct behavior with a multiple taxon set mode |DataSet| object, we need to explicitly pass a |TaxonNamespace| object to each of the :meth:`~dendropy.dataobject.dataset.DataSet.read_from_path()` statements::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus", taxon_namespace=ds.taxon_namespaces[0])
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus", taxon_namespace=ds.taxon_namespaces[0])
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus", taxon_namespace=ds.taxon_namespaces[0])
    >>> ds.write_to_path("pythonidae_combined.nex", "nexus")

In the previous example, the first :meth:`~dendropy.dataobject.dataset.DataSet.read_from_path()` statement results in a new |TaxonNamespace| object, which is added to the :attr:`~dendropy.dataobject.char.DataSet.taxon_namespaces` property of the |DataSet| object ``ds``.
This |TaxonNamespace| object gets passed via the ``taxon_namespace`` keyword to subsequent :meth:`~dendropy.dataobject.dataset.DataSet.read_from_path()` statements, and thus as each of the data sources are processed, the taxon references get mapped to |Taxon| objects in the same, single, |TaxonNamespace| object.

While this approach works to ensure correct taxon mapping across multiple data object reads and instantiation, in this context, it is probably more convenient to use the |DataSet| in "attached taxon set" mode.

Attached (Single) Taxon Set Mode
--------------------------------
In the "attached taxon set" mode, |DataSet| objects ensure that the taxon references of all data objects that are added to them are mapped to the same |TaxonNamespace| object (at least one for each independent read or creation operation).
The "attached taxon set" mode can be set by passing the keyword argument ``attach_taxon_namespace=True`` to the constructor of the |DataSet| when instantiating a new |DataSet| object (in which case a new |TaxonNamespace| object will be created and added to the |DataSet| object as the default), by passing an existing |TaxonNamespace| object to which to attach using the keyword argument ``taxon_namespace``, or by calling :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_namespace()` on an existing |DataSet| object

For example::

    >>> import dendropy
    >>> ds = dendropy.DataSet(attach_taxon_namespace=True)
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")

Or::

    >>> import dendropy
    >>> taxa = dendropy.TaxonNamespace(label="global")
    >>> ds = dendropy.DataSet(taxon_namespace=taxa)
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")

Or::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.attach_taxon_namespace()
    <TaxonNamespace object at 0x5779c0>
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")

All of the above will result in only a single |TaxonNamespace| object that have all the taxa from the four data sources mapped to them.
Note how :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_namespace()` returns the new |TaxonNamespace| object created and attached when called.
If you needed to detach the |TaxonNamespace| object and then later on reattach it again, you would assign the return value to a variable, and pass it to as an argument to the later call to :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_namespace()`.

Switching Between Attached and Detached Taxon Set Modes
-------------------------------------------------------
As noted above, you can use the :meth:`~dendropy.dataobject.dataset.DataSet.attached_taxon_namespace()` method to switch a |DataSet| object to attached taxon set mode.
To restore it to multiple taxon set mode, you would use the :meth:`~dendropy.dataobject.dataset.DataSet.detach_taxon_namespace()` method::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.attach_taxon_namespace()
    <TaxonNamespace object at 0x5779c0>
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")
    >>> ds.detach_taxon_namespace()
    >>> ds.read_from_path("primates.nex", "nexus")

Here, the same |TaxonNamespace| object is used to manage taxon references for data parsed from the first four files, while the data from the fifth and final file gets its own, distinct, |TaxonNamespace| object and associated |Taxon| object references.

Attaching a Particular Taxon Set
--------------------------------

When :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_namespace()` is called without arguments, a new |TaxonNamespace| object is created and added to the :attr:`~dendropy.dataobject.char.DataSet.taxon_namespaces` list of the |DataSet| object, and taxon references of all data subsequently read (or created and added independentally) will be mapped to |Taxon| objects in this new |TaxonNamespace| object.
If you want to use an existing |TaxonNamespace| object instead of a new one, you can pass this object as an argument to the :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_namespace()` method::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("primates.nex", "nexus")
    >>> ds.attach_taxon_namespace(ds.taxon_namespaces[0])
    <TaxonNamespace object at 0x5b8150>
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")
    >>> ds.detach_taxon_namespace()

Here, the first two :meth:`~dendropy.dataobject.dataset.DataSet.read_from_path()` statements result in two distinct |TaxonNamespace| objects, one for each read, each with their own independent |Taxon| objects.
The :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_namespace()` statement is passed the |TaxonNamespace| object from the first read operation, and all data created from the next three :meth:`~dendropy.dataobject.dataset.DataSet.read_from_path()` statements will have their taxon references mapped to this first |TaxonNamespace| object.


.. SCRATCH
    Unattached vs. Attached |TaxonNamespace| Modes
    ========================================

    A |DataSet| object can manage taxon references in one of two modes: unattached or attached.
    The "unattached" taxon set mode is the default.
    In this mode, every time a data source is parsed, at least one new |TaxonNamespace| object will be created to manage taxon references in the data source, and new |Taxon| objects will be created and added to this |TaxonNamespace| for every taxon reference in the data source.
    This means that multiple read statements will result in multiple |TaxonNamespace| objects being created and added to the |DataSet|.
    In contrast, in "attached" taxon set mode, a single |TaxonNamespace| object will be used to manage taxon references across all data source reading operations.

    A |DataSet| can be placed in attached taxon set mode by calling the instance method :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_namespace()`.
    This method optionally takes a |TaxonNamespace| object as an argument that will be used as the |TaxonNamespace| object to manage all subsequent taxon references.
    If not given, a new |TaxonNamespace| object will be created.

    A |DataSet| can be placed in unattached taxon set mode by calling the instance method :meth:`~dendropy.dataobject.dataset.DataSet.unattach_taxon_namespace()`.
    This will restore the default behavior of a multiple taxon set |DataSet|.

    Note that placing a |DataSet| in attached taxon set mode using does not affect existing data: only data parsed while the |DataSet| object has an attached |TaxonNamespace| will have their taxon references mapped to the attached |TaxonNamespace|.
    You can use the instance method :meth:`~dendropy.dataobject.dataset.DataSet.unify_taxa()` to remap all taxon references of existing |TreeList| and |CharacterMatrix| objects to a (new) single |TaxonNamespace| object.
