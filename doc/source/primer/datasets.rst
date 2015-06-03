*********
Data Sets
*********

The |DataSet| class provides for objects that allow you to manage multiple types of phylogenetic data.

It has three primary attributes:

    :attr:`~dendropy.datamodel.datasetmodel.DataSet.taxon_namespaces`
        A list of all |TaxonNamespace| objects in the |DataSet|, in the order that they were added or read, include |TaxonNamespace| objects added implicitly through being associated with added |TreeList| or |CharacterMatrix| objects.

    :attr:`~dendropy.datamodel.datasetmodel.DataSet.tree_lists`
        A list of all |TreeList| objects in the |DataSet|, in the order that they were added or read.

    :attr:`~dendropy.datamodel.datasetmodel.DataSet.char_matrices`
        A list of all |CharacterMatrix| objects in the |DataSet|, in the order that they were added or read.

|DataSet| Creation and Reading
===============================

Reading and Writing |DataSet| Objects
-------------------------------------

You can use the :meth:`~dendropy.datamodel.datasetmodel.DataSet.get` factory class method for simultaneously instantiating and populating |DataSet| object, taking a data source as the first argument and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` ("``nexus``", "``newick``", "``nexml``", "``fasta``", "``phylip``", etc.) as the second::

    >>> import dendropy
    >>> ds = dendropy.DataSet.get(
        path='pythonidae.nex',
        schema='nexus')

The :meth:`~dendropy.datamodel.datasetmodel.DataSet.read` instance method for reading additional data into existing objects are also supported, taking the same arguments (i.e., a data source, a :ref:`schema specification string <Specifying_the_Data_Source_Format>`, as well as optional :keyword arguments to customize the parse behavior):

.. literalinclude:: /examples/ds1.py

.. Note::

    Note how the :meth:`~dendropy.datamodel.datasetmodel.DataSet.attach_taxon_namespace()` method is called before invoking any ":meth:`~dendropy.datamodel.datasetmodel.DataSet.read`" statements, to ensure that all the taxon references in the data sources get mapped to the same |TaxonNamespace| instance.
    It is **HIGHLY** recommended that you do this, i.e., manage all data with the same |DataSet| instance under the same taxonomic namespace, unless you have a special reason to include multiple independent taxon "domains" in the same data set.

The ":meth:`~dendropy.datamodel.datasetmodel.DataSet.write`" method allows you to write the data of a |DataSet| to a file-like object or a file path
The following example aggregates the post-burn in MCMC samples from a series of NEXUS-formatted tree files into a single |TreeList|, then, adds the |TreeList| as well as the original character data into a single |DataSet| object, which is then written out as NEXUS-formatted file:

.. literalinclude:: /examples/dsrw1.py

If you do not want to actually write to a file, but instead simply need a string representing the data in a particular format, you can call the instance method :meth:`as_string()`, passing a :ref:`schema specification string <Specifying_the_Data_Writing_Format>` as the first argument::

    import dendropy
    ds = dendropy.DataSet()
    ds.read_from_path('pythonidae.cytb.fasta', 'dnafasta')
    s = ds.as_string('nexus')

or::

    dna1 = dendropy.DataSet.get(file=open("pythonidae.nex"), schema="nexus")
    s = dna1.as_string(schema="fasta")
    print(s)


In addition, fine-grained control over the reading and writing of data is available through various keyword arguments.
More information on reading operations is available in the :doc:`/primer/reading_and_writing` section.

.. Cloning an Existing |DataSet|
.. -----------------------------
..
.. You can also clone an existing |DataSet| object by passing it as an argument to the |DataSet| constructor::
..
..     >>> import dendropy
..     >>> ds1 = dendropy.DataSet.get(
..     ... path='pythonidae.cytb.fasta',
..     ... schema='dnafasta')
..     >>> ds2 = dendropy.DataSet(ds1)
..
.. Following this, ``ds2`` will be a *full* deep-copy clone of ``ds1``, with distinct and independent, but identical, |Taxon|, |TaxonNamespace|, |TreeList|, |Tree| and |CharacterMatrix| objects.
.. Note that, in distinction to the similar cloning methods of |Tree| and |TreeList|, even the |Taxon| and |TaxonNamespace| objects are cloned, meaning that you manipulate the |Taxon| and |TaxonNamespace| objects of ``ds2`` without in any way effecting those of ``ds1``.

Creating a New |DataSet| from Existing |TreeList| and |CharacterMatrix| Objects
-------------------------------------------------------------------------------

You can add independentally created or parsed data objects to a |DataSet| by passing them as unnamed arguments to the constructor:

.. literalinclude:: /examples/ds4.py

Note how we call the instance method :meth:`~dendropy.datamodel.datasetmodel.DataSet.unify_taxon_namespaces()` after the creation of the |DataSet| object.
This method will remove all existing |TaxonNamespace| objects from the |DataSet|, create and add a new one, and then map all taxon references in all contained |TreeList| and |CharacterMatrix| objects to this new, unified |TaxonNamespace|.

Adding Data to an Exisiting |DataSet|
-------------------------------------

You can add independentally created or parsed data objects to a |DataSet| using the :meth:`~dendropy.datamodel.datasetmodel.DataSet.add()` method::

.. literalinclude:: /examples/ds4.py

Here, again, we call the :meth:`~dendropy.datamodel.datasetmodel.DataSet.unify_taxon_namespaces()` to map all taxon references to the same, common, unified |TaxonNamespace|.

Taxon Management with Data Sets
===============================

The |DataSet| object, representing a meta-collection of phylogenetic data, differs in one important way from all the other phylogenetic data objects discussed so far with respect to taxon management, in that it is not associated with any particular |TaxonNamespace| object.
Rather, it maintains a list (in the property :attr:`~dendropy.datamodel.datasetmodel.DataSet.taxon_namespaces`) of *all* the |TaxonNamespace| objects referenced by its contained |TreeList| objects (in the property :attr:`~dendropy.datamodel.datasetmodel.DataSet.tree_lists`) and |CharacterMatrix| objects (in the property :attr:`~dendropy.datamodel.datasetmodel.DataSet.char_matrices`).

With respect to taxon management, |DataSet| objects operate in one of two modes: "detached taxon set" mode and "attached taxon set" mode.

Detached (Multiple) Taxon Set Mode
----------------------------------

In the "detached taxon set" mode, which is the default, |DataSet| object tracks all |TaxonNamespace| references of their other data members in the property :attr:`~dendropy.datamodel.datasetmodel.DataSet.taxon_namespaces`, but no effort is made at taxon management as such.
Thus, every time a data source is read with a "detached taxon set" mode |DataSet| object, by default, a new  |TaxonNamespace| object will be created and associated with the |Tree|, |TreeList|, or |CharacterMatrix| objects created from each data source, resulting in multiple |TaxonNamespace| independent references.
As such, "detached taxon set" mode |DataSet| objects are suitable for handling data with multiple distinct sets of taxa.

For example::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read(path="primates.nex", schema="nexus")
    >>> ds.read(path="snakes.nex", schema="nexus")

The dataset, ``ds``, will now contain two distinct sets of |TaxonNamespace| objects, one for the taxa defined in "primates.nex", and the other for the taxa defined for "snakes.nex".
In this case, this behavior is correct, as the two files do indeed refer to different sets of taxa.

However, consider the following::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read(path="pythonidae_cytb.fasta", schema="fasta", data_type="dna")
    >>> ds.read(path="pythonidae_aa.nex", schema="nexus")
    >>> ds.read(path="pythonidae_morphological.nex", schema="nexus")
    >>> ds.read(path="pythonidae.mle.tre", schema="nexus")

Here, even though all the data files refer to the same set of taxa, the resulting  |DataSet| object will actually have 4 distinct  |TaxonNamespace| objects, one for each of the independent reads, and a taxon with a particular label in the first file (e.g., "Python regius" of "pythonidae_cytb.fasta") will map to a completely distinct |Taxon| object than a taxon with the same label in the second file (e.g., "Python regius" of "pythonidae_aa.nex").
This is incorrect behavior, and to achieve the correct behavior with a multiple taxon set mode |DataSet| object, we need to explicitly pass a |TaxonNamespace| object to each of the :meth:`~dendropy.datamodel.datasetmodel.DataSet.read_from_path()` statements::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read(path="pythonidae_cytb.fasta", schema="fasta", data_type="dna")
    >>> ds.read(schema="pythonidae_aa.nex", "nexus", taxon_namespace=ds.taxon_namespaces[0])
    >>> ds.read(schema="pythonidae_morphological.nex", "nexus", taxon_namespace=ds.taxon_namespaces[0])
    >>> ds.read(schema="pythonidae.mle.tre", "nexus", taxon_namespace=ds.taxon_namespaces[0])
    >>> ds.write_to_path("pythonidae_combined.nex", "nexus")

In the previous example, the first :meth:`~dendropy.datamodel.datasetmodel.DataSet.read()` statement results in a new |TaxonNamespace| object, which is added to the :attr:`~dendropy.datamodel.datasetmodel.DataSet.taxon_namespaces` property of the |DataSet| object ``ds``.
This |TaxonNamespace| object gets passed via the ``taxon_namespace`` keyword to subsequent :meth:`~dendropy.datamodel.datasetmodel.DataSet.read_from_path()` statements, and thus as each of the data sources are processed, the taxon references get mapped to |Taxon| objects in the same, single, |TaxonNamespace| object.

While this approach works to ensure correct taxon mapping across multiple data object reads and instantiation, in this context, it is probably more convenient to use the |DataSet| in "attached taxon set" mode.
In fact, it is highly recommended that |DataSet| instances *always* use the "attached taxon set" mode, as, conceptually there are very few cases where a collection of data should span multiple independent taxon namespaces.

Attached (Single) Taxon Set Mode
--------------------------------
In the "attached taxon set" mode, |DataSet| objects ensure that the taxon references of all data objects that are added to them are mapped to the same |TaxonNamespace| object (at least one for each independent read or creation operation).
The "attached taxon set" mode is activated by calling the :meth:`~dendropy.datamodel.datasetmodel.DataSet.attach_taxon_namespace` method on a |DataSet| and passing in the |TaxonNamespace| to use::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> taxa = dendropy.TaxonNamespace(label="global")
    >>> ds.attach_taxon_namespace(taxa)
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")

Switching Between Attached and Detached Taxon Set Modes
-------------------------------------------------------
As noted above, you can use the :meth:`~dendropy.datamodel.datasetmodel.DataSet.attached_taxon_namespace()` method to switch a |DataSet| object to attached taxon set mode.
To restore it to multiple taxon set mode, you would use the :meth:`~dendropy.datamodel.datasetmodel.DataSet.detach_taxon_namespace()` method::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> taxa = dendropy.TaxonNamespace(label="global")
    >>> ds.attach_taxon_namespace(taxa)
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")
    >>> ds.detach_taxon_namespace()
    >>> ds.read_from_path("primates.nex", "nexus")

Here, the same |TaxonNamespace| object is used to manage taxon references for data parsed from the first four files, while the data from the fifth and final file gets its own, distinct, |TaxonNamespace| object and associated |Taxon| object references.


