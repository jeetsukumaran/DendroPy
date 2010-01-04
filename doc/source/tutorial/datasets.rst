**********************
Working with Data Sets
**********************

The |DataSet| class provides for objects that allow you to manage multiple types of phylogenetic data.

It has three primary attributes:

    - :attr:`~dendropy.dataobject.dataset.DataSet.taxon_sets`, a list of all |TaxonSet|         objects in the |DataSet|, in the order that they were added or read.
    - :attr:`~dendropy.dataobject.dataset.DataSet.tree_lists`, a list of all |TreeList| objects in the |DataSet|, in the order that they were added or read.
    - :attr:`~dendropy.dataobject.dataset.DataSet.char_matrices`, a list of all |CharacterMatrix| objects in the |DataSet|, in the order that they were added or read.

|DataSet| Creation and Reading
===============================

Creating a new |DataSet| from a Data Source
--------------------------------------------

You can use the :meth:`get_from_stream()`, :meth:`get_from_path()`, and :meth:`get_from_string()` factory class methods for simultaneously instantiating and populating an object, taking a data source as the first argument and a data format or schema specification as the second:

    >>> import dendropy
    >>> ds = dendropy.DataSet.get_from_path('pythonidae.nex', 'nexus')

Valid schema specification strings include: "``nexus``", "``newick``", "``nexml``", "``dnafasta``", "``rnafasta``", "``proteinfasta``" etc.

Reading into an Existing |DataSet| from a Data Source
-----------------------------------------------------

The :meth:`read_from_stream()`, :meth:`read_from_path()`, and :meth:`read_from_string()` instance methods for populating existing objects are also supported, taking the same arguments:

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.attach_taxon_set()
    >>> ds = dendropy.DataSet.read_from_path('pythonidae.cytb.fasta', 'dnafasta')
    >>> ds = dendropy.DataSet.read_from_path('pythonidae.mle.nex', 'nexus')

Note how the :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_set()` method is called before invoking any :meth:`read_from_*()` statements, to ensure that all the taxon references in the data sources get mapped to the same |TaxonSet| instance.

Cloning an Existing |DataSet|
-----------------------------

You can also clone an existing |DataSet| object by passing it as an argument to the |DataSet| constructor:

    >>> import dendropy
    >>> ds1 = dendropy.DataSet.get_from_path('pythonidae.cytb.fasta', 'dnafasta')
    >>> ds2 = dendropy.DataSet(ds1)

Following this, ``ds2`` will be a *full* deep-copy clone of ``ds1``, with distinct and independent, but identical, |Taxon|, |TaxonSet|, |TreeList|, |Tree| and |CharacterMatrix| objects.
Note that, in distinction to the similar cloning methods of |Tree| and |TreeList|, even the |Taxon| and |TaxonSet| objects are cloned, meaning that you manipulate the |Taxon| and |TaxonSet| objects of ``ds2`` without in any way effecting those of ``ds1``.

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
This method will remove all existing |TaxonSet| objects from the |DataSet|, create and add a new one, and then map all taxon references in all contained |TreeList| and |CharacterMatrix| objects to this new, unified |TaxonSet|.

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

Here, again, we call the :meth:`~dendropy.dataobject.dataset.DataSet.unify_taxa()` to map all taxon references to the same, common, unified |TaxonSet|.

.. _Customizing_Data_Set_Creation_and_Reading:

Customizing Data Set Creation and Reading
------------------------------------------

You can control how data is parsed from a data source using the following keywords passed to any :meth:`get_from_*()` or :meth:`read_from_*()` method of a |DataSet| object:

    ``attached_taxon_set``
        If :keyword:`True`, then a new |TaxonSet| object will be created and added to the :attr:`~dendropy.dataobject.dataset.DataSet.taxon_sets` list of the |DataSet| object, and the |DataSet| object will be placed in "attached" (or single) taxon set mode, i.e., all taxa in any data sources parsed or read will be mapped to the same |TaxonSet| object. By default, this is :keyword:`False`, resulting in a multi-taxon set mode |DataSet| object.

    ``taxon_set``
        A |TaxonSet| object that will be used to manage **all** taxon references in the data source.
        Every time a data source is parsed, by default at least one new |TaxonSet| object will be created to manage the taxa defined in the data source.
        If the data source defines multiple collections of taxa (as is possible with, for example, the NEXML schema, or the Mesquite variant of the NEXUS schema), then multiple new |TaxonSet| object will be created.
        By passing a |TaxonSet| object through the ``taxon_set`` keyword, you can force DendroPy to use the same |TaxonSet| object for all taxon references.

    ``exclude_trees``
        A boolean value indicating whether or not tree data should be parsed from the data source.
        Default value is :keyword:`False`, i.e., all tree data will be included.

    ``exclude_chars``
        A boolean value indicating whether or not character data should be parsed from the data source.
        Default value is :keyword:`False`, i.e., all character data will be included.

|DataSet| Saving and Writing
=============================

The :meth:`write_to_stream()`, and :meth:`write_to_path()` instance methods allow you to write the data of a |DataSet| object to a file-like object or a file path respectively.
These methods take a file-like object (in the case of :meth:`write_to_stream()`) or a string specifying a filepath (in the case of :meth:`write_to_path()`) as the first argument, and a format or schema specification string as the second argument.

The following example aggregates the post-burn in MCMC samples from a series of NEXUS-formatted tree files into a single |TreeList|, then, adds the |TreeList| as well as the original character data into a single |DataSet| object, which is then written out as NEXUS-formatted file:

.. literalinclude:: /examples/dsrw1.py
    :linenos:


If you do not want to actually write to a file, but instead simply need a string representing the data in a particular format, you can call the instance method :meth:`as_string()`, passing a schema or format specification string as the first argument:

    >>> import dendropy
    >>> ds = dendropy.DataSet(attached_taxon_set=True)
    >>> ds.read_from_path('pythonidae.cytb.fasta', 'dnafasta')
    >>> s = ds.as_string('nexus')
