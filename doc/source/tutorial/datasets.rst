**********************
Working with Data Sets
**********************

The |DataSet| class provides for objects that allow you to manage multiple types of phylogenetic data.

It has three primary attributes:

    - :attr:`~dendropy.dataobject.dataset.DataSet.taxon_sets`, a list of all |TaxonSet|         objects in the |DataSet|, in the order that they were added or read.
    - :attr:`~dendropy.dataobject.dataset.DataSet.tree_lists`, a list of all |TreeList| objects in the |DataSet|, in the order that they were added or read.
    - :attr:`~dendropy.dataobject.dataset.DataSet.char_matrices`, a list of all |CharacterMatrix| objects in the |DataSet|, in the order that they were added or read.

Basic Instantiation and Reading
-------------------------------

The |DataSet| class provides the :meth:`get_from_*()` factory class methods for simultaneously instantiating and populating |DataSet| objects, taking a data source as the first argument and a data format or schema specification as the second:

    >>> import dendropy
    >>> ds = dendropy.DataSet.get_from_path('pythonidae.nex', 'nexus')

The :meth:`read_from_*()` instance methods for populating existing objects, are also supported, taking the same arguments:

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.attach_taxon_set()
    >>> ds = dendropy.DataSet.read_from_path('pythonidae.cytb.fasta', 'dnafasta')
    >>> ds = dendropy.DataSet.read_from_path('pythonidae.mle.nex', 'nexus')

Note how the :meth:`~dendropy.dataobject.dataset.DataSet.attach_taxon_set()` method is called before invoking any :meth:`read_from_*()` statements, to ensure that all the taxon references in the data sources get mapped to the same |TaxonSet| instance.

Adding Data to a |DataSet|
--------------------------

You can add independentally created or parsed data objects to a |DataSet| by passing them as unnamed arguments to the constructor:

    >>> import dendropy
    >>> treelist1 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run1.t', 'nexus')
    >>> treelist2 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run1.t', 'nexus')
    >>> treelist3 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run1.t', 'nexus')
    >>> treelist4 = dendropy.TreeList.get_from_path('pythonidae_cytb.mb.run1.t', 'nexus')
    >>> cytb = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.fasta', 'dnafasta')
    >>> ds = dendropy.DataSet(cytb, treelist1, treelist2, treelist3, treelist4)
    >>> ds.reindex_taxa()

Note how we call the instance method :meth:`~dendropy.dataobject.dataset.DataSet.reindex_taxa()` after the creation of the |DataSet| object.
This method will remove all existing |TaxonSet| objects from the |DataSet|, create and add a new one, and then map all taxon references in all contained |TreeList| and |CharacterMatrix| objects to this new, unified |TaxonSet|.

.. _Customizing_Data_Set_Creation_and_Reading:

Customizing Data Set Creation and Reading
=========================================

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
