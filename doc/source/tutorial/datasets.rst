**********************
Working with Data Sets
**********************

The |DataSet| class provides for objects that allow you to manage multiple types of phylogenetic data.


.. _Customizing_Data_Set_Creation_and_Reading:

Customizing Data Set Creation and Reading
=========================================

A |DataSet| can be instantiated by:

    - passing a

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
