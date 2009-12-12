**********************
Working with Data Sets
**********************

The :class:`~dendropy.dataobject.dataset.DataSet` class provides for objects that allow you to manage multiple types of phylogenetic data.

Customizing Data Set Creation and Reading
=========================================

You can control how data is parsed from a data source using the following keywords passed to any :meth:`get_from_*()` or :meth:`read_from_*()` method of a :class:`~dendropy.dataobject.dataset.DataSet` object:

    ``taxon_set``
        A :class:`~dendropy.dataobject.taxon.TaxonSet` object that will be used to manage **all** taxon references in the data source.
        Every time a data source is parsed, by default at least one new :class:`~dendropy.dataobject.taxon.TaxonSet` object will be created to manage the taxa defined in the data source.
        If the data source defines multiple collections of taxa (as is possible with, for example, the NEXML format, or the Mesquite variant of the NEXUS format), then multiple new :class:`~dendropy.dataobject.taxon.TaxonSet` object will be created.
        By passing a :class:`~dendropy.dataobject.taxon.TaxonSet` object through the ``taxon_set`` keyword, you can force DendroPy to use the same :class:`~dendropy.dataobject.taxon.TaxonSet` object for all taxon references.

    ``exclude_trees``
        A boolean value indicating whether or not tree data should be parsed from the data source.
        Default value is :keyword:`False`, i.e., all tree data will be included.

    ``exclude_chars``
        A boolean value indicating whether or not character data should be parsed from the data source.
        Default value is :keyword:`False`, i.e., all character data will be included.
