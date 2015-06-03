************************
Partitions of Taxon Sets
************************

A number of different applications require a specification of a partition of a set of taxa.
For example, when calculating population genetic summary statistics for a multi-population sample or numbers of deep coalescences given a particular monophyletic groupings of taxa.
The :class:`~dendropy.datamodel.taxonmodel.TaxonNamespacePartition` object describes a partitioning of a :class:`~dendropy.datamodel.taxonmodel.TaxonNamespace` into an exhaustive set of mutually-exclusive :class:`~dendropy.datamodel.taxonmodel.TaxonNamespace` subsets.

There are four different ways to specify a partitioning scheme: by using a function, attribute, dictionary or list.
The first three of these rely on providing a mapping of a :class:`~dendropy.datamodel.taxonmodel.Taxon` object to a subset membership identifier, i.e., a string, integer or some other type of value that identifies the grouping. The last explicitly describes the grouping as a list of lists.

For example, consider the following::

    >>> seqstr = """\
    ... #NEXUS
    ...
    ... BEGIN TAXA;
    ...     DIMENSIONS NTAX=13;
    ...     TAXLABELS a1 a2 a3 b1 b2 b3 c1 c2 c3 c4 c5 d1 d2;
    ... END;
    ... BEGIN CHARACTERS;
    ...     DIMENSIONS NCHAR=7;
    ...     FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=.;
    ...     MATRIX
    ...         a1 ACCTTTG
    ...         a2 ACCTTTG
    ...         a3 ACCTTTG
    ...         b1 ATCTTTG
    ...         b2 ATCTTTG
    ...         b3 ACCTTTG
    ...         c1 ACCCTTG
    ...         c2 ACCCTTG
    ...         c3 ACCCTTG
    ...         c4 ACCCTTG
    ...         c5 ACCCTTG
    ...         d1 ACAAAAG
    ...         d2 ACCAAAG
    ...     ;
    ... END
    ... """
    >>> seqs = DnaCharacterMatrix.get_from_string(seqstr, 'nexus')
    >>> taxon_namespace = seqs.taxon_namespace

Here we have sequences sampled from four populations, with the population identified by the first character of the taxon label. To create a parition of the |TaxonNamespace| resulting from parsing the file, we call the :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.partition` method. This method takes one of the following four keyword arguments:

        ``membership_func``
            A function that takes a ``Taxon`` object as an argument and
            returns a a population membership identifier or flag
            (e.g., a string, an integer) .

        ``membership_attr_name``
            Name of an attribute of ``Taxon`` objects that serves as an
            identifier for subset membership.

        ``membership_dict``
            A dictionary with ``Taxon`` objects as keys and population
            membership identifier or flag as values (e.g., a string,
            an integer).

        ``membership_lists``
            A container of containers of ``Taxon`` objects, with every
            ``Taxon`` object in ``taxon_namespace`` represented once and only
            once in the sub-containers.

For example, using the membership function approach, we define a function that returns the first character of the taxon label, and pass it to the :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.partition` using the ``membership_func`` keyword argument::

    >>> def mf(t):
    ...     return t.label[0]
    ...
    >>> tax_parts = taxon_namespace.partition(membership_func=mf)

Or, as would probably be done with such a simple membership function in practice::

    >>> tax_parts = taxon_namespace.partition(membership_func=lambda x: x.label[0])

Either way, we would get the following results::

    >>> for s in tax_parts.subsets():
    ...     print(s.description())
    ...
    TaxonNamespace object at 0x101116838 (TaxonNamespace4312885304: 'a'): 3 Taxa
    TaxonNamespace object at 0x101116788 (TaxonNamespace4312885128: 'c'): 5 Taxa
    TaxonNamespace object at 0x101116730 (TaxonNamespace4312885040: 'd'): 2 Taxa
    TaxonNamespace object at 0x1011167e0 (TaxonNamespace4312885216: 'b'): 3 Taxa

We could also add a population identification attribute to each |Taxon| object, and use the ``membership_attr_name`` keyword to specify that subsets should be created based on the value of this attribute::

    >>> for t in taxon_namespace:
    ...     t.population = t.label[0]
    ...
    >>> tax_parts = taxon_namespace.partition(membership_attr_name='population')

The results are identical to that above::

    >>> for s in tax_parts.subsets():
    ...     print(s.description())
    ...
    TaxonNamespace object at 0x1011166d8 (TaxonNamespace4312884952: 'a'): 3 Taxa
    TaxonNamespace object at 0x1011165d0 (TaxonNamespace4312884688: 'c'): 5 Taxa
    TaxonNamespace object at 0x1011169f0 (TaxonNamespace4312885744: 'd'): 2 Taxa
    TaxonNamespace object at 0x101116680 (TaxonNamespace4312884864: 'b'): 3 Taxa

The third approach involves constructing a dictionary that maps |Taxon| objects to their identification label and passing this using the ``membership_dict`` keyword argument::

    >>> tax_pop_label_map = {}
    >>> for t in taxon_namespace:
    ...     tax_pop_label_map[t] = t.label[0]
    ...
    >>> tax_parts = taxon_namespace.partition(membership_dict=tax_pop_label_map)

Again, the results are the same::

    >>> for s in tax_parts.subsets():
    ...     print(s.description())
    ...
    TaxonNamespace object at 0x1011166e8 (TaxonNamespace4312884952: 'a'): 3 Taxa
    TaxonNamespace object at 0x1011165f0 (TaxonNamespace4312884688: 'c'): 5 Taxa
    TaxonNamespace object at 0x1011169f1 (TaxonNamespace4312885744: 'd'): 2 Taxa
    TaxonNamespace object at 0x101116620 (TaxonNamespace4312884864: 'b'): 3 Taxa

Finally, a list of lists can be constructed and passed using the ``membership_lists`` argument::

    >>> pops = []
    >>> pops.append(taxon_namespace[0:3])
    >>> pops.append(taxon_namespace[3:6])
    >>> pops.append(taxon_namespace[6:11])
    >>> pops.append(taxon_namespace[11:13])
    >>> tax_parts = taxon_namespace.partition(membership_lists=pops)

Again, a :class:`~dendropy.datamodel.taxonmodel.TaxonNamespacePartition` object with four |TaxonNamespace| subsets is the result, only this time the subset labels are based on the list indices::

    >>> subsets = tax_parts.subsets()
    >>> print(subsets)
    set([<TaxonNamespace object at 0x10069f838>, <TaxonNamespace object at 0x10069fba8>, <TaxonNamespace object at 0x101116520>, <TaxonNamespace object at 0x1011164c8>])
    >>> for s in subsets:
    ...     print(s.description())
    ...
    TaxonNamespace object at 0x10069f838 (TaxonNamespace4301912120: '0'): 3 Taxa
    TaxonNamespace object at 0x10069fba8 (TaxonNamespace4301913000: '1'): 3 Taxa
    TaxonNamespace object at 0x101116520 (TaxonNamespace4312884512: '3'): 2 Taxa
    TaxonNamespace object at 0x1011164c8 (TaxonNamespace4312884424: '2'): 5 Taxa

