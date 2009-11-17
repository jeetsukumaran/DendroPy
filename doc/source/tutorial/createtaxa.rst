*************************
Taxa and Taxon Management
*************************

Operational taxonomic units in DendroPy are represented by :class:`~dendropy.dataobject.taxon.Taxon` objects, and distinct collections of operational taxonomic units are represented by :class:`~dendropy.dataobject.taxon.TaxonSet` objects.
Every time a definition of taxa is encountered in a data source, for example, a "TAXA" block in a NEXUS file, a new :class:`~dendropy.dataobject.taxon.TaxonSet` object is created and populated with :class:`~dendropy.dataobject.taxon.Taxon` objects corresponding to the taxa defined in the data source.
Some data formats do not have explicit definition of taxa, e.g. a NEWICK tree source.
These nonetheless can be considered to have an implicit definition of a collection of operational taxonomic units given by the aggregate of all operational taxonomic units referenced in the data (i.e., the set of all distinct labels on trees in a NEWICK file).

It is important to recognize that, by default, DendroPy will create new :class:`~dendropy.dataobject.taxon.TaxonSet` object every time a data source is parsed (and, if the data source has multiple taxon objects, there may be more than one :class:`~dendropy.dataobject.taxon.TaxonSet` created).

Consider the following example::

    >>> import dendropy
    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick')
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick')
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

We now have two distinct :class:`~dendropy.dataobject.tree.Tree` objects, each associated with a distinct :class:`~dendropy.dataobject.taxon.TaxonSet` objects, each with its own set of :class:`~dendropy.dataobject.taxon.Taxon` objects that, while having the same labels, are distinct from one another::

    >>> t1.leaf_nodes()[0].taxon == t2.leaf_nodes()[0].taxon
    False
    >>> t1.leaf_nodes()[0].taxon.label == t2.leaf_nodes()[0].taxon.label
    True

This means that even though the tree shape and structure is identical between the two trees, they exist in different universes as far as DendroPy is concerned::

.. Most phylogenetic data objects that reference taxa (such as :class:`~dendropy.dataobject.tree.Tree`, :class:`~dendropy.dataobject.tree.TreeList`, :class:`~dendropy.dataobject.char.CharacterArray`, etc.) have an attribute named `taxon_set` that points to a :class:`~dendropy.dataobject.taxon.TaxonSet` object.


