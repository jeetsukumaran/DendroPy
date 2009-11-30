*************************
Taxa and Taxon Management
*************************

Operational taxonomic units in DendroPy are represented by :class:`~dendropy.dataobject.taxon.Taxon` objects, and distinct collections of operational taxonomic units are represented by :class:`~dendropy.dataobject.taxon.TaxonSet` objects.
Every time a definition of taxa is encountered in a data source, for example, a "TAXA" block in a NEXUS file, a new :class:`~dendropy.dataobject.taxon.TaxonSet` object is created and populated with :class:`~dendropy.dataobject.taxon.Taxon` objects corresponding to the taxa defined in the data source.
Some data formats do not have explicit definition of taxa, e.g. a NEWICK tree source.
These nonetheless can be considered to have an implicit definition of a collection of operational taxonomic units given by the aggregate of all operational taxonomic units referenced in the data (i.e., the set of all distinct labels on trees in a NEWICK file).

Taxon Management with Trees
===========================

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

This means that even though the tree shape and structure is identical between the two trees, they exist in different universes as far as DendroPy is concerned, and many operations that involving comparing trees will fail::

    >>> from dendropy import treecalc
    >>> treecalc.robinson_foulds_distance(t1, t2)
    ------------------------------------------------------------
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>

      File "/Users/jeet/Documents/Projects/dendropy/dendropy/treecalc.py", line 263, in robinson_foulds_distance
        value_type=float)

      File "/Users/jeet/Documents/Projects/dendropy/dendropy/treecalc.py", line 200, in splits_distance
        % (hex(id(tree1.taxon_set)), hex(id(tree2.taxon_set))))

    TypeError: Trees have different TaxonSet objects: 0x101f630 vs. 0x103bf30

The solution is to explicitly specify the same ``taxon_set`` when creating the trees. In DendroPy all phylogenetic data classes that are associated with :class:`~dendropy.dataobject.taxon.TaxonSet` objects have constructors, factory methods, and ``read_from_*`` methods take a specific :class:`TaxonSet` object as an argument using the ``taxon_set`` a keyword. For example::

    >>> taxa = dendropy.TaxonSet()
    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick', taxon_set=taxa)
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick', taxon_set=taxa)
    >>> treecalc.robinson_foulds_distance(t1, t2)
    0.0

Taxon Management with Tree Lists
================================

The :class:`~dendropy.dataobject.tree.TreeList` class is designed to manage collections of :class:`~dendropy.dataobject.tree.Tree` objects that share the same :class:`~dendropy.dataobject.taxon.TaxonSet`.
As :class:`~dendropy.dataobject.tree.Tree` objects are appended to a :class:`~dendropy.dataobject.tree.TreeList` object, the :class:`~dendropy.dataobject.tree.TreeList` object will automatically take care of remapping the :class:`~dendropy.dataobject.taxon.TaxonSet` and associated :class:`~dendropy.dataobject.taxon.Taxon` objects::

    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick')
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick')
    >>> print(repr(t1.taxon_set))
    <TaxonSet object at 0x1243a20>
    >>> repr(t1.taxon_set)
    '<TaxonSet object at 0x1243a20>'
    >>> repr(t2.taxon_set)
    '<TaxonSet object at 0x12439f0>'
    >>> trees = dendropy.TreeList()
    >>> trees.append(t1)
    >>> trees.append(t2)
    >>> repr(t1.taxon_set)
    '<TaxonSet object at 0x1243870>'
    >>> repr(t2.taxon_set)
    '<TaxonSet object at 0x1243870>'
    >>> treecalc.robinson_foulds_distance(t1, t2)
    0.0

The same applies when using the :meth:`read_from_*` method of a :class:`~dendropy.dataobject.tree.TreeList` object: all trees read from the data source will be assigned the same :class:`~dendropy.dataobject.taxon.TaxonSet` object, and the taxa referenced in the tree definition will be mapped to corresponding :class:`~dendropy.dataobject.taxon.Taxon` objects, identified by label, in the :class:`~dendropy.dataobject.taxon.TaxonSet`, with new :class:`~dendropy.dataobject.taxon.Taxon` objects created if no suitable match is found.
