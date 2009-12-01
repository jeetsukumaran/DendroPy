*************************
Taxa and Taxon Management
*************************

Operational taxonomic units in DendroPy are represented by :class:`~dendropy.dataobject.taxon.Taxon` objects, and distinct collections of operational taxonomic units are represented by :class:`~dendropy.dataobject.taxon.TaxonSet` objects.

Every time a definition of taxa is encountered in a data source, for example, a "TAXA" block in a NEXUS file, a new :class:`~dendropy.dataobject.taxon.TaxonSet` object is created and populated with :class:`~dendropy.dataobject.taxon.Taxon` objects corresponding to the taxa defined in the data source.
Some data formats do not have explicit definition of taxa, e.g. a NEWICK tree source.
These nonetheless can be considered to have an implicit definition of a collection of operational taxonomic units given by the aggregate of all operational taxonomic units referenced in the data (i.e., the set of all distinct labels on trees in a NEWICK file).

Every time a reference to a taxon is encountered in a data source, such as a taxon label in a tree or matrix statement in a NEXUS file, the current :class:`~dendropy.dataobject.taxon.TaxonSet` object is searched for corresponding :class:`~dendropy.dataobject.taxon.Taxon` object with a matching label (see below for details on how the match is made). If found, the :class:`~dendropy.dataobject.taxon.Taxon` object is used to represent the taxon. If not, a new :class:`~dendropy.dataobject.taxon.Taxon` object is created, added to the :class:`~dendropy.dataobject.taxon.TaxonSet` object, and used to represent the taxon.

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

While :class:`~dendropy.dataobject.tree.TreeList` objects ensure that all :class:`~dendropy.dataobject.tree.Tree` objects created, read or added using them all have the same :class:`~dendropy.dataobject.taxon.TaxonSet` object reference, if two :class:`~dendropy.dataobject.tree.TreeList` objects are independentally created, they will each have their own, distinct, :class:`~dendropy.dataobject.taxon.TaxonSet` object reference.
For example, if you want to read in two collections of trees and compare trees between the two collections, the following will **not** work:


    >>> import dendropy
    >>> mcmc1 = dendropy.TreeList.get_from_path('pythonidae.mcmc1.nex', 'nexus')
    >>> mcmc2 = dendropy.TreeList.get_from_path('pythonidae.mcmc2.nex', 'nexus')

Of course, reading both data sources into the same  :class:`~dendropy.dataobject.tree.TreeList` object *will* work insofar as ensuring all the :class:`~dendropy.dataobject.tree.Tree` objects have the same :class:`~dendropy.dataobject.taxon.TaxonSet`  reference, but then you will lose the distinction between the two sources, unless you keep track of the indexes of where one source begins and the other ends, which error-prone and tedious.
A better approach would be simply to create a :class:`~dendropy.dataobject.taxon.TaxonSet` object, and pass it to the factory methods of both  :class:`~dendropy.dataobject.tree.TreeList` objects::

    >>> import dendropy
    >>> taxa = dendropy.TaxonSet()
    >>> mcmc1 = dendropy.TreeList.get_from_path('pythonidae.mcmc1.nex', 'nexus', taxon_set=taxa)
    >>> mcmc2 = dendropy.TreeList.get_from_path('pythonidae.mcmc2.nex', 'nexus', taxon_set=taxa)

Now both ``mcmc1`` and ``mcmc2`` share the same :class:`~dendropy.dataobject.taxon.TaxonSet`, and thus so do the :class:`~dendropy.dataobject.tree.Tree` objects created within them, which means the :class:`~dendropy.dataobject.tree.Tree` objects can be compared both within and between the collections.

You can also pass the :class:`~dendropy.dataobject.taxon.TaxonSet` to the constructor of :class:`~dendropy.dataobject.tree.TreeList`.
So, for example, the following is logically identical to the previous::

    >>> import dendropy
    >>> taxa = dendropy.TaxonSet()
    >>> mcmc1 = dendropy.TreeList(taxon_set=taxa)
    >>> mcmc1.read_from_path('pythonidae.mcmc1.nex', 'nexus')
    >>> mcmc2 = dendropy.TreeList(taxon_set=taxa)
    >>> mcmc2.read_from_path('pythonidae.mcmc2.nex', 'nexus')

Taxon Management with Character Arrays
======================================

Taxon management with :class:`~dendropy.dataobject.char.CharacterArray`-derived objects work very much the same as it does with :class:`~dendropy.dataobject.tree.Tree` or :class:`~dendropy.dataobject.tree.TreeList objects`: every time a :class:`~dendropy.dataobject.char.CharacterArray`-derived object is independentally created or read, a new :class:`~dendropy.dataobject.taxon.TaxonSet` is created, unless an existing one is specified.
Thus, again, if you are creating multiple character arrays that refer to the same set of taxa, you will want to make sure to pass each of them a common :class:`~dendropy.dataobject.taxon.TaxonSet` reference::

    >>> import dendropy
    >>> taxa = dendropy.TaxonSet()
    >>> dna1 = dendropy.DnaCharacterArray.get_from_path("pythonidae_cytb.fasta", "dnafasta", taxon_set=taxa)
    >>> std1 = dendropy.ProteinCharacterArray.get_from_path("pythonidae_morph.nex", "nexus", taxon_set=taxa)

Taxon Management with Data Sets
===============================
The :class:`~dendropy.dataobject.dataset.DataSet` object, representing a meta-collection of phylogenetic data, differs in one important way from all the other phylogenetic data objects discussed so far with respect to taxon management, in that it is not associated with any particular :class:`~dendropy.dataobject.taxon.TaxonSet` object.
Rather, it maintains a list (in the property :attr:`~dendropy.dataobject.char.DataSet.taxon_sets`) of *all* the :class:`~dendropy.dataobject.taxon.TaxonSet` objects referenced by its contained :class:`~dendropy.dataobject.tree.TreeList` objects (in the property :attr:`~dendropy.dataobject.char.DataSet.tree_lists`) and :class:`~dendropy.dataobject.char.CharacterArray` objects (in the property :attr:`~dendropy.dataobject.char.DataSet.char_arrays`).

With respect to taxon management, :class:`~dendropy.dataobject.dataset.DataSet` objects operate in one of two modes: "multiple taxon set" mode and "fixed taxon set" mode.

Multiple Taxon Set Mode
-----------------------

In the "multiple taxon set" mode, which is the default, :class:`~dendropy.dataobject.dataset.DataSet` object tracks all :class:`~dendropy.dataobject.taxon.TaxonSet` references of their other data members in the property :attr:`~dendropy.dataobject.char.DataSet.taxon_sets`, but no effort is made at taxon management as such.
Thus, every time a data source is read with a "multiple taxon set" mode :class:`~dendropy.dataobject.dataset.DataSet` object, by deault, a new  :class:`~dendropy.dataobject.taxon.TaxonSet` object will be created and associated with the :class:`~dendropy.dataobject.tree.Tree`, :class:`~dendropy.dataobject.tree.TreeList`, or :class:`~dendropy.dataobject.char.CharacterArray` objects created from each data source, resulting in multiple :class:`~dendropy.dataobject.taxon.TaxonSet` independent references.
As such, "multiple taxon set" mode :class:`~dendropy.dataobject.dataset.DataSet` objects are suitable for handling data with multiple distinct sets of taxa.

For example::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read_from_path("primates.nex", "nexus")
    >>> ds.read_from_path("snakes.nex", "nexus")

The dataset, ``ds``, will now contain two distinct sets of :class:`~dendropy.dataobject.taxon.TaxonSet` objects, one for the taxa defined in "primates.nex", and the other for the taxa defined for "snakes.nex".
In this case, this behavior is correct, as the two files do indeed refer to different sets of taxa.

However, consider the following::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")

Here, even though all the data files refer to the same set of taxa, the resulting  :class:`~dendropy.dataobject.dataset.DataSet` object will actually have 4 distinct  :class:`~dendropy.dataobject.taxon.TaxonSet` objects, one for each of the independent reads, and a taxon with a particular label in the first file (e.g., "Python regius" of "pythonidae_cytb.fasta") will map to a completely distinct :class:`~dendropy.dataobject.taxon.Taxon` object than a taxon with the same label in the second file (e.g., "Python regius" of "pythonidae_aa.nex").
This is incorrect behavior, and to achieve the correct behavior with a multiple taxon set mode :class:`~dendropy.dataobject.dataset.DataSet` object, we need to explicitly pass a :class:`~dendropy.dataobject.taxon.TaxonSet` object to each of the :meth:`~dendropy.dataobject.dataset.DataSet.read_from_path()` statements::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus", taxon_set=ds.taxon_sets[0])
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus", taxon_set=ds.taxon_sets[0])
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus", taxon_set=ds.taxon_sets[0])
    >>> ds.write_to_path("pythonidae_combined.nex", "nexus")

In the previous example, the first :meth:`~dendropy.dataobject.dataset.DataSet.read_from_path()` statement results in a new :class:`~dendropy.dataobject.taxon.TaxonSet` object, which is added to the :attr:`~dendropy.dataobject.char.DataSet.taxon_sets` property of the :class:`~dendropy.dataobject.dataset.DataSet` object ``ds``.
This :class:`~dendropy.dataobject.taxon.TaxonSet` object gets passed via the ``taxon_set`` keyword to subsequent :meth:`~dendropy.dataobject.dataset.DataSet.read_from_path()` statements, and thus as each of the data sources are processed, the taxon references get mapped to :class:`~dendropy.dataobject.taxon.Taxon` objects in the same, single, :class:`~dendropy.dataobject.taxon.TaxonSet` object.

While this approach works to ensure correct taxon mapping across multiple data object reads and instantiation, in this context, it is probably more convenient to use the :class:`~dendropy.dataobject.dataset.DataSet` in "fixed taxon set" mode.

Fixed (Single) Taxon Set Mode
-----------------------------
In the "fixed taxon set" mode, :class:`~dendropy.dataobject.dataset.DataSet` objects ensure that the taxon references of all data objects that are added to them are mapped to the same :class:`~dendropy.dataobject.taxon.TaxonSet` object (at least one for each independent read or creation operation).
The "fixed taxon set" can be set by passing the keyword argument ``fixed_taxon_set=True`` to the constructor of the :class:`~dendropy.dataobject.dataset.DataSet` when instantiating a new :class:`~dendropy.dataobject.dataset.DataSet` object, or by calling :meth:`~dendropy.dataobject.dataset.DataSet.fix_taxon_set()` on an existing :class:`~dendropy.dataobject.dataset.DataSet` object

For example::

    >>> import dendropy
    >>> ds = dendropy.DataSet(fixed_taxon_set=True)
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")

Or alternatively::

    >>> import dendropy
    >>> ds.fix_taxon_set()
    >>> ds = dendropy.DataSet(fixed_taxon_set=True)
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")

Both of the above will result in only a single :class:`~dendropy.dataobject.taxon.TaxonSet` object that have all the taxa from the four data sources mapped to them.

A Word of Caution: Taxon Label Mapping
======================================
DendroPy maps taxon definitions encountered in a data source to :class:`~dendropy.dataobject.taxon.Taxon` objects by the taxon label.
The labels have to match **exactly** for the taxa to be correctly mapped, with the match being **case-sensitive**.
Thus, "Python regius", "PYTHON REGIUS", "python regious", "P. regious", etc. will all be considered as referring to distinct and different taxa.

Further quirks may arise due to some format-specific idiosyncracies.
For example, the NEXUS standard dictates that an underscore ("_") should be substituted for a space in all labels.
Thus, when reading a NEXUS or NEWICK source, the taxon labels "Python_regius" and "Python regius" are exactly equivalent, and will be mapped to the same :class:`~dendropy.dataobject.taxon.Taxon` object.

However, this underscore-to-space mapping does **not** take place when reading, for example, a FASTA format file.
Here, underscores are preserved, and thus "Python_regius" does not map to "Python regius".
This means that if you were to read a NEXUS file with the taxon label, "Python_regius", and later a read a FASTA file with the same taxon label, i.e., "Python_regius", these would map to different taxa!
This is illustrated by the following:

.. literalinclude:: /examples/taxon_labels1.py
    :linenos:

Which produces the following, almost certainly incorrect, result::

    TaxonSet object at 0x43b4e0 (TaxonSet4437216): 4 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python regious'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python sebae'
        [2] Taxon object at 0x22867d0 (Taxon36202448): 'Python_regious'
        [3] Taxon object at 0x2286830 (Taxon36202544): 'Python_sebae'

Even more confusingly, if this file is written out in NEXUS format, it would result in the space/underscore substitution taking place, resulting in two pairs of taxa with the same labels.

As such, if you plan on mixing sources from different formats, it is important to keep in mind the space/underscore substitution, and in data formats that do not have this convention, avoid underscores and use spaces instead:

.. literalinclude:: /examples/taxon_labels2.py
    :linenos:

Which results in the following, correct, behavior::

    TaxonSet object at 0x43b4e0 (TaxonSet4437216): 2 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python regious'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python sebae'

Alternatively, you can wrap the underscore-bearing labels in the NEXUS/NEWICK source in quotes, which preserves them from being substituted for spaces:

.. literalinclude:: /examples/taxon_labels3.py
    :linenos:

Which results in the following, also correct, behavior::

    TaxonSet object at 0x43b4e0 (TaxonSet4437216): 2 Taxa
        [0] Taxon object at 0x22867b0 (Taxon36202416): 'Python_regious'
        [1] Taxon object at 0x2286810 (Taxon36202512): 'Python_sebae'




