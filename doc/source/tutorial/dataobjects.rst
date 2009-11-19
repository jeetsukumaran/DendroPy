#############################
Phylogenetic Data in DendroPy
#############################

Introduction
============

Phylogenetic data in DendroPy is represented by one or more objects of the following classes:

    :class:`~dendropy.dataobject.taxon.Taxon`
        A representation of an operational taxonomic unit, with an attribute, :attr:`~dendropy.dataobject.taxon.Taxon.label`, corresponding to the taxon label.

    :class:`~dendropy.dataobject.taxon.TaxonSet`
        A collection of :class:`~dendropy.dataobject.taxon.Taxon` objects representing a distinct definition of taxa (for example, as specified explicitly in a NEXUS "TAXA" block, or implicitly in the set of all taxon labels used across a NEWICK tree file).

    :class:`~dendropy.dataobject.tree.Tree`
        A collection of :class:`~dendropy.dataobject.tree.Node` and :class:`~dendropy.dataobject.tree.Edge` objects representing a phylogenetic tree.
        Each :class:`~dendropy.dataobject.tree.Tree` object maintains a reference to a :class:`~dendropy.dataobject.taxon.TaxonSet` object in its attribute, :attr:`~dendropy.dataobject.tree.Tree.taxon_set`, which specifies the set of taxa that are referenced by the tree and its nodes. Each :class:`~dendropy.dataobject.tree.Node` object has a :attr:`~dendropy.dataobject.tree.Node.taxon` attribute, which points to a particular :class:`~dendropy.dataobject.taxon.Taxon` object if there is an operational taxonomic unit associated with this node, or is :keyword:`None` if not.

    :class:`~dendropy.dataobject.tree.TreeList`
        A :class:`list` of :class:`~dendropy.dataobject.tree.Tree` objects. A :class:`~dendropy.dataobject.tree.TreeList` object has an attribute, :attr:`~dendropy.dataobject.tree.TreeList.taxon_set`, which specifies the set of taxa that are referenced by all member :class:`~dendropy.dataobject.tree.Tree` elements. This is enforced when a :class:`~dendropy.dataobject.tree.Tree` object is added to a :class:`~dendropy.dataobject.tree.TreeList`, with the :class:`~dendropy.dataobject.taxon.TaxonSet` of the :class:`~dendropy.dataobject.tree.Tree` object and all :class:`~dendropy.dataobject.taxon.Taxon` references of the :class:`~dendropy.dataobject.tree.Node` objects in the :class:`~dendropy.dataobject.tree.Tree` mapped to the :class:`~dendropy.dataobject.taxon.TaxonSet` of the :class:`~dendropy.dataobject.tree.TreeList`.

    :class:`~dendropy.dataobject.char.CharacterArray`
        Representation of character data, with specializations for different data types: :class:`~dendropy.dataobject.char.DnaCharacterArray`, :class:`~dendropy.dataobject.char.RnaCharacterArray`, :class:`~dendropy.dataobject.char.ProteinCharacterArray`, :class:`~dendropy.dataobject.char.StandardCharacterArray`, :class:`~dendropy.dataobject.char.ContinuousCharacterArray`, etc. A :class:`~dendropy.dataobject.char.CharacterArray` can treated very much like a :class:`dict` object, with
        :class:`~dendropy.dataobject.taxon.Taxon` objects as keys and character data as values associated with those keys.

    :class:`~dendropy.dataobject.dataset.DataSet`
        A meta-collection of phylogenetic data, consisting of lists of multiple :class:`~dendropy.dataobject.taxon.TaxonSet` objects (:attr:`~dendropy.dataobject.DataSet.taxon_sets`), :class:`~dendropy.dataobject.tree.TreeList` objects (:attr:`~dendropy.dataobject.DataSet.tree_lists`), and :class:`~dendropy.dataobject.CharacterArray` objects (:attr:`~dendropy.dataobject.DataSet.char_arrays`).

Creating New (Empty) Objects
============================

All of the above names are imported into the the the :mod:`dendropy` namespace, and so to instantiate new, empty objects of these classes, you would need to import :mod:`dendropy`::

    >>> import dendropy
    >>> tree1 = dendropy.Tree()
    >>> tree_list11 = dendropy.TreeList()
    >>> dataset1 = dendropy.DataSet()

Or import the names directly::

    >>> from dendropy import Tree, TreeList, DataSet
    >>> tree1 = Tree()
    >>> tree_list1 = TreeList()
    >>> dataset1 = DataSet()

Creating and Populating New Objects
===================================

The :class:`~dendropy.dataobject.tree.Tree`, :class:`~dendropy.dataobject.tree.TreeList`, and :class:`~dendropy.dataobject.dataset.DataSet` classes all support ":meth:`get_from_*()`" factory methods that allow for the simultaneous instantiation and population of the objects from a data source:

    - :meth:`get_from_stream`
        Takes a file or file-like object opened for reading the data source as the first argument, and a string specifying the format as the second.

    - :meth:`get_from_path`
        Takes a string specifying the path to the the data source file as the first argument, and a string specifying the format as the second.

    - :meth:`get_from_string`
        Takes a string specifying containing the source data as the first argument, and a string specifying the format as the second.

All these methods minimally take a source and format reference as arguments and return a new object of the given type populated from the given sourcee::

    >>> import dendropy
    >>> tree1 = dendropy.Tree.get_from_string("((A,B),(C,D))", format="newick")
    >>> tree_list1 = dendropy.TreeList.get_from_path("pythonidae.mcmc.nex", format="nexus")
    >>> dataset1 = dendropy.TreeList.get_from_stream(open("pythonidae.fasta", format="fasta")

The format specification can be one of: "nexus", "newick", "nexml", "dnafasta", "rnafasta", "proteinfasta" etc. Not all formats are supported for reading, and not all formats make sense for particular objects (for example, it would not make sense to try and instantiate a :class:`~dendropy.dataobject.tree.Tree` or :class:`~dendropy.dataobject.tree.TreeList` object from a FASTA-formatted data source).

All these methods take a `taxon_set` keyword that specifies the :class:`~dendropy.dataobject.taxon.TaxonSet` to use to manage the operational taxonomic units defined or referenced in the data source. If not given, a new :class:`~dendropy.dataobject.taxon.TaxonSet` will be created and used.

Depending on the particular type being created and data source, these factory methods may take additional keyword arguments.
For example, if a data source contains multiple trees, you may specify a particular tree to be parsed by passing the 0-based index of the tree to the ":meth:`get_from_*()`" of :class:`~dendropy.dataobject.tree.Tree`::

    >>> tree2 = dendropy.TreeList.get_from_path("pythonidae.mcmc.nex", format="nexus", from_index=199)

The object `tree2` is now a DendroPy representation of the 200th tree found in the specified :download:`file </examples/pythonidae.mcmc.nex>`. :class:`~dendropy.dataobject.tree.TreeList`, also takes a `from_index` keyword argument, and specifying this will result in all trees starting from the given index being parsed and added to the collection.

With a :class:`~dendropy.dataobject.dataset.DataSet`, you can request that only trees or only characters are parsed::

    >>> ds1 = dendropy.DataSet.get_from_path('data1.nex', 'nexus', exclude_chars=True)
    >>> ds2 = dendropy.DataSet.get_from_path('data1.nex', 'nexus', exclude_trees=True)

In addition to the factory methods, you can specify a data source to the constructor of the objects directly using the `stream` and `format` keywords::

    >>> tree2 = dendropy.Tree(stream="primates.tre", format="newick", from_index=2)
    >>> tree_list2 = dendropy.TreeList(stream="primates.tre", format="newick")
    >>> dataset2 = dendropy.DataSet(stream="primates.nex", format="nexus")

You can also clone existing objects (i.e., create a deep-copy of everything but the taxon references)::

    >>> tree3 = dendropy.Tree(tree2)
    >>> tree_list3 = dendropy.TreeList(tree_list2)
    >>> dataset3 = dendropy.DataSet(dataset2)

Reading and Populating (or Repopulating) Existing Objects
=========================================================

The :class:`~dendropy.dataobject.tree.Tree`, :class:`~dendropy.dataobject.tree.TreeList`, and :class:`~dendropy.dataobject.dataset.DataSet` classes all support a suite of ":meth:`read_from_*()`" instance methods that parallels the ":meth:`get_from_*()`" factory methods described above:

    - :meth:`read_from_stream`
        Takes a file or file-like object opened for reading the data source as the first argument, and a string specifying the format as the second.

    - :meth:`read_from_path`
        Takes a string specifying the path to the the data source file as the first argument, and a string specifying the format as the second.

    - :meth:`read_from_string`
        Takes a string specifying containing the source data as the first argument, and a string specifying the format as the second.

When called on an existing :class:`~dendropy.dataobject.tree.TreeList`, or :class:`~dendropy.dataobject.dataset.DataSet` object, these methods *add* the data from the data source to the object, whereas when called on an existing :class:`~dendropy.dataobject.tree.Tree` object, they *replace* the :class:`~dendropy.dataobject.tree.Tree` with the tree specified in the data source.
As with the ":meth:`get_from_*()`" methods, the format specification can be any supported and type-apppropriate format, such as "nexus", "newick", "nexml", "dnafasta", "rnafasta", "proteinfasta" etc.

For example, the following accumulates post-burn-in trees from a several different files into a single :class:`~dendropy.dataobject.tree.TreeList` object::

    >>> import dendropy
    >>> post_trees = dendropy.TreeList()
    >>> post_trees.read_from_path("pythonidae.nex.run1.t", "nexus", from_index=200)
    >>> post_trees.read_from_path("pythonidae.nex.run2.t", "nexus", from_index=200)
    >>> post_trees.read_from_path("pythonidae.nex.run3.t", "nexus", from_index=200)
    >>> post_trees.read_from_path("pythonidae.nex.run4.t", "nexus", from_index=200)
    >>> print(post_trees.description())
    TreeList object at 0x5508a0 (TreeList5572768):  3204 Trees

The :class:`~dendropy.dataobject.tree.TreeList` objects automatically handles taxon management, and esnures that all appended :class:`~dendropy.dataobject.tree.Tree` objects share the same :class:`~dendropy.dataobject.taxon.TaxonSet` reference. Thus all the :class:`~dendropy.dataobject.tree.Tree` o

    >>> t = dendropy.Tree()
    >>> t.read_from_path('pythonidae.mcmc.nex', 'nexus', from_index=200)
    <Tree object at 0x6e71470>
    >>> print(t.description())
    Tree object at 0x6e71470 (Tree150000016: 'rep.200000'): ('Morelia boeleni':0.106115,((('Bothrochilus boa':0.092919,'Python timoriensis':0.180712):0.020687,(((('Antaresia perthensis':0.167512,'Antaresia stimsoni':0.059787):0.033053,'Antaresia maculosa':0.146173):0.016954,(('Morelia carinata':0.100305,'Morelia bredli':0.114501):0.015794,'Morelia viridis':0.130131):0.004453):0.033047,'Liasis fuscus':0.166956):0.026128):0.004973,('Morelia oenpelliensis':0.084937,'Python brongersmai':0.245248):0.017803):0.030474,'Aspidites ramsayi':0.121686)
    >>> t.read_from_path('pythonidae.mcmc.nex', 'nexus', from_index=201)
    <Tree object at 0x6e71470>
    >>> print(t.description())
    Tree object at 0x6e71470 (Tree148637488: 'rep.201000'): (((('Morelia oenpelliensis':0.097643,'Python brongersmai':0.296685):0.041638,(((('Morelia carinata':0.124252,'Morelia bredli':0.11777):0.013171,'Morelia viridis':0.134126):0.016724,('Antaresia maculosa':0.141482,('Antaresia perthensis':0.150992,'Antaresia stimsoni':0.072758):0.026511):0.027869):0.025264,'Liasis fuscus':0.215957):0.007182):0.001742,('Python timoriensis':0.257266,'Bothrochilus boa':0.165363):0.028538):0.027447,'Morelia boeleni':0.163616,'Aspidites ramsayi':0.071978)



.. toctree::
    :hidden:

    createtaxa.rst
    createtrees.rst
    createtreelists.rst
