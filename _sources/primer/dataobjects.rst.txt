*****************************************
Introduction to Phylogenetic Data Objects
*****************************************

Primary Phylogenetic Data Objects
==================================

Phylogenetic data in DendroPy is represented by one or more objects of the following classes:

    |Taxon|
        A representation of an operational taxonomic unit, with an attribute, :attr:`~dendropy.datamodel.taxonmodel.Taxon.label`, corresponding to the taxon label.

    |TaxonNamespace|
        A collection of |Taxon| objects representing a distinct definition of taxa (for example, as specified explicitly in a NEXUS "TAXA" block, or implicitly in the set of all taxon labels used across a Newick tree file).

    |Tree|
        A collection of |Node| and |Edge| objects representing a phylogenetic tree.
        Each |Tree| object maintains a reference to a |TaxonNamespace| object in its attribute, :attr:`~dendropy.datamodel.treemodel.Tree.taxon_namespace`, which specifies the set of taxa that are referenced by the tree and its nodes. Each |Node| object has a :attr:`~dendropy.datamodel.treemodel.Node.taxon` attribute (which points to a particular |Taxon| object if there is an operational taxonomic unit associated with this node, or is |None| if not), a :attr:`~dendropy.datamodel.treemodel.Node.parent_node` attribute (which will be |None| if the |Node| has no parent, e.g., a root node), a |Edge| attribute, as well as a list of references to child nodes, a copy of which can be obtained by calling :meth:`~dendropy.datamodel.treemodel.Node.child_nodes()`.
        In addition, advanced operations with tree data often make use of a |Bipartition| object associated with each |Edge| on a |Tree| (see ":doc:`/primer/bipartitions`" for more information).


    |TreeList|
        A :class:`list` of |Tree| objects. A |TreeList| object has an attribute, :attr:`~dendropy.datamodel.treemodel.TreeList.taxon_namespace`, which specifies the set of taxa that are referenced by all member |Tree| elements. This is enforced when a |Tree| object is added to a |TreeList|, with the |TaxonNamespace| of the |Tree| object and all |Taxon| references of the |Node| objects in the |Tree| mapped to the |TaxonNamespace| of the |TreeList|.

    |CharacterMatrix|
        Representation of character data, with specializations for different data types: |DnaCharacterMatrix|, |RnaCharacterMatrix|, |ProteinCharacterMatrix|, |StandardCharacterMatrix|, |ContinuousCharacterMatrix|, etc. A |CharacterMatrix| can treated very much like a :class:`dict` object, with
        |Taxon| objects as keys and character data as values associated with those keys.

    |DataSet|
        A meta-collection of phylogenetic data, consisting of lists of multiple |TaxonNamespace| objects (:attr:`DataSet.taxon_namespaces`), |TreeList| objects (:attr:`DataSet.tree_lists`), and |CharacterMatrix| objects (:attr:`DataSet.char_matrices`).


    |TreeArray|
        A high-performance container designed to efficiently store and manage (potentially) large collections of structures of (potentially) large trees for processing.


Creating New (Empty) Objects
============================

All of the above names are imported into the the the :mod:`dendropy` namespace, and so to instantiate new, empty objects of these classes, you would need to import :mod:`dendropy`::

    >>> import dendropy
    >>> tree1 = dendropy.Tree()
    >>> tree_list11 = dendropy.TreeList()
    >>> dna1 = dendropy.DnaCharacterMatrix()
    >>> dataset1 = dendropy.DataSet()

Or import the names directly::

    >>> from dendropy import Tree, TreeList, DnaCharacterMatrix, DataSet
    >>> tree1 = Tree()
    >>> tree_list1 = TreeList()
    >>> dna1 = DnaCharacterMatrix()
    >>> dataset1 = DataSet()

More details on how to create and populate new objects of various kinds programmatically are given in later chapters (e.g., ":doc:`trees`", ":doc:`chars`", ":doc:`datasets`").

Reading, Writing, and Annotating Phylogenetic Data
==================================================

DendroPy provides a rich set of tools for reading and writing phylogenetic data
in various formats, such as NEXUS, Newick, PHYLIP, etc., with *many* options to
customize and control how the data is ingested and parsed, as well as formatted
and written-out.
For example::

    >>> import dendropy
    >>> tree_list1 = dendropy.TreeList()
    >>> tree_list1.read_from_path("pythonidae.mcmc1.nex",
    ...     schema="nexus",
    ...     collection_offset=0,
    ...     tree_offset=100)
    >>> tree_list1.read_from_path("pythonidae.mcmc2.nex",
    ...     schema="nexus",
    ...     collection_offset=0,
    ...     tree_offset=100)
    >>> tree_list1.write_to_path("combined.newick",
    ...     suppress_edge_lengths=True,
    ...     schema="newick")

These are covered in detail in the next chapter, ":doc:`/primer/reading_and_writing`".


Support is also available for adding, accessing, and managing rich and
expressive metadata annotations to many of the above objects and components of
those objects. This is covered in detail in the
":doc:`/primer/working_with_metadata_annotations`" chapter.


