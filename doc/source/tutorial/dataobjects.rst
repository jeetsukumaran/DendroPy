*****************************************
Introduction to Phylogenetic Data Objects
*****************************************

Types of Phylogenetic Data Objects
==================================

Phylogenetic data in DendroPy is represented by one or more objects of the following classes:

    |Taxon|
        A representation of an operational taxonomic unit, with an attribute, :attr:`~dendropy.dataobject.taxon.Taxon.label`, corresponding to the taxon label.

    |TaxonSet|
        A collection of |Taxon| objects representing a distinct definition of taxa (for example, as specified explicitly in a NEXUS "TAXA" block, or implicitly in the set of all taxon labels used across a Newick tree file).

    |Tree|
        A collection of |Node| and |Edge| objects representing a phylogenetic tree.
        Each |Tree| object maintains a reference to a |TaxonSet| object in its attribute, :attr:`~dendropy.dataobject.tree.Tree.taxon_set`, which specifies the set of taxa that are referenced by the tree and its nodes. Each |Node| object has a :attr:`~dendropy.dataobject.tree.Node.taxon` attribute (which points to a particular |Taxon| object if there is an operational taxonomic unit associated with this node, or is |None| if not), a :attr:`~dendropy.dataobject.tree.Node.parent_node` attribute (which will be |None| if the |Node| has no parent, e.g., a root node), a |Edge| attribute, as well as a list of references to child nodes, a copy of which can be obtained by calling :meth:`~dendropy.dataobject.tree.Node.child_nodes()`.

    |TreeList|
        A :class:`list` of |Tree| objects. A |TreeList| object has an attribute, :attr:`~dendropy.dataobject.tree.TreeList.taxon_set`, which specifies the set of taxa that are referenced by all member |Tree| elements. This is enforced when a |Tree| object is added to a |TreeList|, with the |TaxonSet| of the |Tree| object and all |Taxon| references of the |Node| objects in the |Tree| mapped to the |TaxonSet| of the |TreeList|.

    |CharacterMatrix|
        Representation of character data, with specializations for different data types: |DnaCharacterMatrix|, |RnaCharacterMatrix|, |ProteinCharacterMatrix|, |StandardCharacterMatrix|, |ContinuousCharacterMatrix|, etc. A |CharacterMatrix| can treated very much like a :class:`dict` object, with
        |Taxon| objects as keys and character data as values associated with those keys.

    |DataSet|
        A meta-collection of phylogenetic data, consisting of lists of multiple |TaxonSet| objects (:attr:`~dendropy.dataobject.DataSet.taxon_sets`), |TreeList| objects (:attr:`~dendropy.dataobject.DataSet.tree_lists`), and |CharacterMatrix| objects (:attr:`~dendropy.dataobject.DataSet.char_matrices`).

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

Reading and Writing Phylogenetic Data
=====================================

DendroPy provides a rich set of tools for reading and writing phylogenetic data in various formats, such as NEXUS, Newick, PHYLIP, etc. These are covered in detail in the following ":doc:`/tutorial/reading`" and ":doc:`/tutorial/writing`" chapters respectively.


