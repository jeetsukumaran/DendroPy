#############################
Phylogenetic Data in DendroPy
#############################

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

.. toctree::
    :hidden:

    createtaxa.rst
    createtrees.rst
    createtreelists.rst
