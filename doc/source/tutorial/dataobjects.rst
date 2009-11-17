#############################
Phylogenetic Data in DendroPy
#############################

Phylogenetic data in DendroPy is represented by one or more objects of the following classes:

    :class:`~dendropy.dataobject.Taxon`
        A representation of an operational taxonomic unit, with an attribute, :attr:`~dendropy.dataobject.Taxon.label`, corresponding to the taxon label.

    :class:`~dendropy.dataobject.TaxonSet`
        A collection of :class:`~dendropy.dataobject.Taxon` objects representing a distinct definition of taxa (for example, as specified explicitly in a NEXUS "TAXA" block, or implicitly in the set of all taxon labels used across a NEWICK tree file).

    :class:`~dendropy.dataobject.Tree`
        A collection of :class:`~dendropy.dataobject.Node` and :class:`~dendropy.dataobject.Edge` objects representing a phylogenetic tree.
        Each :class:`~dendropy.dataobject.Tree` object maintains a reference to a :class:`~dendropy.dataobject.TaxonSet` object in its attribute, :attr:`~dendropy.dataobject.Tree.taxon_set`, which specifies the set of taxa that are referenced by the tree and its nodes. Each :class:`~dendropy.dataobject.Node` object has a :attr:`~dendropy.dataobject.Node.taxon` attribute, which points to a particular :class:`~dendropy.dataobject.Taxon` object if there is an operational taxonomic unit associated with this node, or is :keyword:`None` if not.

    :class:`~dendropy.dataobject.TreeList`
        A :class:`list` of :class:`~dendropy.dataobject.Tree` objects. A :class:`~dendropy.dataobject.TreeList` object has an attribute, :attr:`~dendropy.dataobject.TreeList.taxon_set`, which specifies the set of taxa that are referenced by all member :class:`~dendropy.dataobject.Tree` elements. This is enforced when a :class:`~dendropy.dataobject.Tree` object is added to a :class:`~dendropy.dataobject.TreeList`, with the :class:`~dendropy.dataobject.TaxonSet` of the :class:`~dendropy.dataobject.Tree` object and all :class:`~dendropy.dataobject.Taxon` references of the :class:`~dendropy.dataobject.Node` objects in the :class:`~dendropy.dataobject.Tree` mapped to the :class:`~dendropy.dataobject.TaxonSet` of the :class:`~dendropy.dataobject.TreeList`.

    :class:`~dendropy.dataobject.CharacterArray`
        Representation of character data, with specializations for different data types: :class:`~dendropy.dataobject.DnaCharacterArray`, :class:`~dendropy.dataobject.RnaCharacterArray`, :class:`~dendropy.dataobject.ProteinCharacterArray`, :class:`~dendropy.dataobject.StandardCharacterArray`, :class:`~dendropy.dataobject.ContinuousCharacterArray`, etc. A :class:`~dendropy.dataobject.CharacterArray` can treated very much like a :class:`dict` object, with
        :class:`~dendropy.dataobject.Taxon` objects as keys and character data as values associated with those keys.

    :class:`~dendropy.dataobject.DataSet`
        A meta-collection of phylogenetic data, consisting of lists of multiple :class:`~dendropy.dataobject.TaxonSet` objects (:attr:`~dendropy.dataobject.DataSet.taxon_sets`), :class:`~dendropy.dataobject.TreeList` objects (:attr:`~dendropy.dataobject.DataSet.tree_lists`), and :class:`~dendropy.dataobject.CharacterArray` objects (:attr:`~dendropy.dataobject.DataSet.char_arrays`).

.. toctree::
    :maxdepth: 3

    createtrees.rst
