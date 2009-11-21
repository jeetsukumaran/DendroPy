*****************
Trees in DendroPy
*****************

Trees in are represented by the class :class:`dendropy.dataobject.tree.Tree`.

Each :class:`dendropy.dataobject.tree.Tree` object has an attribute, :attr:`dendropy.dataobject.tree.Tree.taxon_set`, which is a ``TaxaBlock`` object, and manages all the :class:`dendropy.dataobject.taxon.Taxon` objects associated with the tree.
The ``TaxaBlock`` object referenced by a :class:`dendropy.dataobject.tree.Tree` object's :attr:`dendropy.dataobject.tree.Tree.taxon_set` might be shared by many other elements of the dataset, including other :class:`dendropy.dataobject.tree.Tree` objects and :class:`dendropy.dataobject.char.CharacterArray` objects, so any modification of elements of a :class:`dendropy.dataobject.tree.Tree` object's :attr:`dendropy.dataobject.tree.Tree.taxon_set` will probably have dataset-wide effects.
That is, if you were to change the label of a :class:`dendropy.dataobject.taxon.Taxon` object maintained by a particular :class:`dendropy.dataobject.tree.Tree` object's :attr:`dendropy.dataobject.tree.Tree.taxon_set`, all other :class:`dendropy.dataobject.tree.Tree` objects in the dataset referencing the same ``TaxaBlock`` will be effected.

Every :class:`dendropy.dataobject.tree.Tree` object has a :attr:`dendropy.dataobject.tree.Tree.seed_node` attribute. If the tree is rooted, then this is the root node. If the tree is not rooted, however, then this is an artificial node that serves as the "starting point" for the tree.

The :attr:`dendropy.dataobject.tree.Tree.seed_node`, like every other node on the tree, is a :class:`dendropy.dataobject.tree.Node` object.
Every :class:`dendropy.dataobject.tree.Node` object maintains a list of its immediate child :class:`dendropy.dataobject.tree.Node` objects as well as a reference to its parent :class:`dendropy.dataobject.tree.Node` object.
You can request a shallow-copy :func:`list` of child :class:`dendropy.dataobject.tree.Node` objects using the :meth:`dendropy.dataobject.tree.child_nodes()` method, and you can access the parent :class:`dendropy.dataobject.tree.Node` object directly through the ``Node.parent_node`` attribute.
By definition, the :attr:`dendropy.dataobject.tree.Tree.seed_node` has no parent node, leaf nodes have no child nodes, and internal nodes have both parent nodes and child nodes.

Every :class:`dendropy.dataobject.tree.Node` object also has an :attr:`dendropy.dataobject.tree.Node.edge` attribute, which points to an :class:`dendropy.dataobject.tree.Edge` object representing the branch subtending the node. :class:`dendropy.dataobject.tree.Edge` objects have a :attr:`dendropy.dataobject.tree.Edge.length` attribute, which is typically either a ``float`` or ``int`` value, representing the weight or length of the branch.
If branch lengths have not been specified, then the value of :attr:`dendropy.dataobject.tree.Edge.length` is :keyword:`None`.
Even if the source tree has had branch lengths specified, if the tree is unrooted, then the edge of the :attr:`dendropy.dataobject.tree.Tree.seed_node` is usually :keyword:`None`.

:class:`dendropy.dataobject.tree.Node` objects also have a :attr:`dendropy.dataobject.tree.Node.label` and :attr:`dendropy.dataobject.tree.Node.taxon` attribute. Leaf nodes usually have their :attr:`dendropy.dataobject.tree.Node.taxon` attribute set, pointing to :class:`dendropy.dataobject.taxon.Taxon` object associated with that tip of the tree. The :attr:`dendropy.dataobject.tree.Node.label` attribute will be set if the source tree has internal node labels, though, of course, you can also assign a value to this programmatically.

.. toctree::
    :maxdepth: 2

    traversing_trees.rst
