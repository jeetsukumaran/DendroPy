*****************
Trees in DendroPy
*****************

Trees in |DendroPy|_ are represented by the class ``Tree``. All trees (generally) belong to a particular ``TreesBlock``, which is derived from a Python ``list``. 

Each ``Tree`` object has an attribute, ``taxa_block``, which is a ``TaxaBlock`` object, and manages all the ``Taxon`` objects associated with the tree.
The ``TaxaBlock`` object referenced by a ``Tree`` object's ``taxa_block`` might be shared by many other elements of the dataset, including other ``Tree`` objects and ``CharactersBlock`` objects, so any modification of elements of a ``Tree`` object's ``taxa_block`` will probably have dataset-wide effects.
That is, if you were to change the label of a ``Taxon`` object maintained by a particular ``Tree`` object's ``taxa_block``, all other ``Tree`` objects in the dataset referencing the same ``TaxaBlock`` will be effected.

Every ``Tree`` object has a ``seed_node`` attribute. If the tree is rooted (``<tree>.is_rooted==True``), then this is the root node. If the tree is not rooted, however, then this is an artificial node that serves as the "starting point" for the tree. 

The ``seed_node``, like every other node on the tree, is a ``Node`` object. 
Every ``Node`` object maintains a list of its immediate child ``Node`` objects as well as a reference to its parent ``Node`` object. 
You can request a shallow-copy ``list`` of child ``Node`` objects using the ``Node.child_nodes()`` method, and you can access the parent ``Node`` object directly through the ``Node.parent_node`` attribute.
By definition, the ``seed_node`` has no parent node (``parent_node==None``), leaf nodes have no child nodes, and internal nodes have both parent nodes and child nodes.

Every ``Node`` object also has an ``edge`` attribute, which points to an ``Edge`` object representing the branch subtending the node. ``Edge`` objects have a ``length`` attribute, which is typically either a ``float`` or ``int`` value, representing the weight or length of the branch.
If branch lengths have not been specified, then the value of ``length`` is ``None``.
Even if the source tree has had branch lengths specified, if the tree is unrooted, then the edge of the ``seed_node`` is usually ``None``.

``Node`` objects also have a ``label`` and ``taxon`` attribute. These are, by default, set to ``None``, but leaf nodes almost always have their ``taxon`` attribute set, pointing to a ``Taxon`` object associated with that tip of the tree. The ``label`` attribute will be set if the source tree has internal node labels, though, of course, you can also assign a value to this programmatically.

.. toctree::
    :maxdepth: 2
    
    traversing_trees.rst