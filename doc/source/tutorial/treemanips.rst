***********************************
Tree Manipulation and Restructuring
***********************************

The |Tree| class provides both low-level and high-level methods for manipulating tree structure.

Low-level methods are associated with |Node| objects, and allow to restructure the relationships between nodes at a fine level: :meth:`~dendropy.dataobject.tree.Node.add_child`, :meth:`~dendropy.dataobject.tree.Node.new_child`, :meth:`~dendropy.dataobject.tree.Node.remove_child`, etc.

In most cases, however, you will be using high-level methods to restructure |Tree| objects.

In all cases, if any part of the |Tree| object's structural relations change, *and* you are interested in calculating any metrics or statistics on the tree or comparing the tree to another tree, you need to call :meth:`~dendropy.dataobject.tree.Tree.update_splits()` on the object to update the internal splits hash representation.

Rooting, Derooting and Rerooting
================================

Setting the Rooting State
-------------------------

All |Tree| objects have a boolean property, :attr:`~dendropy.dataobject.tree.Tree.is_rooted` that DendroPy uses to track whether or not the tree should be treated as rooted. The property :attr:`~dendropy.dataobject.tree.Tree.is_unrooted` is also defined, and these two properties are synchronized. Thus setting :attr:`~dendropy.dataobject.tree.Tree.is_rooted` to :keyword:`True` will result in :attr:`~dendropy.dataobject.tree.Tree.is_rooted` being set to :keyword:`False` and vice versa.

The state of a |Tree| object's rootedness flag does not modify any internal structural relationship between nodes. It simply determines how its splits hashes are calculated, which in turn affects a broad range of comparison and metric operations. Thus you probably want to update the splits hashes after modifying the :attr:`~dendropy.dataobject.tree.Tree.is_rooted` property by calling the :meth:`~dendropy.dataobject.tree.Tree.update_splits()`. Note that calling :meth:`~dendropy.dataobject.tree.Tree.update_splits()` on an unrooted tree will force the basal split to be a trifurcation. So if the original tree was bifurcating, the end result will be a tree with a trifurcation at the root. This can be prevented by passing in the keyword argument ``delete_outdegree_one=False`` to :meth:`~dendropy.dataobject.tree.Tree.update_splits()`.

::

    #! /usr/bin/env python

    import dendropy

    tree_str = "[&R] (A, (B, (C, (D, E))));"

    tree = dendropy.Tree.get_from_string(
            tree_str,
            "newick")

    print("Original:")
    print(tree.as_ascii_plot())

    tree.is_rooted = False
    print("After `is_rooted=False`:")
    print(tree.as_ascii_plot())

    tree.update_splits()
    print("After `update_splits()`:")
    print(tree.as_ascii_plot())

    tree2 = dendropy.Tree.get_from_string(
            tree_str,
            "newick")
    tree2.is_rooted = False
    tree2.update_splits(delete_outdegree_one=False)
    print("After `update_splits(delete_outdegree_one=False)`:")
    print(tree2.as_ascii_plot())

will result in::


    Original:
    /---------------------------------------------------- A
    +
    |            /--------------------------------------- B
    \------------+
                 |            /-------------------------- C
                 \------------+
                              |            /------------- D
                              \------------+
                                           \------------- E


    After `is_rooted=False`:
    /---------------------------------------------------- A
    +
    |            /--------------------------------------- B
    \------------+
                 |            /-------------------------- C
                 \------------+
                              |            /------------- D
                              \------------+
                                           \------------- E


    After `update_splits()`:
    /---------------------------------------------------- A
    |
    +---------------------------------------------------- B
    |
    |                /----------------------------------- C
    \----------------+
                     |                 /----------------- D
                     \-----------------+
                                       \----------------- E


    After `update_splits(delete_outdegree_one=False)`:
    /---------------------------------------------------- A
    +
    |            /--------------------------------------- B
    \------------+
                 |            /-------------------------- C
                 \------------+
                              |            /------------- D
                              \------------+
                                           \------------- E

Derooting
---------

To deroot a rooted |Tree|, you can also call the :meth:`~dendropy.dataobject.tree.Tree.deroot()` method, which collapses the root to a trifurcation if it is bifurcation *and* sets the :attr:`~dendropy.dataobject.tree.Tree.is_rooted` to :keyword:`False`. The :meth:`~dendropy.dataobject.tree.Tree.deroot()` method has the same structural and semantic affect of :attr:`~dendropy.dataobject.tree.Tree.is_rooted` to :keyword:`False` and then calling :meth:`~dendropy.dataobject.tree.Tree.update_splits()`. You would use the former if you are *not* going to be doing any tree comparisons or calculating tree metrics, and thus do not want to calculate the splits hashes.

Rerooting
---------

To reroot a |Tree| along an existing edge, you can use the :meth:`~dendropy.dataobject.tree.Tree.reroot_at_node()` method. This method takes an |Edge| object as as its first argument. This rerooting is a structural change that will require the splits hashes to be updated before performing any tree comparisons or calculating tree metrics. You can do this yourself by calling :meth:`~dendropy.dataobject.tree.Tree.update_splits()` later, or you can pass in :keyword:`True` as the second argument to the  :meth:`~dendropy.dataobject.tree.Tree.reroot_at_edge()` method call, which instructs DendroPy to automatically update the splits for you. For example, the following reroots the tree along an internal edge:

.. literalinclude:: /examples/reroot_at_internal_edge.py

and results in::

    Before:
    [&R] (A,(B,(C,(D,E))));

    /---------------------------------------------------- A
    +
    |            /--------------------------------------- B
    \------------+
                 |            /-------------------------- C
                 \------------+
                              |            /------------- D
                              \------------+
                                           \------------- E


    After:
    [&R] ((D,E),(C,(B,A)));

                                       /----------------- D
    /----------------------------------+
    |                                  \----------------- E
    +
    |                /----------------------------------- C
    \----------------+
                     |                 /----------------- B
                     \-----------------+
                                       \----------------- A

Another example, this time rerooting along an edge subtending a tip instead of an internal edge:

.. literalinclude:: /examples/reroot_at_external_edge.py

which results in::

    Before:
    [&R] (A,(B,(C,(D,E))));

    /---------------------------------------------------- A
    +
    |            /--------------------------------------- B
    \------------+
                |             /-------------------------- C
                 \------------+
                              |            /------------- D
                              \------------+
                                           \------------- E


    After:
    [&R] (D,(E,(C,(B,A))));

    /---------------------------------------------------- D
    +
    |            /--------------------------------------- E
    \------------+
                 |            /-------------------------- C
                 \------------+
                              |            /------------- B
                              \------------+
                                           \------------- A

To reroot a |Tree| at a node instead, you can use the :meth:`~dendropy.dataobject.tree.Tree.reroot_at_node()` method:

.. literalinclude:: /examples/reroot_at_node.py

which results in::

    Before:
    [&R] (A,(B,(C,(D,E))));

    /---------------------------------------------------- A
    +
    |            /--------------------------------------- B
    \------------+
                 |            /-------------------------- C
                 \------------+
                              |            /------------- D
                              \------------+
                                           \------------- E


    After:
    [&R] (D,E,(C,(B,A)));

    /---------------------------------------------------- D
    |
    +---------------------------------------------------- E
    |
    |                /----------------------------------- C
    \----------------+
                     |                 /----------------- B
                     \-----------------+
                                       \----------------- A


You can also reroot the tree such that a particular node is moved to the outgroup position using the :meth:`~dendropy.dataobject.tree.Tree.to_outgroup_position()`, which takes a |Node| as the first argument. Again, you can update the splits hashes *in situ* by passing :keyword:`True` to the second argument. For example:

.. literalinclude:: /examples/to_outgroup_position.py

which will result in::

    Before:
    [&R] (A,(B,(C,(D,E))));

    /---------------------------------------------------- A
    +
    |            /--------------------------------------- B
    \------------+
                 |            /-------------------------- C
                 \------------+
                              |            /------------- D
                              \------------+
                                           \------------- E


    After:
    [&R] (C,(D,E),(B,A));

    /---------------------------------------------------- C
    |
    |                         /-------------------------- D
    +-------------------------+
    |                         \-------------------------- E
    |
    |                         /-------------------------- B
    \-------------------------+
                               \-------------------------- A

