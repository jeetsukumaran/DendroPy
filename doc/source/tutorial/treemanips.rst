***********************************
Tree Manipulation and Restructuring
***********************************

The |Tree| class provides both low-level and high-level methods for manipulating tree structure.

Low-level methods are associated with |Node| objects, and allow to restructure the relationships between nodes at a fine level: :meth:`~dendropy.dataobject.tree.Node.add_child`, :meth:`~dendropy.dataobject.tree.Node.new_child`, :meth:`~dendropy.dataobject.tree.Node.remove_child`, etc.

In most cases, however, you will be using high-level methods to restructure |Tree| objects.

In all cases, if any part of the |Tree| object's structural relations change, *and* you are interested in calculating any metrics or statistics on the tree or comparing the tree to another tree, you need to call :meth:`~dendropy.dataobject.tree.Tree.update_splits()` on the object to update the internal splits hash representation.

Rooting, Derooting and Rerooting
================================

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

To deroot a rooted |Tree|, you can also call the :meth:`~dendropy.dataobject.tree.Tree.deroot()` method, which collapses the root to a trifurcation if it is bifurcation *and* sets the :attr:`~dendropy.dataobject.tree.Tree.is_rooted` to :keyword:`False`. The :meth:`~dendropy.dataobject.tree.Tree.deroot()` method has the same structural and semantic affect of :attr:`~dendropy.dataobject.tree.Tree.is_rooted` to :keyword:`False` and then calling :meth:`~dendropy.dataobject.tree.Tree.update_splits()`. You would use the former if you are *not* going to be doing any tree comparisons or calculating tree metrics, and thus do not want to calculate the splits hashes.

To *reroot* a |Tree| at a given node, you can use the :meth:`~dendropy.dataobject.tree.Tree.reroot_at()` method. This method takes a |Node| object as an argument::





