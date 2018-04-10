***********************************
Tree Manipulation and Restructuring
***********************************

The |Tree| class provides both low-level and high-level methods for manipulating tree structure.

Low-level methods are associated with |Node| objects, and allow to restructure the relationships between nodes at a fine level: :meth:`~dendropy.datamodel.treemodel.Node.add_child`, :meth:`~dendropy.datamodel.treemodel.Node.new_child`, :meth:`~dendropy.datamodel.treemodel.Node.remove_child`, etc.

In most cases, however, you will be using high-level methods to restructure |Tree| objects.

In all cases, if any part of the |Tree| object's structural relations change, *and* you are interested in calculating any metrics or statistics on the tree or comparing the tree to another tree, you need to call :meth:`~dendropy.datamodel.treemodel.Tree.update_bipartitions()` on the object to update the internal splits hash representation.
This is not done for you automatically because there is a computational cost associated with the operation, and the splits hashes are not always needed. Furthermore, even when needed, if there are a number of structural changes to be made to a |Tree| object before calculations/comparisions, it makes sense to postpone the splits rehashing until there all the tree manipulations are completed.
Most methods that affect the tree structure that require the splits hashes to updated take a ``update_bipartitions`` argument. By specifying |True| for this, the |Tree| object will recalculate the splits hashes after the changes have been made.

Rooting, Derooting and Rerooting
================================

The Rooting of Tree(s) Read from External Sources
-------------------------------------------------

|Newick| and |Nexus| formats have a convention where the rooting of the tree is specified by a special comment token preceding the tree statement: "``[&R]``" to indicate a rooted tree::

    [&R] ((A,B),(C,D));

and : "``[&U]``" to indicate an unrooted tree::

    [&U] ((A,B),(C,D));

These rooting comment tokens are respected when tree data is read. If no such comment token is given, then the tree is assumed to be **unrooted** by default.

You can control the behavior of trees read by using the "``rooting`` keyword argument when using the "|get|" or "|read|" methods of the |Tree| or |TreeList| classes. This takes one of four string values which determines how the rooting states of the tree(s) will be set:

    "default-unrooted" [default]
        All trees are interpreted as unrooted unless a "[&R]"
        comment token explicitly specifies them as rooted.
    "default-rooted"
        All trees are interpreted as rooted unless a "[&U]"
        comment token explicitly specifies them as unrooted.
    "force-unrooted"
        All trees are unconditionally interpreted as unrooted.
    "force-rooted"
        All trees are unconditionally interpreted as rooted.

The behavior of this can be be summarized by the following:

+-----------------------------------+----------------------------+----------------------------+------------------------------------+
| Keyword Argument                  | ``[&U]`` in Tree Statement | ``[&R]`` in Tree Statement | No Rooting Given in Tree Statement |
+-----------------------------------+----------------------------+----------------------------+------------------------------------+
| ``rooting=None`` (or unspecified) | unrooted                   | rooted                     | None                               |
+-----------------------------------+----------------------------+----------------------------+------------------------------------+
| ``rooting="default-unrooted"``    | unrooted                   | rooted                     | unrooted                           |
+-----------------------------------+----------------------------+----------------------------+------------------------------------+
| ``rooting="default-rooted"``      | unrooted                   | rooted                     | rooted                             |
+-----------------------------------+----------------------------+----------------------------+------------------------------------+
| ``rooting="force-unrooted"``      | unrooted                   | unrooted                   | unrooted                           |
+-----------------------------------+----------------------------+----------------------------+------------------------------------+
| ``rooting="force-rooted"``        | rooted                     | rooted                     | rooted                             |
+-----------------------------------+----------------------------+----------------------------+------------------------------------+

As an example:

.. literalinclude:: /examples/readroot2.py

Setting the Rooting State
-------------------------

All |Tree| objects have a boolean property, :attr:`~dendropy.datamodel.treemodel.Tree.is_rooted` that DendroPy uses to track whether or not the tree should be treated as rooted. The property :attr:`~dendropy.datamodel.treemodel.Tree.is_unrooted` is also defined, and these two properties are synchronized. Thus setting :attr:`~dendropy.datamodel.treemodel.Tree.is_rooted` to |True| will result in :attr:`~dendropy.datamodel.treemodel.Tree.is_rooted` being set to |False| and vice versa.

The state of a |Tree| object's rootedness flag does not modify any internal structural relationship between nodes. It simply determines how its splits hashes are calculated, which in turn affects a broad range of comparison and metric operations. Thus you need to update the splits hashes after modifying the :attr:`~dendropy.datamodel.treemodel.Tree.is_rooted` property by calling the :meth:`~dendropy.datamodel.treemodel.Tree.update_bipartitions()` before carrying out any calculations on or with the |Tree| object. Note that calling :meth:`~dendropy.datamodel.treemodel.Tree.update_bipartitions()` on an unrooted tree will force the basal split to be a trifurcation. So if the original tree was bifurcating, the end result will be a tree with a trifurcation at the root. This can be prevented by passing in the keyword argument ``suppress_unifurcations=False`` to :meth:`~dendropy.datamodel.treemodel.Tree.update_bipartitions()`.

For example, the following:

.. literalinclude:: /examples/setroot1.py


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


    After `update_bipartitions()`:
    /---------------------------------------------------- A
    |
    +---------------------------------------------------- B
    |
    |                /----------------------------------- C
    \----------------+
                     |                 /----------------- D
                     \-----------------+
                                       \----------------- E


    After `update_bipartitions(suppress_unifurcations=False)`:
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

To deroot a rooted |Tree|, you can also call the :meth:`~dendropy.datamodel.treemodel.Tree.deroot()` method, which collapses the root to a trifurcation if it is bifurcation *and* sets the :attr:`~dendropy.datamodel.treemodel.Tree.is_rooted` to |False|. The :meth:`~dendropy.datamodel.treemodel.Tree.deroot()` method has the same structural and semantic affect of :attr:`~dendropy.datamodel.treemodel.Tree.is_rooted` to |False| and then calling :meth:`~dendropy.datamodel.treemodel.Tree.update_bipartitions()`. You would use the former if you are *not* going to be doing any tree comparisons or calculating tree metrics, and thus do not want to calculate the splits hashes.

Rerooting
---------

To reroot a |Tree| along an existing edge, you can use the :meth:`~dendropy.datamodel.treemodel.Tree.reroot_at_edge()` method. This method takes an |Edge| object as as its first argument. This rerooting is a structural change that will require the splits hashes to be updated before performing any tree comparisons or calculating tree metrics. If needed, you can do this yourself by calling :meth:`~dendropy.datamodel.treemodel.Tree.update_bipartitions()` later, or you can pass in |True| as the second argument to the  :meth:`~dendropy.datamodel.treemodel.Tree.reroot_at_edge()` method call, which instructs DendroPy to automatically update the splits for you.

As an example, the following reroots the tree along an internal edge (note that we do not recalculate the splits hashes, as we are not carrying out any calculations or comparisons with the |Tree|):

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

To reroot a |Tree| at a node instead, you can use the :meth:`~dendropy.datamodel.treemodel.Tree.reroot_at_node()` method:

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


You can also reroot the tree such that a particular node is moved to the outgroup position using the :meth:`~dendropy.datamodel.treemodel.Tree.to_outgroup_position()`, which takes a |Node| as the first argument. Again, you can update the splits hashes *in situ* by passing |True| to the second argument, and again, here we do not because we are not carrying out any calculations. For example:

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

If you have a tree with edge lengths specified, you can reroot it at the midpoint, using the :meth:`~dendropy.datamodel.treemodel.Tree.reroot_at_midpoint()` method:

.. literalinclude:: /examples/reroot_at_midpoint.py

which results in::

    Before:
    [&R] (A:0.55,(B:0.82,(C:0.74,(D:0.42,E:0.64):0.24):0.15):0.2):0.3;

              /------------------- A
              +
              |      /---------------------------- B
              \------+
                     |    /-------------------------- C
                     \----+
                          |        /-------------- D
                          \--------+
                                   \---------------------- E


    After:
    [&R] ((C:0.74,(D:0.42,E:0.64):0.24):0.045,(B:0.82,A:0.75):0.105):0.3;

                   /------------------------------- C
                 /-+
                 | |         /------------------ D
                 | \---------+
                 +           \---------------------------- E
                 |
                 |   /------------------------------------ B
                 \---+
                     \-------------------------------- A



Pruning Subtrees and Tips
=========================

To remove a set of tips from a |Tree|, you cna use either the :meth:`~dendropy.datamodel.treemodel.Tree.prune_taxa()` or the :meth:`~dendropy.datamodel.treemodel.Tree.prune_taxa_with_labels()` methods. The first takes a container of |TaxonNamespace| objects as an argument, while the second takes container of strings. In both cases, nodes associated with the specified taxa (as given by the |TaxonNamespace| objects directly in the first case, or |TaxonNamespace| objects with labels given in the list of string in the second case) will e removed from the tree. For example:

.. literalinclude:: /examples/prune_taxa_with_labels.py

which results in::

    Before:
    [&R] ((A,(B,(C,(D,E)))),(F,(G,H)));

              /------------------------------------------- A
    /---------+
    |         |          /-------------------------------- B
    |         \----------+
    |                    |          /--------------------- C
    |                    \----------+
    +                               |          /---------- D
    |                               \----------+
    |                                          \---------- E
    |
    |                               /--------------------- F
    \-------------------------------+
                                    |          /---------- G
                                    \----------+
                                               \---------- H


    After:
    [&R] ((B,(D,E)),(F,H));

                      /----------------------------------- B
    /-----------------+
    |                 |                 /----------------- D
    |                 \-----------------+
    +                                   \----------------- E
    |
    |                                   /----------------- F
    \-----------------------------------+
                                        \----------------- H

Alternatively, the tree can be pruned based on a set of taxa that you want to *keep*. This can be affected through the use of the counterpart "retain" methods, :meth:`~dendropy.datamodel.treemodel.Tree.retain_taxa()` and :meth:`~dendropy.datamodel.treemodel.Tree.retain_taxa_with_labels()`. For example:

.. literalinclude:: /examples/retain_taxa_with_labels.py

which results in::

    Before:
    [&R] ((A,(B,(C,(D,E)))),(F,(G,H)));

              /------------------------------------------- A
    /---------+
    |         |          /-------------------------------- B
    |         \----------+
    |                    |          /--------------------- C
    |                    \----------+
    +                               |          /---------- D
    |                               \----------+
    |                                          \---------- E
    |
    |                               /--------------------- F
    \-------------------------------+
                                    |          /---------- G
                                    \----------+
                                               \---------- H


    After:
    [&R] ((A,C),G);

                               /-------------------------- A
    /--------------------------+
    +                          \-------------------------- C
    |
    \----------------------------------------------------- G

Again, it should be noted that, as these operations modify the structure of the tree, you need to call :meth:`~dendropy.datamodel.treemodel.Tree.update_bipartitions()` to update the internal splits hashes, before carrying out any calculations, comparisons, or metrics.

Extracting Trees and Subtrees from an Existing Tree
===================================================

If you just need a tree *structure*, i.e. the nodes and edges, and minimal other attributes, such as taxon associations, node and edge labels, and edge lengths, the following methods provide fast and powerful ways to give you copies of the current tree or parts of the current tree:

-   :meth:`~dendropy.datamodel.treemodel.Tree.extract_tree`
-   :meth:`~dendropy.datamodel.treemodel.Tree.extract_tree_with_taxa`
-   :meth:`~dendropy.datamodel.treemodel.Tree.extract_tree_without_taxa`
-   :meth:`~dendropy.datamodel.treemodel.Tree.extract_tree_with_taxa_labels`
-   :meth:`~dendropy.datamodel.treemodel.Tree.extract_tree_without_taxa_labels`

All these methods create a "thin" or "light" copy of the tree, optionally pruning out nodes based on criteria.

The basic :meth:`~dendropy.datamodel.treemodel.Tree.extract_tree` method returns a full clone of the structure of the source tree:

.. literalinclude:: /examples/tree_extract1.py

As can be seen, the nodes on the extracted tree have an attribute set on them, "``extraction_source``" that point to the corresponding node on the original tree.
Note how while labels, taxon associations and edge lengths are copied, the metadata annotations are not.

The ``node_filter_fn`` argument provides a powerful and flexible way to selective pull out parts of the tree that might interest you:

.. literalinclude:: /examples/tree_extract2.py

If you do not need the annotations, then this approach can be dramatically more performant than cloning and then pruning the tree. For example, the following:

.. literalinclude:: /examples/tree_clone_and_prune_vs_extract.py

results in::

    0.00191879272461
    0.000579190254211

Some convenience wrappers around the :meth:`~dendropy.datamodel.treemodel.Tree.extract_tree` method allow you to more easily pull out or drop taxa of interest, either by passing in the |Taxon| instances directly or their labels:

.. literalinclude:: /examples/tree_extract3.py

You can also use the :meth:`~dendropy.datamodel.treemodel.Tree.extract_tree` method to created a "casted" copy, i.e. a (light) copy of the tree using your own custom classes instead of DendroPy's native |Tree| and |Node| classes:

.. literalinclude:: /examples/tree_extract4.py

Rotating
========

You can ladderize trees (sort the child nodes in order of the number of their children) by calling the :meth:`~dendropy.datamodel.treemodel.Tree.ladderize()` method. This method takes one argument, ``ascending``. If ``ascending=True``, which is the default, then the nodes are sorted in ascending order (i.e., nodes with fewer children sort before nodes with more children). If ``ascending=False``, then the nodes are sorted in descending order (i.e., nodes with more children sorting before nodes with fewer children). For example:

.. literalinclude:: /examples/ladderize.py


results in::

    Before:
    [&R] ((A,(B,(C,(D,E)))),(F,(G,H)));

              /------------------------------------------- A
    /---------+
    |         |          /-------------------------------- B
    |         \----------+
    |                    |          /--------------------- C
    |                    \----------+
    +                               |          /---------- D
    |                               \----------+
    |                                          \---------- E
    |
    |                               /--------------------- F
    \-------------------------------+
                                    |          /---------- G
                                    \----------+
                                               \---------- H


    Ladderize, ascending=True:
    [&R] ((F,(G,H)),(A,(B,(C,(D,E)))));

                                    /--------------------- F
    /-------------------------------+
    |                               |          /---------- G
    |                               \----------+
    +                                          \---------- H
    |
    |         /------------------------------------------- A
    \---------+
              |          /-------------------------------- B
              \----------+
                         |          /--------------------- C
                         \----------+
                                    |          /---------- D
                                    \----------+
                                               \---------- E


    Ladderize, ascending=False:
    [&R] (((((D,E),C),B),A),((G,H),F));

                                               /---------- D
                                    /----------+
                         /----------+          \---------- E
                         |          |
              /----------+          \--------------------- C
              |          |
    /---------+          \-------------------------------- B
    |         |
    |         \------------------------------------------- A
    +
    |                                          /---------- G
    |                               /----------+
    \-------------------------------+          \---------- H
                                    |
                                    \--------------------- F

Tree rotation operations do not actually change the tree structure, at least in so far as splits are concerned, so it is not neccessary to update the splits hashes.
