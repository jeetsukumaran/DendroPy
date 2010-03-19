***********************************
Tree Manipulation and Restructuring
***********************************

The :mod:`~dendropy.treemanip` module provides high-level functions for manipulating the structure of trees: :func:`~dendropy.treemanip.collapse_clade`, :func:`~dendropy.treemanip.collapse_edge`, :func:`~dendropy.treemanip.prune_subtree`, :func:`~dendropy.treemanip.prune_taxa`, :func:`~dendropy.treemanip.retain_taxa`,  :func:`~dendropy.treemanip.randomly_rotate`, etc.

In addition, the native methods of a |Node| object of a |Tree|, such as :meth:`~dendropy.dataobject.tree.Node.add_child`, :meth:`~dendropy.dataobject.tree.Node.new_child`, :meth:`~dendropy.dataobject.tree.Node.remove_child`, etc. can be used for low-level manipulation of the tree structure.

In all cases, if any the structure of the tree changes in any way, it is important to call the native |Tree| method, :meth:`~dendropy.dataobject.tree.Tree.update_splits()`, before calculating any metrics or statistics on the tree.


MORE COMING SOON ...
