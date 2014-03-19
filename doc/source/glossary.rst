Tree Terminology
----------------

.. glossary::
    :sorted:

    node
        A node is a structure which may contain a value or condition, or represent
        a separate data structure (which could be a tree of its own). Each node in
        a tree has zero or more child nodes, which are below it in the tree (by
        convention, trees are drawn growing downwards). A node that has a
        child is called the child's parent node (or ancestor node, or superior).
        A node has at most one parent.

    parent node
        A :term:`node` that is the immediate ancestor of a given node.

    child node
    child nodes
    child
    children
    children nodes
        A :term:`node` that is the immediate descendent of a given node.

    descendent nodes
        The full set of nodes that are descended from a given node (i.e., the
        given node's children and children's children and children's children's
        children, and so on, until the :term:`leaf` nodes of the tree,
        inclusive).

    ancestor nodes
        The full set of nodes from which a given node has descended (i.e.,
        given node's parent, parent's parent, parent's parent's parent, and so
        on until the :term:`root` or the :term:`seed node` of the tree,
        inclusive).

    internal node
        An internal node (also known as an inner node, inode for short, or branch
        node) is any :term:`node` of a tree that has child nodes.

    leaf node
    leaf
    tip
    terminal node
    external node
    outer node
        An leaf :term:`node` (also known as a tip, outer node, external node, or
        terminal node) is any :term:`node` that does not have :term:`children`.

    seed node
    root
        The first or topmost :term:`node` in a tree is called the seed node.
        This is also called the "root" or "root node" or the tree, though, in
        the strictest sense, this equivalence is only valid when the tree is
        explicitly :term:`rooted`. Both rooted trees and unrooted trees have
        seed nodes. In rooted trees, the seed node is the root of the tree.

        By definition, the seed node does not have a :term:`parent node`.  It
        is the node at which algorithms on the tree begin, since as a data
        structure, one can only pass from :term:`parent node` to :term:`child
        nodes`.  If the tree is :term:`rooted`, then the seed node is
        equivalent to the root of the tree.

    height
        The height of a node is the length of the longest downward path to a
        leaf from that node. The height of the root is the height of the tree.
        The depth of a node is the length of the path to its root (i.e., its
        root path). This is commonly needed in the manipulation of the various
        self-balancing trees, AVL Trees in particular. The root node has depth
        zero, leaf nodes have height zero, and a tree with only a single node
        (hence both a root and leaf) has depth and height zero. Conventionally,
        an empty tree (tree with no nodes, if such are allowed) has depth and
        height 1.

    subtree
        A subtree of a tree T is a tree consisting of a node in T and all of
        its descendants in T.[c][1] Nodes thus correspond to subtrees (each
        node corresponds to the subtree of itself and all its descendants)  the
        subtree corresponding to the root node is the entire tree, and each
        node is the root node of the subtree it determines; the subtree
        corresponding to any other node is called a proper subtree (in analogy
        to the term proper subset).
