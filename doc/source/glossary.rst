Glossary and Terminological Reference
-------------------------------------

.. glossary::
    :sorted:

    ancestor nodes
        On a :term:`tree`, the full set of :term:`nodes <node>` from which a
        given :term:`nodes <node>` has descended (i.e., given node's parent,
        parent's parent, parent's parent's parent, and so on until the
        :term:`root` or the :term:`seed node` of the tree, inclusive).

    bipartition
    split
        On a :term:`tree`, a partition of the leaf set of a tree into two
        mutually-exclusive and collectively-comprehensive subsets. It
        corresponds to an edge of a tree: if we imagine "splitting" or cutting
        a tree into two trees at a given edge, the leaf sets of each of the new
        trees form the two subsets of the partitioning. A bipartition is often
        referred to as a split, especially in the context of :term:`unrooted`
        trees.

    child node
        On a :term:`tree`, a :term:`node` or the set of :term:`nodes <node>`
        that immediately descend from a given :term:`node`.

    descendent nodes
        The full set of nodes that are descended from a given node (i.e., the
        given node's children and children's children and children's children's
        children, and so on, until the :term:`leaf nodes <leaf node>` of the tree,
        inclusive).

    edge
    branch
        A connection between two :term:`nodes <node>` on a :term:`tree`,
        modeled in DendroPy by the |Edge| class. A synonym for "branch" in the
        context of phylogenetic trees. An edge connects a "tail node" (also
        called an origin or source node) to a "head node" (also called a target
        or destination node). The tail node is called the ":term:`parent
        <parent node>`" of the head node, while the head node is called the
        ":term:`child <child node>`" of the tail node. Edges can have any number
        of properties or attributes associated with them, representing a
        diverse range of phenomena, but the most important one is the edge
        :term:`length or weight <edge length>`.

    edge length
    edge weight
        A (typically) quantitative attribute of an :term:`edge`.

    node
        A node is a structure which may contain a value or condition, or
        represent a separate data structure (which could be a tree of its own).
        Each node in a tree has zero or more :term:`child nodes <child node>`,
        which are below it in the tree (by convention, trees are drawn growing
        downwards). A node that has a child is called the child's parent node
        (or ancestor node, or superior).  A node has at most one parent, to
        which it is connected by its subtending :term:`edge`.  A node has
        exactly one subtending edge, and this is typically accessed as an
        attribute of the node. A node may have zero or more outgoing edges,
        which connect it to its :term:`child nodes <child node>`.

    node depth
        On a :term:`tree`, the depth of a node is the length of the
        :term:`path` to its :term:`root` (i.e., its root path). The root node
        has a depth zero.

    internal node
        An internal node (also known as an inner node, inode for short, or branch
        node) is any :term:`node` of a tree that has :term:`child nodes <child node>`.

    leaf node
    tip node
    terminal node
    external node
    outer node
        An leaf :term:`node` (also known as a tip, outer node, external node, or
        terminal node) is any :term:`node` that does not have :term:`children
        <child node>`.

    parent node
    ancestor node
        On a :term:`tree`, a :term:`node` from which a given node
        immediately descends.

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
        nodes <child node>`.  If the tree is :term:`rooted`, then the seed node
        is equivalent to the root of the tree.

    node height
        The height of a node is the length of the longest downward path to a
        leaf from that node. The height of the root is the height of the tree.
        The depth of a node is the length of the path to its root (i.e., its
        root path). The root node has depth zero,
        leaf nodes have height zero, and a tree with only a single node
        (hence both a root and leaf) has depth and height zero. Conventionally,
        an empty tree (tree with no nodes, if such are allowed) has depth and
        height 1.

    path
    path length
    path weight
    unweighted path
    weighted path
        In the context of :term:`trees <tree>`, the number or sum of lengths of
        :term:`edges <edge>` connecting two :term:`nodes <node>`. An
        *unweighted* path length is just the number of :term:`edges:, while a
        *weighted* path length or path weight is the sum of :term:`edge lengths
        <edge length>`.

    rooted
        A state of a :term:`tree` in which its :term:`seed node` represents the
        most-recent common ancestor of all the :term:`leaf nodes <leaf node>` on the tree.

    schema
        The format or syntax of serialized phylogenetic or related data.
        Examples are NEXUS, NEWICk, Phylip, NeXML, etc. A "schema" is
        DendroPy-speak for "format" (we cannot use the argument name "format"
        because this is a Python built-in, and hence we adopted this
        terminology for consistency), and is typicallly specified using one of
        a set of predefined string values, known as "schema specification
        strings". Supported reading (input) schemas are described :ref:`here
        <Specifying_the_Data_Source_Format>` while supported writing (output)
        schemas are described :ref:`here <Specifying_the_Data_Writing_Format>`.


    subtree
        A subtree of a tree T is a tree consisting of a node in T and all of
        its descendants in T.[c][1] Nodes thus correspond to subtrees (each
        node corresponds to the subtree of itself and all its descendants)  the
        subtree corresponding to the root node is the entire tree, and each
        node is the root node of the subtree it determines; the subtree
        corresponding to any other node is called a proper subtree (in analogy
        to the term proper subset).

    tree
        An `arborescence
        <http://en.wikipedia.org/wiki/Arborescence_(graph_theory)>`_, or a
        fully-connected `directed acylic graph
        <http://en.wikipedia.org/wiki/Directed_acyclic_graph>`_ in which the
        directionality is from the :term:`root` (or ":term:`seed node`" in
        DendroPy's parlance) in which the direction to the :term:`tips <leaf node>`.

    unrooted
        A state of a :term:`tree` in which its :term:`seed node` is an
        algorithmic artifact, and not necessarily represents the most-recent
        common ancestor of all the :term:`leaf nodes <leaf node>` on the tree.
