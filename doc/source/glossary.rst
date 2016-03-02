*************************************
Glossary and Terminological Reference
*************************************

.. glossary::
    :sorted:

    ancestor nodes
        On a :term:`tree`, the full set of :term:`nodes <node>` from which a
        given :term:`node <node>` has descended (i.e., the given node's parent,
        parent's parent, parent's parent's parent, and so on until the
        :term:`root <seed node>` or the :term:`seed node` of the tree, inclusive).

    basal
    basal node
        On a :term:`tree`, a :term:`node` that is ancestral with respect to
        another node. :term:`Leaf or tip nodes <leaf node>` can never be basal,
        except if the tree consists of only a single node (i.e., the
        :term:`seed node` is a :term:`leaf node`). The basal-most node of a
        tree is the :term:`seed or root node <seed node>`.

    basal bifurcation
        On a :term:`tree`, a :term:`seed or root <seed node>` which has exactly
        two :term:`child nodes <child node>`. On a :term:`rooted tree <rooted
        tree>` this (can) reflect actual information, but on an :term:`unrooted
        tree <unrooted tree>`, this is actually an artifact, as the :term:`seed
        or root <seed node>` does not actually exist as it is just an
        algorithmic contrivance. In practical terms, this means that
        :term:`bipartition <bipartition>` calculations and operations on
        :term:`unrooted trees <unrooted tree>` with basal bifurcations may
        result in undetermined and errorenous behavior. Thus, typically, on
        unrooted trees the basal bifurcation is collapsed to yield a basal
        trifurcation.

    basal trifurcation
        On a :term:`tree`, a :term:`seed or root <seed node>` which has exactly
        three :term:`child nodes <child node>`. While this might occur in both
        :term:`rooted <rooted tree>` and :term:`unrooted <unrooted tree>` trees,
        this is typically the result of :term:`collapsing a basal bifurcation
        <basal bifurcation>` on :term:`unrooted trees <unrooted tree>`.

    bifurcation
    bifurcating node
        On a :term:`tree`, a :term:`node` with exactly two :term:`child nodes
        <child node>`. Also known as an "out-degree two" node.

    bifurcating tree
        A :term:`tree` in which all :term:`nodes <node>` are :term:`bifurcating
        <bifurcation>`.

    bipartition
    split
        On a :term:`tree`, a partition of the leaf set of a tree into two
        mutually-exclusive and collectively-comprehensive subsets. It
        corresponds to an edge of a tree: if we imagine "splitting" or cutting
        a tree into two trees at a given edge, the leaf sets of each of the new
        trees form the two subsets of the partitioning. A bipartition is often
        referred to as a split, especially in the context of :term:`unrooted trees
        <unrooted tree>`. Bipartitions are described in detail in
        the :doc:`DendroPy Primer <primer/bipartitions>`.

    child node
        On a :term:`tree`, a :term:`node` descending from another :term:`node`.
        A node on a tree may have zero or more child nodes. A node with zero
        child nodes is termed a :term:`leaf node`.

    descendent nodes
        On a :term:`tree`, the full set of nodes that are descended from a
        given node (i.e., the given node's children and children's children and
        children's children's children, and so on, until the :term:`leaf nodes
        <leaf node>` of the tree, inclusive).

    edge
    branch
        A connection between two :term:`nodes <node>` on a :term:`tree`,
        modeled in DendroPy by the |Edge| class. A synonym for "branch" in the
        context of phylogenetic trees.

        An edge connects a "tail node" (also called an origin or source node)
        to a "head node" (also called a target or destination node). The tail
        node is called the ":term:`parent <parent node>`" of the head node,
        while the head node is called the ":term:`child <child node>`" of the
        tail node.

        An edge is said to subtend or be incident to the node *to* which an
        edge connects, i.e., the head node.
        Conversely, the edges connecting a node to its children are called
        outgoing edges with respect to the tail node.

        On a :term:`tree`, every node has one and exactly one incident edge,
        and every edge has a :term:`head node`. On the other hand, not every
        node may have outgoing edges (e.g, :term:`leaf nodes <leaf node>`),
        and not every edge has a :term:`tail node` (e.g., :term:`root nodes
        <seed node>`). As such, edges can be thought of "belonging" to their
        head nodes, due to this one-to-one relationship.

        Edges can have any number of properties or attributes associated with
        them, representing a diverse range of phenomena, but the most important
        one is the edge :term:`length or weight <edge length>`.

    edge length
    edge weight
        A (typically) quantitative value associated with an :term:`edge`. This
        value may represent any number of things, but most typically is used to
        model time, evolutionary distance, or expected rates of substitution on
        a phylogenetic tree. An :term:`edge` may have many values, quantitative
        or otherwise, associated with it, but the length or weight is special
        as it usually denotes the relationship between the :term:`tail node`
        and :term:`head node` related by an :term:`edge`.

    incident edge
    subtending edge
        An :term:`edge` that connects *to* a particular :term:`node` is termed
        the incident or subtending edge of that node.

    internal edge
        An :term:`edge` that has an :term:`internal node` as a :term:`head
        node`.

    internal node
        A :term:`node` that has :term:`child nodes <child node>`. Also known as
        an inner node or branch node.

    head node
    target node
    destination node
        On an :term:`edge` connecting two :term:`nodes <node>`, the node *to*
        which the edge extends to link *from* the other node, termed the
        :term:`tail node`. The head node is the :term:`child node` of the
        :term:`tail node`,  and the :term:`tail node` is the :term:`parent
        node` of the head node. The :term:`edge` is said to subtend, or be
        incident, to the head node.

    leaf edge
    terminal edge
    external edge
    outer edge
        An :term:`edge` that has an :term:`leaf node` as a :term:`head
        node`.

    leaf node
    tip node
    terminal node
    external node
    outer node
        A :term:`node` that does not have any :term:`child nodes <child node>`
        descending from it. Also known as a tip, outer node, external node, or
        terminal node.

    node
        An fundamental element of information or data on a :term:`tree`,
        connected to other such elements in a parent-child relationshop by
        :term:`edges <edge>`.
        In DendroPy, a node is modeled by the |Node| class.
        A node has at most one :term:`parent <parent node>`, to which it is
        connected by its :term:`incident or subtending <incident edge>` edge.
        A node may have zero or more :term:`children <child node>`, to each of
        which it is connected by an independent :term:`outgoing edge <outgoing
        edge>` edge.
        A node can be associated with zero or more informational or data
        values. In a phylogenetic :term:`tree`, one of these values is often a
        :term:`taxon`, but many other aspects of information can be modeled.

    node depth
        On a :term:`tree`, the depth of a node is the length of the
        :term:`path` to its :term:`root <seed node>` (i.e., its root path). The
        root node has a depth zero.

    outgoing edge
        An :term:`edge` that connects *from* a particular :term:`node` (to,
        e.g., its :term:`children <child node>` is said to be an outgoing edge
        for that node.

    parent node
    ancestor node
        On a :term:`tree`, a :term:`node` from which a given node
        immediately descends.

    seed node
    root node
        The first or topmost :term:`node` in a tree. This is also more commonly
        called the "root" or "root node" of the tree, though, in the strictest
        sense, this equivalence is only valid when the tree is explicitly
        :term:`rooted <rooted tree>`. Both :term:`rooted trees <rooted tree>`
        and :term:`unrooted trees <unrooted tree>` have seed nodes. In rooted
        trees, the seed node is the root of the tree.

        By definition, the seed node does not have a :term:`parent node`.  It
        is the node at which algorithms on the tree begin, since as a data
        structure, one can only pass from :term:`parent node` to :term:`child
        nodes <child node>`.  If the tree is :term:`rooted <rooted tree>`, then
        the seed node is equivalent to the root of the tree.

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
        The sequence of :term:`edges <edge>` connecting two :term:`nodes
        <node>` on a :term:`tree`.

    path length
    path weight
    unweighted path
    weighted path
        The number or sum of lengths of the :term:`edges <edge>` connecting two
        :term:`nodes <node>` on a :term:`tree`. An *unweighted* path length is
        just the number of :term:`edges:, while a *weighted* path length or
        path weight is the sum of :term:`edge lengths <edge length>`.

    postorder traversal
        A process by which all nodes of a tree are visited, with a particular
        node being visited only after its children.

    preorder traversal
        A process by which all nodes of a tree are visited, with a particular
        node being visited before its children.

    rooted tree
        A state of a :term:`tree` in which its :term:`seed node` represents the
        most-recent common ancestor of all the :term:`leaf nodes <leaf node>`
        on the tree. Rooted trees have a distinct directionality, and
        ancestor-descendent relationships are not invertible.

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

    sibling node
    sister node
        Two or more :term:`nodes <node>` that share the same :term:`parent
        node`, i.e., are descended from the same node, are said to be siblings
        or sisters with respect to each other.

    subtree
        A subtree of a :term:`tree` is a tree consisting of a :term:`node` in
        the tree and all its :term:`descendents <descendent nodes>`. Each
        :term:`node` on a :term:`tree` thus corresponds to the :term:`root
        <seed node>` of the subtree that it determines.

    tail node
    source node
    origin node
        On an :term:`edge` connecting two :term:`nodes <node>`, the node *from*
        which the edge extends to link *to* the other node, termed the :term:`head node`.
        The tail node is the :term:`parent node` of the :term:`head node`,  and
        the :term:`head node` is the :term:`child node` of the tail node.
        The edge is said to be an outgoing node with respect to the tail node.

    taxon
    operational taxonomic unit
    taxon concept
    taxon name
        A class of organisms being modeled represented by a string label or
        namethat is guaranteed to be unique within a particular :term:`taxon
        namespace`.

    taxon namespace
        A set of distinct and unique labels, with each label mapping to one and
        exactly one or names that is used to relate data from across different
        data sources to each other by reference to a :term:`taxon concept`.

    tree
    arborescence
        A tree is a set of :term:`nodes <node>` connected to each other in
        parent-child relationships given by a set of :term:`edges <edge>`. In
        DendroPy, a tree is modeled by the |Tree| class. A tree is a
        specialization of a `graph
        <http://en.wikipedia.org/wiki/Graph_%28mathematics%29>`_, constrained
        such that:

            1. All its :term:`edges <edge>` are directional.
            2. It has no `directed cycles <http://en.wikipedia.org/wiki/Cycle_graph#Directed_cycle_graph>`_ .
            3. The directionality is from the :term:`root <seed node>` (or
               ":term:`seed node`" in DendroPy's parlance) to the
               :term:`tips <leaf node>`.

        The first and second constraints alone result in a `directed acylic graph
        <http://en.wikipedia.org/wiki/Directed_acyclic_graph>`_ .
        The addition of the third constraint results in an `arborescence
        <http://en.wikipedia.org/wiki/Arborescence_(graph_theory)>`_, which is
        strictly synonymous with "tree".

    unifurcation
    unifurcating node
        On a :term:`tree`, a :term:`node` with exactly one :term:`child node`.
        Also known as an "out-degree one" node. In some cases, unifurcations
        may be used to represent information (e.g., a change in some value
        associated with edges, such as population size or a rate of some kind),
        but they more typically arise as side-effect of tree manipulation
        operations, such as re-rooting or pruning. Though DendroPy has no
        problem in handling unifurcations, trees with unifurcating nodes are
        considered pathological in many contexts and operations, and DendroPy
        thus provides facilities for suppressing unifurcations, either on
        existing trees or as they occur as a side-effect of other operations.

    unrooted tree
        A state of a :term:`tree` in which its :term:`seed node` is an
        algorithmic artifact, and not necessarily represents the most-recent
        common ancestor of all the :term:`leaf nodes <leaf node>` on the tree.
        In an unrooted trees, ancestor-descendent relationships are also
        algorithmic artifacts and can be (conceptually) inverted without
        changing the information represented by the tree, though this operation
        usually requires a fundamental restructuring of the computational
        representation of the tree.
