#! /usr/bin/env python
# -*- coding: utf-8 -*-

from io import StringIO
from dendropy.utility import deprecate
from dendropy.utility import error
from dendropy.datamodel import basemodel
from dendropy.datamodel.treemodel import _edge

class Node(basemodel.DataObject, basemodel.Annotable):
    """
    A :term:|Node| on a :term:|Tree|.
    """

    @classmethod
    def edge_factory(cls, **kwargs):
        """
        Creates and returns a |Edge| object.

        Derived classes can override this method to provide support for
        specialized or different types of edges on the tree.

        Parameters
        ----------

        kwargs : keyword arguments
            Passed directly to constructor of |Edge|.

        Returns
        -------
        |Edge|
            A new |Edge| object.

        """
        return _edge.Edge(**kwargs)

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------
        taxon : |Taxon|, optional
            The |Taxon| instance representing the operational taxonomic
            unit concept associated with this Node.
        label : string, optional
            A label for this node.
        edge_length : numeric, optional
            Length or weight of the edge subtending this node.

        """
        basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
        self.taxon = kwargs.pop("taxon", None)
        self.age = None
        self._edge = None
        self._child_nodes = []
        self._parent_node = None
        self.edge = self.edge_factory(
            head_node=self, length=kwargs.pop("edge_length", None)
        )
        if kwargs:
            raise TypeError("Unsupported keyword arguments: {}".format(kwargs))
        self.comments = []

    def __copy__(self, memo=None):
        raise TypeError("Cannot directly copy Edge")

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot directly copy Node")

    def __deepcopy__(self, memo=None):
        return basemodel.Annotable.__deepcopy__(self, memo=memo)
        # if memo is None:
        #     memo = {}
        # other = basemodel.Annotable.__deepcopy__(self, memo=memo)
        # memo[id(self._child_nodes)] = other._child_nodes
        # for ch in self._child_nodes:
        #     try:
        #         och = memo[id(ch)]
        #         if och not in other._child_nodes:
        #             other._child_nodes.append(och)
        #     except KeyError:
        #         och = copy.deepcopy(ch, memo)
        #         memo[id(chd)] = och
        #         if och not in other._child_nodes:
        #             other._child_nodes.append(och)
        # return other
        # return super(Node, self).__deepcopy__(memo=memo)

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        # IMPORTANT LESSON LEARNED: if you define __hash__, you *must* define __eq__
        return self is other

    def __repr__(self):
        return "<{} object at {}: '{}' ({})>".format(
            self.__class__.__name__, hex(id(self)), self._label, repr(self.taxon)
        )

    def __iter__(self, *args, **kwargs):
        return self.preorder_iter(*args, **kwargs)

    def preorder_iter(self, filter_fn=None):
        """
        Pre-order iterator over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited before its
        children. Nodes can optionally be filtered by ``filter_fn``: only nodes
        for which ``filter_fn`` returns |True| when called with the node as an
        argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes of the subtree rooted at this node in
            pre-order sequence.
        """
        stack = [self]
        while stack:
            node = stack.pop()
            if filter_fn is None or filter_fn(node):
                yield node
            stack.extend(n for n in reversed(node._child_nodes))

    def preorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes of subtree rooted at this node.

        Visits self and all internal descendant nodes, with each node visited
        before its children. In DendroPy, "internal nodes" are nodes that have
        at least one child node, and thus the root or seed node is typically included
        unless ``exclude_seed_node`` is |True|. Nodes can optionally be filtered
        by ``filter_fn``: only nodes for which ``filter_fn`` returns |True| when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If |False| (default), then the seed node or root is visited. If
            |True|, then the seed node is skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the internal nodes of the subtree rooted at
            this node in pre-order sequence.
        """
        if exclude_seed_node:
            froot = lambda x: x._parent_node is not None
        else:
            froot = lambda x: True
        if filter_fn:
            f = lambda x: (froot(x) and x._child_nodes and filter_fn(x)) or None
        else:
            f = lambda x: (x and froot(x) and x._child_nodes) or None
        return self.preorder_iter(filter_fn=f)

    def postorder_iter(self, filter_fn=None):
        """
        Post-order iterator over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited after its
        children. Nodes can optionally be filtered by ``filter_fn``: only nodes
        for which ``filter_fn`` returns |True| when called with the node as an
        argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the nodes of the subtree rooted at
            this node in post-order sequence.
        """
        # if self._child_nodes:
        #     for nd in self._child_nodes:
        #         for ch in nd.postorder_iter(filter_fn=filter_fn):
        #             yield ch
        # if filter_fn is None or filter_fn(self):
        #     yield self
        # return

        # stack = [(self, False)]
        # while stack:
        #     node, state = stack.pop(0)
        #     if state:
        #         if filter_fn is None or filter_fn(node):
        #             yield node
        #     else:
        #         stack.insert(0, (node, True))
        #         child_nodes = [(n, False) for n in node._child_nodes]
        #         child_nodes.extend(stack)
        #         stack = child_nodes

        ## Prefer `pop()` to `pop(0)`.
        ## Thanks to Mark T. Holder
        ## From peyotl commits: d1ffef2 + 19fdea1
        stack = [(self, False)]
        while stack:
            node, state = stack.pop()
            if state:
                if filter_fn is None or filter_fn(node):
                    yield node
            else:
                stack.append((node, True))
                stack.extend([(n, False) for n in reversed(node._child_nodes)])

    def postorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes of subtree rooted at this node.

        Visits self and all internal descendant nodes, with each node visited
        after its children. In DendroPy, "internal nodes" are nodes that have
        at least one child node, and thus the root or seed node is typically
        included unless ``exclude_seed_node`` is |True|. Nodes can optionally be
        filtered by ``filter_fn``: only nodes for which ``filter_fn`` returns
        |True| when passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If |False| (default), then the seed node or root is visited. If
            |True|, then the seed node is skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding the internal nodes of the subtree rooted at
            this node in post-order sequence.
        """
        if exclude_seed_node:
            froot = lambda x: x._parent_node is not None
        else:
            froot = lambda x: True
        if filter_fn:
            f = lambda x: (froot(x) and x._child_nodes and filter_fn(x)) or None
        else:
            f = lambda x: (x and froot(x) and x._child_nodes) or None
        return self.postorder_iter(filter_fn=f)

    def levelorder_iter(self, filter_fn=None):
        """
        Level-order iteration over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node and other nodes at
        the same level (distance from root) visited before their children.
        Nodes can optionally be filtered by ``filter_fn``: only nodes for which
        ``filter_fn`` returns |True| when called with the node as an argument are
        visited.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes of the subtree rooted at this node in
            level-order sequence.
        """
        if filter_fn is None or filter_fn(self):
            yield self
        remaining = self.child_nodes()
        while len(remaining) > 0:
            node = remaining.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = node.child_nodes()
            remaining.extend(child_nodes)

    def level_order_iter(self, filter_fn=None):
        """
        DEPRECATED: Use :meth:`Node.levelorder_iter()` instead.
        """
        deprecate.dendropy_deprecation_warning(
            message=(
                "Deprecated since DendroPy 4: 'level_order_iter()' will no longer be"
                " supported in future releases; use 'levelorder_iter()' instead"
            ),
            stacklevel=3,
        )
        return self.levelorder_iter(filter_fn=filter_fn)

    def inorder_iter(self, filter_fn=None):
        """
        In-order iteration over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited in-between
        its children. Only valid for strictly-bifurcating trees. Nodes can
        optionally be filtered by ``filter_fn``: only nodes for which ``filter_fn``
        returns |True| when called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes of the subtree rooted at this node in
            infix or in-order sequence.
        """
        if len(self._child_nodes) == 0:
            if filter_fn is None or filter_fn(self):
                yield self
        elif len(self._child_nodes) == 2:
            for nd in self._child_nodes[0].inorder_iter(filter_fn=filter_fn):
                yield nd
            if filter_fn is None or filter_fn(self):
                yield self
            for nd in self._child_nodes[1].inorder_iter(filter_fn=filter_fn):
                yield nd
        else:
            raise TypeError("In-order traversal only supported for binary trees")

    def leaf_iter(self, filter_fn=None):
        """
        Iterate over all tips or leaves that ultimately descend from this node.

        Visits all leaf or tip nodes descended from this node. Nodes can
        optionally be filtered by ``filter_fn``: only nodes for which ``filter_fn``
        returns |True| when called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding leaf nodes of the subtree rooted at this node.
        """
        if filter_fn:
            ff = lambda x: x.is_leaf() and filter_fn(x) or None
        else:
            ff = lambda x: x.is_leaf() and x or None
        for node in self.postorder_iter(ff):
            yield node

    def child_node_iter(self, filter_fn=None):
        """
        Iterator over all nodes that are the (immediate) children of this node.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            An iterator yielding nodes that have this node as a parent.
        """
        for node in self._child_nodes:
            if filter_fn is None or filter_fn(node):
                yield node

    def child_edge_iter(self, filter_fn=None):
        """
        Iterator over all edges that are the (immediate) children of this edge.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Edge| object as an argument
            and returns |True| if the |Edge| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all edges visited will be yielded.

        Returns
        -------
        :py:class:`collections.Iterator` [|Edge|]
            An iterator yielding edges that have this edge as a parent.
        """
        for node in self._child_nodes:
            if filter_fn is None or filter_fn(node.edge):
                yield node.edge

    def ancestor_iter(self, filter_fn=None, inclusive=False):
        """
        Iterator over all ancestors of this node.

        Visits all nodes that are the ancestors of this node.  If ``inclusive``
        is |True|, ``self`` is returned as the first item of the sequence;
        otherwise ``self`` is skipped. Nodes can optionally be filtered by
        ``filter_fn``: only nodes for which ``filter_fn`` returns |True| when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (default), then all nodes visited will be yielded.
        inclusive : boolean, optional
            If |True|, includes this node in the sequence. If |False|, this is
            skipped.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            Iterator over all predecessor/ancestor nodes of this node.
        """
        if inclusive and (filter_fn is None or filter_fn(self)):
            yield self
        node = self
        while node is not None:
            node = node._parent_node
            if node is not None and (filter_fn is None or filter_fn(node)):
                yield node

    def ageorder_iter(self, filter_fn=None, include_leaves=True, descending=False):
        """
        Iterator over nodes of subtree rooted at this node in order of the age
        of the node (i.e., the time since the present).

        Iterates over nodes in order of age ('age' is as given by the ``age``
        attribute, which is usually the sum of edge lengths from tips
        to node, i.e., time since present).
        If ``include_leaves`` is |True| (default), leaves are included in the
        iteration; if ``include_leaves`` is |False|, leaves will be skipped.
        If ``descending`` is |False| (default), younger nodes will be returned
        before older ones; if |True|, older nodes will be returned before
        younger ones.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a |Node| object as an argument
            and returns |True| if the |Node| object is to be yielded by
            the iterator, or |False| if not. If ``filter_fn`` is |None|
            (defau
        include_leaves : boolean, optional
            If |True| (default), then leaf nodes are included in the iteration.
            If |False|, then leaf nodes are skipped.lt), then all nodes visited will be yielded.
        descending : boolean, optional
            If |False| (default), then younger nodes are visited before older
            ones. If |True|, then older nodes are visited before younger ones.

        Returns
        -------
        :py:class:`collections.Iterator` [|Node|]
            Iterator over age-ordered sequence of nodes in subtree rooted at
            this node.
        """
        # if not descending:
        #     leaves = [nd for nd in self.leaf_iter()]
        #     queued_pairs = []
        #     in_queue = set()
        #     for leaf in leaves:
        #         age_nd_tuple = (leaf.age, leaf)
        #         queued_pairs.insert(bisect.bisect(queued_pairs, age_nd_tuple), age_nd_tuple)
        #         in_queue.add(leaf)
        #     while queued_pairs:
        #         next_el = queued_pairs.pop(0)
        #         age, nd = next_el
        #         in_queue.remove(nd)
        #         p = nd._parent_node
        #         if p and p not in in_queue:
        #             age_nd_tuple = (p.age, p)
        #             queued_pairs.insert(bisect.bisect(queued_pairs, age_nd_tuple), age_nd_tuple)
        #             in_queue.add(p)
        #         if include_leaves or nd.is_internal():
        #             yield nd
        # else:
        #     nds = [(nd.age, nd) for nd in self.preorder_iter()]
        #     nds.sort(reverse=True)
        #     for nd in nds:
        #         if include_leaves or nd[1].is_internal():
        #             yield nd[1]
        nds = [nd for nd in self.preorder_iter()]
        if descending:
            reverse = True
        else:
            reverse = False
        nds.sort(key=lambda x: x.age, reverse=reverse)
        for nd in nds:
            if (include_leaves or nd._child_nodes) and (
                filter_fn is None or filter_fn(nd)
            ):
                yield nd

    def age_order_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Deprecated: use :meth:`Node.ageorder_iter()` instead.
        """
        deprecate.dendropy_deprecation_warning(
            message=(
                "Deprecated since DendroPy 4: 'age_order_iter()' will no longer be"
                " supported in future releases; use 'ageorder_iter()' instead"
            ),
            stacklevel=3,
        )
        return self.ageorder_iter(
            include_leaves=include_leaves, filter_fn=filter_fn, descending=descending
        )

    def apply(self, before_fn=None, after_fn=None, leaf_fn=None):
        r"""
        Applies function ``before_fn`` and ``after_fn`` to all internal nodes and
        ``leaf_fn`` to all terminal nodes in subtree starting with ``self``, with
        nodes visited in pre-order.

        Given a tree with preorder sequence of nodes of
        [a,b,i,e,j,k,c,g,l,m,f,n,h,o,p,]::

                           a
                          / \
                         /   \
                        /     \
                       /       \
                      /         \
                     /           \
                    /             c
                   b             / \
                  / \           /   \
                 /   e         /     f
                /   / \       /     / \
               /   /   \     g     /   h
              /   /     \   / \   /   / \
             i   j       k l   m n   o   p


        the following order of function calls results:

            before_fn(a)
            before_fn(b)
            leaf_fn(i)
            before_fn(e)
            leaf_fn(j)
            leaf_fn(k)
            after_fn(e)
            after_fn(b)
            before_fn(c)
            before_fn(g)
            leaf_fn(l)
            leaf_fn(m)
            after_fn(g)
            before_fn(f)
            leaf_fn(n)
            before_fn(h)
            leaf_fn(o)
            leaf_fn(p)
            after_fn(h)
            after_fn(f)
            after_fn(c)
            after_fn(a)

        Parameters
        ----------
        before_fn : function object or |None|
            A function object that takes a |Node| as its argument.
        after_fn : function object or |None|
            A function object that takes a |Node| as its argument.
        leaf_fn : function object or |None|
            A function object that takes a |Node| as its argument.

        Notes
        -----
        Adapted from work by Mark T. Holder (the ``peyotl`` module of the Open
        Tree of Life Project):

            https://github.com/OpenTreeOfLife/peyotl.git

        """
        stack = [self]
        while stack:
            node = stack.pop()
            if not node._child_nodes:
                if leaf_fn:
                    leaf_fn(node)
                # (while node is the last child of parent ...)
                while (node._parent_node is None) or (
                    node._parent_node._child_nodes[-1] is node
                ):
                    node = node._parent_node
                    if node is not None:
                        if after_fn is not None:
                            after_fn(node)
                    else:
                        break
            else:
                if before_fn is not None:
                    before_fn(node)
                stack.extend([i for i in reversed(node._child_nodes)])
        return

    def set_child_nodes(self, child_nodes):
        """
        Assigns the set of child nodes for this node.

        Results in the ``parent_node`` attribute of each |Node| in ``nodes``
        as well as the ``tail_node`` attribute of corresponding |Edge|
        objects being assigned to ``self``.

        Parameters
        ----------
        child_nodes : collections.Iterable[|Node|]
            The (iterable) collection of child nodes to be assigned this node
            as a parent.
        """
        self.clear_child_nodes()
        # Go through add to ensure book-keeping
        # (e.g. avoiding multiple adds) takes
        # place.
        for nd in child_nodes:
            self.add_child(nd)

    def set_children(self, child_nodes):
        """Deprecated: use :meth:`Node.set_child_nodes()` instead."""
        return self.set_child_nodes(child_nodes)

    def add_child(self, node):
        """
        Adds a child node to this node if it is not already a child.

        Results in the ``parent_node`` attribute of ``node`` as well as the
        ``tail_node`` attribute of ``node.edge`` being assigned to ``self``.

        Parameters
        ----------
        node : |Node|
            The node to be added as a child of this node.

        Returns
        -------
        |Node|
            The node that was added.
        """
        assert node is not self, "Cannot add node as child of itself"
        assert self._parent_node is not node, (
            "Cannot add a node's parent as its child: remove the node from its parent's"
            " child set first"
        )
        node._parent_node = self
        if node not in self._child_nodes:
            self._child_nodes.append(node)
        return node

    def insert_child(self, index, node):
        """
        Adds a child node to this node.

        If the node is already a child of this node, then it is moved
        to the specified position.
        Results in the ``parent_node`` attribute of ``node`` as well as the
        ``tail_node`` attribute of ``node.edge`` being assigned to ``self``.

        Parameters
        ----------
        index : integer
            The index before which to insert the new node.
        node : |Node|
            The node to be added as a child of this node.

        Returns
        -------
        |Node|
            The node that was added.
        """
        node._parent_node = self
        try:
            cur_index = self._child_nodes.index(node)
        except ValueError:
            pass
        else:
            if cur_index == index:
                return
            self._child_nodes.remove(node)
        self._child_nodes.insert(index, node)
        return node

    def new_child(self, **kwargs):
        """
        Create and add a new child to this node.

        Parameters
        ----------
        kwargs : keyword arguments
            Keyword arguments will be passed directly to the |Node|
            constructor (:meth:`Node.__init()__`).

        Returns
        -------
        |Node|
            The new child node that was created and added.
        """
        node = self.__class__(**kwargs)
        return self.add_child(node=node)

    def insert_new_child(self, index, **kwargs):
        """
        Create and add a new child to this node at a particular position.

        Results in the ``parent_node`` attribute of ``node`` as well as the
        ``tail_node`` attribute of ``node.edge`` being assigned to ``self``.

        Parameters
        ----------
        index : integer
            The index before which to insert the new node.
        kwargs : keyword arguments, optional
            Keyword arguments will be passed directly to the |Node|
            constructor (:meth:`Node.__init()__`).

        Returns
        -------
        |Node|
            The new child node that was created and added.
        """
        node = self.__class__(**kwargs)
        return self.insert_child(index=index, node=node)

    def remove_child(self, node, suppress_unifurcations=False):
        """
        Removes a node from the child set of this node.

        Results in the parent of the node being removed set to |None|.  If
        ``suppress_unifurcations`` is |True|, if this node ends up having only one
        child after removal of the specified node, then this node will be
        removed from the tree, with its single child added to the child node
        set of its parent and the edge length adjusted accordingly.
        ``suppress_unifurcations`` should only be |True| for unrooted trees.

        Parameters
        ----------
        node : |Node|
            The node to be removed.
        suppress_unifurcations : boolean, optional
            If |False| (default), no action is taken. If |True|, then if the
            node removal results in a node with degree of two (i.e., a single
            parent and a single child), then it will be removed from
            the tree and its (sole) child will be added as a child of its
            parent (with edge lengths adjusted accordingly).

        Returns
        -------
        |Node|
            The node removed.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        children = self._child_nodes
        if node in children:
            node._parent_node = None
            node.edge.tail_node = None
            index = children.index(node)
            children.remove(node)
            if suppress_unifurcations:
                if self._parent_node:
                    if len(children) == 1:
                        child = children[0]
                        pos = self._parent_node._child_nodes.index(self)
                        self._parent_node.insert_child(pos, child)
                        self._parent_node.remove_child(
                            self, suppress_unifurcations=False
                        )
                        try:
                            child.edge.length += self.edge.length
                        except:
                            pass
                        self._child_nodes = []
                else:
                    to_remove = None
                    if len(children) == 2:
                        if children[0].is_internal():
                            to_remove = children[0]
                            other = children[1]
                        elif children[1].is_internal():
                            to_remove = children[1]
                            other = children[0]
                    if to_remove is not None:
                        try:
                            other.edge.length += to_remove.edge.length
                        except:
                            pass
                        pos = self._child_nodes.index(to_remove)
                        self.remove_child(to_remove, suppress_unifurcations=False)
                        tr_children = to_remove._child_nodes
                        tr_children.reverse()
                        for c in tr_children:
                            self.insert_child(pos, c)
                        to_remove._child_nodes = []
        else:
            raise ValueError("Tried to remove a node that is not listed as a child")
        return node

    def clear_child_nodes(self):
        """
        Removes all child nodes.
        """
        self._child_nodes.clear()

    def reversible_remove_child(self, node, suppress_unifurcations=False):
        """
        This function is a (less-efficient) version of remove_child that also
        returns the data needed by reinsert_nodes to "undo" the removal.

        Returns a list of tuples.  The first element of each tuple is the
        node removed, the other elements are the information needed by
        ``reinsert_nodes`` in order to restore the tree to the same topology as
        it was before the call to ``remove_child.`` If ``suppress_unifurcations`` is False
        then the returned list will contain only one item.

        ``suppress_unifurcations`` should only be called on unrooted trees.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        children = self._child_nodes
        try:
            pos = children.index(node)
        except:
            raise ValueError("Tried to remove a node that is not listed as a child")
        removed = [(node, self, pos, [], None)]
        node._parent_node = None
        node.edge.tail_node = None
        children.remove(node)
        if suppress_unifurcations:
            p = self._parent_node
            if p:
                if len(children) == 1:
                    child = children[0]
                    pos = p._child_nodes.index(self)
                    p.insert_child(pos, child)
                    self._child_nodes = []
                    p.remove_child(self, suppress_unifurcations=False)
                    e = child.edge
                    try:
                        e.length += self.edge.length
                    except:
                        e = None
                    t = (self, p, pos, [child], e)
                    removed.append(t)
            else:
                to_remove = None
                if len(children) == 2:
                    if children[0].is_internal():
                        to_remove = children[0]
                        other = children[1]
                    elif children[1].is_internal():
                        to_remove = children[1]
                        other = children[0]
                if to_remove is not None:
                    e = other.edge
                    try:
                        e.length += to_remove.edge.length
                    except:
                        e = None
                    pos = self._child_nodes.index(to_remove)
                    self.remove_child(to_remove, suppress_unifurcations=False)
                    tr_children = to_remove._child_nodes
                    to_remove._child_nodes = []
                    for n, c in enumerate(tr_children):
                        new_pos = pos + n
                        self.insert_child(pos, c)
                    t = (to_remove, self, pos, tr_children, e)
                    removed.append(t)

        return removed

    def reinsert_nodes(self, nd_connection_list):
        """
        This function should be used to "undo" the effects of
        Node.reversible_remove_child NOTE: the behavior is only
        guaranteed if the tree has not been modified between the
        remove_child and reinsert_nodes calls! (or the tree has been
        restored such that the node/edge identities are identical to the
        state before the remove_child call.

        The order of info in each tuple is:

            0 - node removed
            1 - parent of node removed
            2 - pos in parent matrix
            3 - children of node removed that were "stolen"
            4 - edge that was lengthened by "stealing" length from node's edge
        """
        # we unroll the stack of operations
        for blob in nd_connection_list[-1::-1]:
            # _LOG.debug(blob)
            n, p, pos, children, e = blob
            for c in children:
                cp = c._parent_node
                if cp:
                    cp.remove_child(c)
                n.add_child(c)
            p.insert_child(pos, n)
            if e is not None:
                e.length -= n.edge.length

    def collapse_neighborhood(self, dist):
        if dist < 1:
            return
        children = self.child_nodes()
        for ch in children:
            if not ch.is_leaf():
                ch.edge.collapse()
        if self._parent_node:
            p = self._parent_node
            self.edge.collapse()
            p.collapse_neighborhood(dist - 1)
        else:
            self.collapse_neighborhood(dist - 1)

    def collapse_clade(self):
        """Collapses all internal edges that are descendants of self."""
        if self.is_leaf():
            return
        leaves = [i for i in self.leaf_iter()]
        self.set_child_nodes(leaves)

    def collapse_conflicting(self, bipartition):
        """
        Collapses every edge in the subtree that conflicts with the given
        bipartition. This can include the edge subtending subtree_root.
        """
        to_collapse_head_nodes = []
        for nd in self.postorder_iter():
            if nd._child_nodes and nd.edge.bipartition.is_incompatible_with(
                bipartition
            ):
                to_collapse_head_nodes.append(nd)
        for nd in to_collapse_head_nodes:
            e = nd.edge
            e.collapse()

    def _get_edge(self):
        """
        Returns the edge subtending this node.
        """
        return self._edge

    def _set_edge(self, new_edge):
        """
        Sets the edge subtending this node, and sets head_node of
        ``edge`` to point to self.
        """
        # if edge is None:
        #     raise ValueError("A Node cannot have 'None' for an edge")
        if new_edge is self._edge:
            return
        if self._parent_node is not None:
            try:
                self._parent_node._child_nodes.remove(self)
            except ValueError:
                pass

        ## Minimal management
        self._edge = new_edge
        if self._edge:
            self._edge._head_node = self

    edge = property(_get_edge, _set_edge)

    def _get_edge_length(self):
        """
        Returns the length of the edge subtending this node.
        """
        return self._edge.length

    def _set_edge_length(self, v=None):
        """
        Sets the edge subtending this node, and sets head_node of
        ``edge`` to point to self.
        """
        self._edge.length = v

    edge_length = property(_get_edge_length, _set_edge_length)

    def _get_bipartition(self):
        """
        Returns the bipartition for the edge subtending this node.
        """
        return self._edge.bipartition

    def _set_bipartition(self, v=None):
        """
        Sets the bipartition for the edge subtending this node.
        """
        self._edge.bipartition = v

    bipartition = property(_get_bipartition, _set_bipartition)

    def _get_split_bitmask(self):
        return self._edge.bipartition._split_bitmask

    def _set_split_bitmask(self, h):
        self._edge.bipartition._split_bitmask = h

    split_bitmask = property(_get_split_bitmask, _set_split_bitmask)

    def _get_leafset_bitmask(self):
        return self._edge.bipartition._leafset_bitmask

    def _set_leafset_bitmask(self, h):
        self._edge.bipartition._leafset_bitmask = h

    leafset_bitmask = property(_get_leafset_bitmask, _set_leafset_bitmask)

    def _get_tree_leafset_bitmask(self):
        return self._edge.bipartition._tree_leafset_bitmask

    def _set_tree_leafset_bitmask(self, h):
        self._edge.bipartition._tree_leafset_bitmask = h

    tree_leafset_bitmask = property(
        _get_tree_leafset_bitmask, _set_tree_leafset_bitmask
    )

    def split_as_bitstring(self):
        return self._edge.bipartition.split_as_bitstring()

    def leafset_as_bitstring(self):
        return self._edge.bipartition.leafset_as_bitstring()

    def _get_parent_node(self):
        """Returns the parent node of this node."""
        return self._parent_node

    def _set_parent_node(self, parent):
        """Sets the parent node of this node."""
        if self._parent_node is not None:
            try:
                self._parent_node._child_nodes.remove(self)
            except ValueError:
                pass
        self._parent_node = parent
        if self._parent_node is not None:
            if self not in self._parent_node._child_nodes:
                self._parent_node._child_nodes.append(self)

    parent_node = property(_get_parent_node, _set_parent_node)

    def is_leaf(self):
        """
        Returns |True| if the node is a tip or a leaf node, i.e. has no child
        nodes.

        Returns
        -------
        boolean
            |True| if the node is a leaf, i.e., has no child nodes. |False|
            otherwise.
        """
        return bool(not self._child_nodes)

    def is_internal(self):
        """
        Returns |True| if the node is *not* a tip or a leaf node.

        Returns
        -------
        boolean
            |True| if the node is not a leaf. |False| otherwise.
        """
        return bool(self._child_nodes)

    def leaf_nodes(self):
        """
        Returns list of all leaf_nodes descended from this node (or just
        list with ``self`` as the only member if ``self`` is a leaf).

        Note
        ----
        Usage of  `leaf_iter()` is preferable for efficiency reasons unless
        actual list is required.

        Returns
        -------
        :py:class:`list` [|Node|]
           A ``list`` of |Node| objects descended from this node
           (inclusive of ``self``) that are the leaves.
        """
        return [
            node
            for node in self.postorder_iter(lambda x: bool(len(x.child_nodes()) == 0))
        ]

    def num_child_nodes(self):
        """
        Returns number of child nodes.

        Returns
        -------
        int
            Number of children in ``self``.
        """
        return len(self._child_nodes)

    def child_nodes(self):
        """
        Returns a shallow-copy list of all child nodes of this node.

        Note
        ----
        Unless an actual ``list`` is needed, iterating over the child nodes using
        :meth:`Node.child_node_iter()` is preferable to avoid the overhead of
        list construction.

        Returns
        -------
        :py:class:`list` [|Node|]
           A ``list`` of |Node| objects that have ``self`` as a parent.
        """
        return list(self._child_nodes)

    def child_edges(self):
        """
        Returns a shallow-copy list of all child edges of this node.

        Note
        ----
        Unless an actual ``list`` is needed, iterating over the child edges using
        :meth:`Node.child_edge_iter()` is preferable to avoid the overhead of
        list construction.

        Returns
        -------
        :py:class:`list` [|Edge|]
           A ``list`` of |Edge| objects that have ``self`` as a tail node.
        """
        return list(ch.edge for ch in self._child_nodes)

    def incident_edges(self):
        """
        Return parent and child edges.

        Returns
        -------
        :py:class:`list` [|Edge|]
            A list of edges linking to this node, with outgoing edges (edges
            connecting to child nodes) followed by the edge connecting
            this node to its parent.
        """
        e = [c.edge for c in self._child_nodes]
        e.append(self.edge)
        return e

    def get_incident_edges(self):
        """Legacy synonym for :meth:`Node.incident_edges()`."""
        return self.incident_edges()

    def adjacent_nodes(self):
        """
        Return parent and child nodes.

        Returns
        -------
        :py:class:`list` [|Node|]
            A list with all child nodes and parent node of this node.
        """
        n = [c for c in self._child_nodes]
        if self._parent_node:
            n.append(self._parent_node)
        return n

    def get_adjacent_nodes(self):
        """Legacy synonym for :meth:`Node.adjacent_edges()`"""
        return self.adjacent_nodes()

    def sibling_nodes(self):
        """
        Return all other children of parent, excluding self.

        Returns
        -------
        :py:class:`list` [|Node|]
            A list of all nodes descended from the same parent as ``self``,
            excluding ``self``.
        """
        p = self._parent_node
        if not p:
            return []
        sisters = [nd for nd in p.child_nodes() if nd is not self]
        return sisters

    def sister_nodes(self):
        """Legacy synonym for :meth:`Node.sister_nodes()`"""
        return self.sibling_nodes()

    def extract_subtree(
        self,
        extraction_source_reference_attr_name="extraction_source",
        node_filter_fn=None,
        suppress_unifurcations=True,
        is_apply_filter_to_leaf_nodes=True,
        is_apply_filter_to_internal_nodes=False,
        node_factory=None,
    ):
        """
        Returns a clone of the structure descending from this node.

        Parameters
        ----------
        extraction_source_reference_attr_name : str
            Name of attribute to set on cloned nodes that references
            corresponding original node. If ``None``, then attribute (and
            reference) will not be created.
        node_filter_fn : None or function object
            If ``None``, then entire tree structure is cloned.
            If not ``None``, must be a function object that returns ``True``
            if a particular |Node| instance on the original tree should
            be included in the cloned tree, or ``False`` otherwise.
        is_apply_filter_to_leaf_nodes : bool
            If ``True`` then the above filter will be applied to leaf nodes. If
            ``False`` then it will not (and all leaf nodes will be
            automatically included, unless excluded by an ancestral node being
            filtered out).
        is_apply_filter_to_internal_nodes : bool
            If ``True`` then the above filter will be applied to internal nodes. If
            ``False`` then it will not (internal nodes without children will
            still be filtered out).
        node_factory : function
            If not ``None``, must be a function that takes no arguments and
            returns a new |Node| (or equivalent) instance.

        Returns
        -------
        nd : |Node|
            A node with descending subtree mirroring this one.

        """
        memo = {}
        is_excluded_nodes = False
        start_node = None
        start_node_to_match = self
        if node_factory is None:
            node_factory = self.__class__
        nd1 = None  # verbosity to mollify linter
        for nd0 in self.postorder_iter():
            if node_filter_fn is not None:
                if nd0._child_nodes:
                    if is_apply_filter_to_internal_nodes:
                        is_apply_filter = True
                    else:
                        is_apply_filter = False
                else:
                    if is_apply_filter_to_leaf_nodes:
                        is_apply_filter = True
                    else:
                        is_apply_filter = False
                if is_apply_filter and not node_filter_fn(nd0):
                    is_excluded_nodes = True
                    continue
            original_node_has_children = False
            children_to_add = []
            for ch_nd0 in nd0.child_node_iter():
                original_node_has_children = True
                ch_nd1 = memo.get(ch_nd0, None)
                if ch_nd1 is not None:
                    children_to_add.append(ch_nd1)
            if not children_to_add and original_node_has_children:
                # filter removes all descendents of internal node,
                # so this internal node is not added
                if nd0.parent_node is None:
                    raise error.SeedNodeDeletionException(
                        "Attempting to remove seed node or node without parent"
                    )
                if nd0 is self:
                    start_node_to_match = nd0.parent_node
                continue
            elif len(children_to_add) == 1 and suppress_unifurcations:
                if nd0.edge.length is not None:
                    if children_to_add[0].edge.length is None:
                        children_to_add[0].edge.length = nd0.edge.length
                    else:
                        children_to_add[0].edge.length += nd0.edge.length
                else:
                    nd1.edge.length = children_to_add[0].edge.length
                if nd0.parent_node is None:
                    start_node = children_to_add[0]
                    break
                if nd0 is self:
                    start_node_to_match = nd0.parent_node
                memo[nd0] = children_to_add[0]
            else:
                nd1 = node_factory()
                nd1.label = nd0.label
                nd1.taxon = nd0.taxon
                nd1.edge.length = nd0.edge.length
                nd1.edge.label = nd0.edge.label
                for ch_nd1 in children_to_add:
                    nd1.add_child(ch_nd1)
                if nd0 is start_node_to_match:
                    start_node = nd1
                memo[nd0] = nd1
                if extraction_source_reference_attr_name:
                    setattr(nd1, extraction_source_reference_attr_name, nd0)
        if start_node is not None:
            return start_node
        else:
            ## TODO: find a replacement node
            raise ValueError

    def level(self):
        """
        Returns the number of nodes between ``self`` and the seed node of the tree.

        Returns
        -------
        integer
            The number of nodes between ``self`` and the seed node of the tree,
            or 0 if ``self`` has no parent.
        """
        if self._parent_node:
            return self._parent_node.level() + 1
        else:
            return 0

    def distance_from_root(self):
        """
        Weighted path length of ``self`` from root.

        Returns
        -------
        numeric
            Total weight of all edges connecting ``self`` with the root of the
            tree.
        """
        if self._parent_node and self.edge.length != None:
            if self._parent_node.distance_from_root == None:
                return float(self.edge.length)
            else:
                distance_from_root = float(self.edge.length)
                parent_node = self._parent_node
                # The root is identified when a node with no
                # parent is encountered. If we want to use some
                # other criteria (e.g., where a is_root property
                # is True), we modify it here.
                while parent_node:
                    if parent_node.edge.length != None:
                        distance_from_root = distance_from_root + float(
                            parent_node.edge.length
                        )
                    parent_node = parent_node._parent_node
                return distance_from_root
        elif not self._parent_node and self.edge.length != None:
            return float(self.edge.length)
        elif self._parent_node and self.edge.length == None:
            # what do we do here: parent node exists, but my
            # length does not?
            return float(self._parent_node.edge.length)
        elif not self._parent_node and self.edge.length == None:
            # no parent node, and no edge length
            return 0.0
        else:
            # WTF????
            return 0.0

    def distance_from_tip(self):
        """
        Maximum weighted length of path of ``self`` to tip.

        If tree is not ultrametric (i.e., descendent edges have different
        lengths), then count the maximum of edge lengths. Note that
        :meth:`Tree.calc_node_ages()` is a more efficient way of doing this
        over the whole tree if this value is need for many or all the nodes on
        the tree.

        Returns
        -------
        numeric
            Maximum weight of edges connecting ``self`` to tip.
        """
        if not self._child_nodes:
            return 0.0
        else:
            distance_from_tips = []
            for ch in self._child_nodes:
                if ch.edge.length is not None:
                    curr_edge_length = ch.edge_length
                else:
                    curr_edge_length = 0.0
                if not hasattr(ch, "_distance_from_tip"):
                    ch._distance_from_tip = ch.distance_from_tip()
                distance_from_tips.append(ch._distance_from_tip + curr_edge_length)
            self._distance_from_tip = float(max(distance_from_tips))
            return self._distance_from_tip

    def description(
        self, depth=1, indent=0, itemize="", output=None, taxon_namespace=None
    ):
        """
        Returns description of object, up to level ``depth``.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        label = str(self)
        output_strio.write(
            "%s%sNode object at %s%s" % (indent * " ", itemize, hex(id(self)), label)
        )
        if depth >= 1:
            leader1 = " " * (indent + 4)
            leader2 = " " * (indent + 8)
            output_strio.write("\n%s[Edge]" % leader1)
            if self.edge is not None:
                edge_desc = self.edge.description(0)
            else:
                edge_desc = "None"
            output_strio.write("\n%s%s" % (leader2, edge_desc))

            output_strio.write("\n%s[Taxon]" % leader1)
            if self.taxon is not None:
                taxon_desc = self.taxon.description(0)
            else:
                taxon_desc = "None"
            output_strio.write("\n%s%s" % (leader2, taxon_desc))

            output_strio.write("\n%s[Parent]" % leader1)
            if self._parent_node is not None:
                parent_node_desc = self._parent_node.description(0)
            else:
                parent_node_desc = "None"
            output_strio.write("\n%s%s" % (leader2, parent_node_desc))
            output_strio.write("\n%s[Children]" % leader1)
            if len(self._child_nodes) == 0:
                output_strio.write("\n%sNone" % leader2)
            else:
                for i, cnd in enumerate(self._child_nodes):
                    output_strio.write("\n%s[%d] %s" % (leader2, i, cnd.description(0)))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    ## For debugging we build-in a full-fledged NEWICK composition independent
    ## of the nexus/newick family of modules. Client code should prefer to
    ## use Newick/Nexus readers/writers, or Tree.write(), TreeList.write(),
    ## DataSet.write() etc.

    def _as_newick_string(self, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.
        For production purposes, use the the full-fledged 'as_string()'
        method of the object.
        """
        out = StringIO()
        self._write_newick(out, **kwargs)
        return out.getvalue()

    def _write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.  For
        production purposes, use the the full-fledged 'write_to_stream()'
        method of the object.
        """
        edge_lengths = not kwargs.get("suppress_edge_lengths", False)
        edge_lengths = kwargs.get("edge_lengths", edge_lengths)
        child_nodes = self.child_nodes()
        if child_nodes:
            out.write("(")
            f_child = child_nodes[0]
            for child in child_nodes:
                if child is not f_child:
                    out.write(",")
                child._write_newick(out, **kwargs)
            out.write(")")

        out.write(self._get_node_token(**kwargs))
        if edge_lengths:
            e = self.edge
            if e:
                sel = e.length
                if sel is not None:
                    fmt = kwargs.get("edge_length_formatter", None)
                    if fmt:
                        out.write(":%s" % fmt(sel))
                    else:
                        s = ""
                        try:
                            s = float(sel)
                            s = str(s)
                        except ValueError:
                            s = str(sel)
                        if s:
                            out.write(":%s" % s)

    def _get_node_token(self, **kwargs):
        """returns a string that is an identifier for the node.  This is called
        by the newick-writing functions, so the kwargs that affect how node
        labels show up in a newick string are the same ones used here:
        ``suppress_internal_labels`` is a Boolean, and defaults to False.
        """
        is_leaf = len(self._child_nodes) == 0
        if not is_leaf:
            if kwargs.get("suppress_internal_labels", False) or not kwargs.get(
                "include_internal_labels", True
            ):
                return ""
        if self.taxon is not None:
            if self.taxon.label:
                label = self.taxon.label
            else:
                # return "_" # taxon, but no label: anonymous
                label = (  # "_" is not anonoymous/unnamed, but a name of <blank>; so we return nothing instead
                    ""
                )
        else:
            if self.label:
                label = self.label
            else:
                label = ""
        if not label or kwargs.get("raw_labels", False):
            return label
        elif " " in label and "_" in label:
            if "'" in label:
                label.replace("'", "''")
            return "'{}'".format(label)
        elif " " in label and not kwargs.get("preserve_spaces"):
            return label.replace(" ", "_")
        else:
            return label

    def _get_indented_form(self, **kwargs):
        out = StringIO()
        self._write_indented_form(out, **kwargs)
        return out.getvalue()

    def _write_indented_form(self, out, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        level = kwargs.get("level", 0)
        ancestors = []
        siblings = []
        n = self
        while n is not None:
            n._write_indented_form_line(out, level, **kwargs)
            n, lev = n._preorder_list_manip(siblings, ancestors)
            level += lev

    def _get_indented_form_line(self, level, **kwargs):
        out = StringIO()
        self._write_indented_form_line(out, level, **kwargs)
        return out.getvalue()

    def _write_indented_form_line(self, out, level, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        label = self._format_node(**kwargs)
        if kwargs.get("bipartitions"):
            cm = "%s " % self.edge.bipartition._format_bipartition(**kwargs)
        else:
            cm = ""
        out.write("%s%s%s\n" % (cm, indentation * level, label))

    def _format_node(self, **kwargs):
        nf = kwargs.get('node_formatter', None)
        if nf:
            return nf(self)
        if self.taxon is not None:
            return str(self.taxon)
        if self.label is not None:
            return self.label
        return ""

    def _preorder_list_manip(self, siblings, ancestors):
        """
        Helper function for recursion free preorder traversal, that does
        not rely on attributes of the node other than child_nodes() (thus it
        is useful for debuggging).

        Returns the next node (or None) and the number of levels toward the
        root the function "moved".
        """
        levels_moved = 0
        c = self.child_nodes()
        if c:
            levels_moved += 1
            ancestors.append(list(siblings))
            del siblings[:]
            siblings.extend(c[1:])
            return c[0], levels_moved
        while not siblings:
            if ancestors:
                levels_moved -= 1
                del siblings[:]
                siblings.extend(ancestors.pop())
            else:
                return None, levels_moved
        return siblings.pop(0), levels_moved

    def _convert_node_to_root_polytomy(self):
        """If ``self`` has two children and at least one of them is an internal node,
        then it will be converted to an out-degree three node (with the edge length
        added as needed).

        Returns a tuple of child nodes that were detached (or() if the tree was not
        modified). This can be useful for removing the deleted node from the split_edge_map
        dictionary.
        """
        nd_children = self.child_nodes()
        if len(nd_children) > 2:
            return ()
        try:
            left_child = nd_children[0]
        except:
            return ()
        if not left_child:
            return ()
        if len(nd_children) == 1:
            right_child = None
            dest_edge_head = self
        else:
            right_child = nd_children[1]
            dest_edge_head = right_child
        curr_add = None
        if right_child and right_child.is_internal():
            try:
                left_child.edge.length += right_child.edge.length
            except:
                pass
            self.remove_child(right_child)
            grand_kids = right_child.child_nodes()
            for gc in grand_kids:
                self.add_child(gc)
            curr_add = right_child
        elif left_child.is_internal():
            try:
                dest_edge_head.edge.length += left_child.edge.length
            except:
                pass
            self.remove_child(left_child)
            grand_kids = left_child.child_nodes()
            for gc in grand_kids:
                self.add_child(gc)
            curr_add = left_child
        if curr_add:
            ndl = [curr_add]
            t = self._convert_node_to_root_polytomy()
            ndl.extend(t)
            return tuple(ndl)
        return ()
