#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
This module handles the core definition of tree data structure class,
as well as all the structural classes that make up a tree.
"""

import warnings
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.utility import terminal
from dendropy.datamodel import base
from dendropy.datamodel import taxon
from dendropy import dataio

##############################################################################
## Edge

class Edge(base.DataObject, base.Annotable):
    """
    An :term:`edge` on a :term:`tree`.
    """

    def __init__(self,
            tail_node=None,
            head_node=None,
            length=None,
            rootedge=False,
            label=None):
        """
        Parameters
        ----------
        tail_node : :class:`Node`
            Node from which this edge originates, i.e., the parent node of this
            edge's `head_node`.

        head_node : :class:`Node`
            Node from to which this edge links, i.e., the child node of this
            node `tail_node`.

        length : numerical
            A value representing the weight of the edge.

        rootedge : boolean
            Is the child node of this edge the root or seed node of the tree?

        label : string
            Label for this edge.

        """
        base.DataObject.__init__(self, label=label)
        self.tail_node = tail_node
        self.head_node = head_node
        self.rootedge = rootedge
        self.length = length
        self.split_bitmask = None
        self.comments = []

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def __copy__(self):
        """
        All member objects are references, except for :attr:`annotation_set` of
        top-level object and member :class:`Annotation` objects: these are
        full, independent instances (though any complex objects in the `value`
        field of :class:`Annotation` objects are also just references).
        """
        other = self.__class__.__new__(self.__class__)
        memo = {}
        memo[id(self)] = other
        for k in self.__dict__:
            if k == "_annotations":
                continue
            other.__dict__[k] = self.__dict__[k]
            memo[id(k)] = k
            memo[id(self.__dict__[k])] = other.__dict__[k]
        self.deep_copy_annotations_from(other, memo=memo)

    def __deepcopy__(self):
        """
        Exhaustive deep-copy: all objects are cloned.
        """
        raise NotImplementedError

    def is_leaf(self):
        "Returns True if the head node has no children"
        return self.head_node and self.head_node.is_leaf()

    def is_terminal(self):
        return self.is_leaf()

    def is_internal(self):
        "Returns True if the head node has children"
        return self.head_node and not self.head_node.is_leaf()

    def get_adjacent_edges(self):
        """
        Returns a list of all edges that "share" a node with `self`.
        """
        he = [i for i in self.head_node.get_incident_edges() if i is not self]
        te = [i for i in self.tail_node.get_incident_edges() if i is not self]
        he.extend(te)
        return he
    adjacent_edges = property(get_adjacent_edges)

    ###########################################################################
    ## Structural

    def collapse(self):
        """
        Inserts all children of the head_node of self as children of the
        tail_node of self in the same place in the child_node list that
        head_node had occupied. The edge length and head_node will no longer be
        part of the tree.
        """
        to_del = self.head_node
        p = self.tail_node
        if not p:
            return
        children = to_del.child_nodes()
        if not children:
            raise ValueError('collapse_self called with a terminal.')
        pos = p.child_nodes().index(to_del)
        p.remove_child(to_del)
        for child in children:
            p.add_child(child, pos=pos)
            pos += 1

    def invert(self):
        """
        Changes polarity of edge.
        """
        self.head_node, self.tail_node = self.tail_node, self.head_node

    ###########################################################################
    ## Representation

    def description(self,
            depth=1,
            indent=0,
            itemize="",
            output=None,
            taxon_namespace=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s, Length=%s)" % (self.oid, str(self.length))
        else:
            label = " (%s: '%s', Length=%s)" % (self.oid, self.label, str(self.length))
        output_strio.write('%s%sEdge object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            leader1 = ' ' * (indent + 4)
            leader2 = ' ' * (indent + 8)
            output_strio.write('\n%s[Length]' % leader1)
            if self.length is not None:
                length = self.length
            else:
                length = "None"
            output_strio.write('\n%s%s' % (leader2, length))
            output_strio.write('\n%s[Tail Node]' % leader1)
            if self.tail_node is not None:
                tn = self.tail_node.description(0)
            else:
                tn = "None"
            output_strio.write('\n%s%s' % (leader2, tn))
            output_strio.write('\n%s[Head Node]' % leader1)
            if self.head_node is not None:
                hn = self.head_node.description(0)
            else:
                hn = "None"
            output_strio.write('\n%s%s' % (leader2, hn))
            if hasattr(self, 'split_bitmask'):
                output_strio.write('\n%s[Clade Mask]' % leader1)
                if taxon_namespace is None:
                    output_strio.write('\n%s%s' % (leader2, self.split_bitmask))
                else:
                    output_strio.write('\n%s%s' % (leader2, taxon_namespace.split_bitmask_string(self.split_bitmask)))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

##############################################################################
## Node

class Node(base.DataObject, base.Annotable):
    """
    A :term:`node` on a :term:`tree`.
    """

    ###########################################################################
    ## Life-cycle

    def __init__(self,
            taxon=None,
            label=None,
            edge_length=None):
        """
        Parameters
        ----------

        taxon : :class:`Taxon`
            The :class:`Taxon` instance representing the operational taxonomic
            unit concept associated with this Node.
        label : string
            A label for this node.
        edge_length : numeric
            Length or weight of the edge subtending this node.

        """
        base.DataObject.__init__(self, label=label)
        self.taxon = taxon
        self.age = None
        self._edge = None
        self._child_nodes = []
        self._parent_node = None
        self.edge = Edge(head_node=self, length=edge_length)
        self.comments = []

    # def __deepcopy__(self, memo):
    #     memo[id(self._child_nodes)] = []
    #     o = TaxonLinked.__deepcopy__(self, memo)
    #     memo[id(self._child_nodes)] = o._child_nodes
    #     for c in self._child_nodes:
    #         o.add_child(copy.deepcopy(c, memo))
    #     return o

    ###########################################################################
    ## Identity

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        ### IMPORTANT LESSON LEARNED: if you define __hash__, you *must* define __eq__
        return self is other

    def __repr__(self):
        return "<Node object at {}: '{}' ({})>".format(hex(id(self)), self.label, repr(self.taxon))

    ###########################################################################
    ## Iterators

    def preorder_iter(self, filter_fn=None):
        """
        Pre-order iterator over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited before its
        children. Nodes can optionally be filtered by `filter_fn`: only nodes
        for which `filter_fn` returns `True` when called with the node as an
        argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding nodes of the subtree rooted at this node in
            pre-order sequence.
        """
        stack = [self]
        while stack:
            node = stack.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = list(node._child_nodes)
            child_nodes.extend(stack)
            stack = child_nodes

    def preorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes of subtree rooted at this node.

        Visits self and all internal descendant nodes, with each node visited
        before its children. In DendroPy, "internal nodes" are nodes that have
        at least one child node, and thus the root or seed node is typically included
        unless `exclude_seed_node` is `True`. Nodes can optionally be filtered
        by `filter_fn`: only nodes for which `filter_fn` returns `True` when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If `False` (default), then the seed node or root is visited. If
            `True`, then the seed node is skipped.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
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

        Visits self and all descendant nodes, with each node visited first
        followed by its children. Nodes can optionally be filtered by
        `filter_fn`: only nodes for which `filter_fn` returns `True` when
        called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
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
        stack = [(self, False)]
        while stack:
            node, state = stack.pop(0)
            if state:
                if filter_fn is None or filter_fn(node):
                    yield node
            else:
                stack.insert(0, (node, True))
                child_nodes = [(n, False) for n in node._child_nodes]
                child_nodes.extend(stack)
                stack = child_nodes

    def postorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes of subtree rooted at this node.

        Visits self and all internal descendant nodes, with each node visited
        after its children. In DendroPy, "internal nodes" are nodes that have
        at least one child node, and thus the root or seed node is typically
        included unless `exclude_seed_node` is `True`. Nodes can optionally be
        filtered by `filter_fn`: only nodes for which `filter_fn` returns
        `True` when passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If `False` (default), then the seed node or root is visited. If
            `True`, then the seed node is skipped.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
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
        Nodes can optionally be filtered by `filter_fn`: only nodes for which
        `filter_fn` returns `True` when called with the node as an argument are
        visited.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
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
        warnings.warn("Use 'levelorder_iter()' instead of 'level_order_iter()'",
                FutureWarning, stacklevel=2)
        return self.levelorder_iter(filter_fn=filter_fn)

    def inorder_iter(self, filter_fn=None):
        """
        In-order iteration over nodes of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited in-between
        its children. Only valid for strictly-bifurcating trees. Nodes can
        optionally be filtered by `filter_fn`: only nodes for which `filter_fn`
        returns `True` when called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
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
        optionally be filtered by `filter_fn`: only nodes for which `filter_fn`
        returns `True` when called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
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
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
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
            A function object that takes a :class:`Edge` object as an argument
            and returns `True` if the :class:`Edge` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all edges visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Edge`]
            An iterator yielding edges that have this edge as a parent.
        """
        for node in self._child_nodes:
            if filter_fn is None or filter_fn(node.edge):
                yield node.edge

    def ancestor_iter(self, filter_fn=None, inclusive=False):
        """
        Iterator over all ancestors of this node.

        Visits all nodes that are the ancestors of this node.  If `inclusive`
        is `True`, `self` is returned as the first item of the sequence;
        otherwise `self` is skipped. Nodes can optionally be filtered by
        `filter_fn`: only nodes for which `filter_fn` returns `True` when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.
        inclusive : boolean, optional
            If `True`, includes this node in the sequence. If `False`, this is
            skipped.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            Iterator over all predecessor/ancestor nodes of this node.
        """
        if inclusive and (filter_fn is None or filter_fn(self)):
            yield self
        node = self
        while node is not None:
            node = node.parent_node
            if node is not None \
                   and (filter_fn is None or filter_fn(node)):
                yield node

    def ageorder_iter(self, filter_fn=None, include_leaves=True, descending=False):
        """
        Iterator over nodes of subtree rooted at this node in order of the age
        of the node (i.e., the time since the present).

        Iterates over nodes in order of age ('age' is as given by the `age`
        attribute, which is usually the sum of edge lengths from tips
        to node, i.e., time since present).
        If `include_leaves` is `True` (default), leaves are included in the
        iteration; if `include_leaves` is `False`, leaves will be skipped.
        If `descending` is `False` (default), younger nodes will be returned
        before older ones; if `True`, older nodes will be returned before
        younger ones.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (defau
        include_leaves : boolean, optional
            If `True` (default), then leaf nodes are included in the iteration.
            If `False`, then leaf nodes are skipped.lt), then all nodes visited will be yielded.
        descending : boolean, optional
            If `False` (default), then younger nodes are visited before older
            ones. If `True`, then older nodes are visited before younger ones.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
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
        #         p = nd.parent_node
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
            if (include_leaves or nd._child_nodes) and (filter_fn is None or filter_fn(nd)):
                yield nd

    def age_order_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Deprecated: use :meth:`Node.ageorder_iter()` instead.
        """
        warnings.warn("Use 'ageorder_iter()' instead of 'age_order_iter()'",
                FutureWarning, stacklevel=2)
        return self.ageorder_iter(include_leaves=include_leaves,
                filter_fn=filter_fn,
                descending=descending)

    ###########################################################################
    ## Child Node Access and Manipulation

    def set_child_nodes(self, child_nodes):
        """
        Assigns the set of child nodes for this node.

        Results in the `parent_node` attribute of each :class:`Node` in `nodes
        as well as the `tail_node` attribute of corresponding :class:`Edge`
        objects being assigned to `self`.

        Parameters
        ----------
        child_nodes : collections.Iterable[:class:`Node`]
            The (iterable) collection of child nodes to be assigned this node
            as a parent.
        """
        self._child_nodes = list(child_nodes)
        for nd in self._child_nodes:
            nd.parent_node = self
            nd.edge.tail_node = self

    def set_children(self, child_nodes):
        """Deprecated: use :meth:`Node.set_child_nodes()` instead."""
        return self.set_child_nodes(child_nodes)

    def add_child(self, node, pos=None):
        """
        Adds a child node to this node.

        Results in the `parent_node` attribute of `node` as well as the
        `tail_node` attribute of `node.edge` being assigned to `self`.

        Parameters
        ----------
        node : :class:`Node`
            The node to be added as a child of this node.
        pos : integer, optional
            If not `None`, the position in the the sequence of children that
            the new child node should occupy.

        Returns
        -------
        node : :class:`Node`
            The node that was added.
        """
        node.parent_node = self
        ## Support for this was removed, due to unclear expected behavior when
        ## `None` is passed: does the client code expect the edge length to be set
        ## to `None` or simply left untouched? Better approach: client code
        ## explictly sets the edge.
        # if edge_length != None:
        #     node.edge_length = edge_length
        if pos is None:
            self._child_nodes.append(node)
        else:
            self._child_nodes.insert(pos, node)
        return node

    def new_child(self, **kwargs):
        """
        Create and add a new child to this node.

        Parameters
        ----------
        \*\*kwargs : keyword arguments
            Keyword arguments will be passed directly to the :class:`Node`
            constructor (:meth:`Node.__init()__`).

        Returns
        -------
        node : :class:`Node`
            The new child node that was created and added.
        """
        node = self.__class__(**kwargs)
        return self.add_child(node=node)

    def insert_new_child(self, pos, **kwargs):
        """
        Create and add a new child to this node at a particular position.

        Results in the `parent_node` attribute of `node` as well as the
        `tail_node` attribute of `node.edge` being assigned to `self`.

        Parameters
        ----------
        pos : integer
            The position in the the sequence of children of `self` that
            the new child node should occupy.
        \*\*kwargs : keyword arguments, optional
            Keyword arguments will be passed directly to the :class:`Node`
            constructor (:meth:`Node.__init()__`).

        Returns
        -------
        node : :class:`Node`
            The new child node that was created and added.
        """
        node = self.__class__(**kwargs)
        return self.add_child(node=node, pos=pos)

    def remove_child(self, node, suppress_deg_two=False):
        """
        Removes a node from the child set of this node.

        Results in the parent of the node being removed set to `None`.  If
        `suppress_deg_two` is `True`, if this node ends up having only one
        child after removal of the specified node, then this node will be
        removed from the tree, with its single child added to the child node
        set of its parent and the edge length adjusted accordingly.
        `suppress_deg_two` should only be `True` for unrooted trees.

        Parameters
        ----------
        node : :class:`Node`
            The node to be removed.
        suppress_deg_two : boolean, optional
            If `False` (default), no action is taken. If `True`, then if the
            node removal results in a node with degree of two (i.e., a single
            parent and a single child), then it will be removed from
            the tree and its (sole) child will be added as a child of its
            parent (with edge lengths adjusted accordingly).

        Returns
        -------
        node : :class:`Node`
            The node removed.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        children = self._child_nodes
        if node in children:
            node.parent_node = None
            node.edge.tail_node = None
            index = children.index(node)
            children.remove(node)
            if suppress_deg_two:
                if self.parent_node:
                    if len(children) == 1:
                        child = children[0]
                        pos = self.parent_node._child_nodes.index(self)
                        self.parent_node.add_child(child, pos=pos)
                        self.parent_node.remove_child(self, suppress_deg_two=False)
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
                        self.remove_child(to_remove, suppress_deg_two=False)
                        tr_children = to_remove._child_nodes
                        tr_children.reverse()
                        for c in tr_children:
                            self.add_child(c, pos=pos)
                        to_remove._child_nodes = []
        else:
            raise ValueError("Tried to remove a node that is not listed as a child")
        return node

    def reversible_remove_child(self, node, suppress_deg_two=False):
        """
        This function is a (less-efficient) version of remove_child that also
        returns the data needed by reinsert_nodes to "undo" the removal.

        Returns a list of tuples.  The first element of each tuple is the
        node removed, the other elements are the information needed by
        `reinsert_nodes' in order to restore the tree to the same topology as
        it was before the call to `remove_child.` If `suppress_deg_two` is False
        then the returned list will contain only one item.

        `suppress_deg_two` should only be called on unrooted trees.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        children = self._child_nodes
        try:
            pos = children.index(node)
        except:
            raise ValueError("Tried to remove a node that is not listed as a child")
        removed = [(node, self, pos, [], None)]
        node.parent_node = None
        node.edge.tail_node = None
#             if index > 0:
#                 self._child_nodes[index-1].next_sib = None
        children.remove(node)
        if suppress_deg_two:
            p = self.parent_node
            if p:
                if len(children) == 1:
                    child = children[0]
                    pos = p._child_nodes.index(self)
                    p.add_child(child, pos=pos)
                    self._child_nodes = []
                    p.remove_child(self, suppress_deg_two=False)
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
                    self.remove_child(to_remove, suppress_deg_two=False)
                    tr_children = to_remove._child_nodes
                    to_remove._child_nodes = []
                    for n, c in enumerate(tr_children):
                        new_pos = pos + n
                        self.add_child(c, pos=new_pos)
                    t = (to_remove, self, pos, tr_children, e)
                    removed.append(t)

        return removed

    def reinsert_nodes(self, nd_connection_list):
        """This function should be used to "undo" the effects of
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
            #_LOG.debug(blob)
            n, p, pos, children, e = blob
            for c in children:
                cp = c.parent_node
                if cp:
                    cp.remove_child(c)
                n.add_child(c)
            p.add_child(n, pos=pos)
            if e is not None:
                e.length -= n.edge.length

    def collapse_neighborhood(self, dist):
        if dist < 1:
            return
        children = self.child_nodes()
        for ch in children:
            if not ch.is_leaf():
                ch.edge.collapse()
        if self.parent_node:
            p = self.parent_node
            self.edge.collapse()
            p.collapse_neighborhood(dist -1)
        else:
            self.collapse_neighborhood(dist - 1)

    def collapse_clade(self):
        """Collapses all internal edges that are descendants of self."""
        if self.is_leaf():
            return
        leaves = [i for i in self.leaf_iter()]
        self.set_child_nodes(leaves)

    ###########################################################################
    ## Edge Access and Manipulation

    def _get_edge(self):
        """
        Returns the edge subtending this node.
        """
        return self._edge
    def _set_edge(self, edge=None):
        """
        Sets the edge subtending this node, and sets head_node of
        `edge` to point to self.
        """
        self._edge = edge
        if edge:
            edge.head_node = self
    edge = property(_get_edge, _set_edge)

    def _get_edge_length(self):
        """
        Returns the length of the edge subtending this node.
        """
        return self._edge.length
    def _set_edge_length(self, v=None):
        """
        Sets the edge subtending this node, and sets head_node of
        `edge` to point to self.
        """
        self._edge.length = v
    edge_length = property(_get_edge_length, _set_edge_length)

    def _get_split_bitmask(self):
        """
        Returns the split hash bitmask for this node.
        """
        return self._edge.split_bitmask
    def _set_split_bitmask(self, v=None):
        """
        Sets the split hash bitmask for htis node.
        """
        self._edge.split_bitmask = v
    split_bitmask = property(_get_split_bitmask, _set_split_bitmask)

    ###########################################################################
    ## Parent Access and Manipulation

    def _get_parent_node(self):
        """Returns the parent node of this node."""
        return self._parent_node
    def _set_parent_node(self, parent):
        """Sets the parent node of this node."""
        self._parent_node = parent
        self.edge.tail_node = parent
    parent_node = property(_get_parent_node, _set_parent_node)

    ###########################################################################
    ## General Structural Access and Information

    def is_leaf(self):
        """
        Returns `True` if the node is a tip or a leaf node, i.e. has no child
        nodes.

        Returns
        -------
        leafness : boolean
            `True` if the node is a leaf, i.e., has no child nodes. `False`
            otherwise.
        """
        return bool(not self._child_nodes)

    def is_internal(self):
        """
        Returns `True` if the node is *not* a tip or a leaf node.

        Returns
        -------
        internalness : boolean
            `True` if the node is not a leaf. `False` otherwise.
        """
        return bool(self._child_nodes)

    def leaf_nodes(self):
        """
        Returns list of all leaf_nodes descended from this node (or just
        list with `self` as the only member if `self` is a leaf).

        Note
        ----
        Usage of  `leaf_iter()` is preferable for efficiency reasons unless
        actual list is required.

        Returns
        -------
        leaves : :py:class:`list` [:class:`Node`]
           A `list` of :class:`Node` objects descended from this node
           (inclusive of `self`) that are the leaves.
        """
        return [node for node in \
                self.postorder_iter(lambda x: bool(len(x.child_nodes())==0))]

    def num_child_nodes(self):
        """
        Returns number of child nodes.

        Returns
        -------
        n : int
            Number of children in `self`.
        """
        return len(self._child_nodes)

    def child_nodes(self):
        """
        Returns a shallow-copy list of all child nodes of this node.

        Note
        ----
        Unless an actual `list` is needed, iterating over the child nodes using
        :meth:`Node.child_node_iter()` is preferable to avoid the overhead of
        list construction.

        Returns
        -------
        children : :py:class:`list` [:class:`Node`]
           A `list` of :class:`Node` objects that have `self` as a parent.
        """
        return list(self._child_nodes)

    def child_edges(self):
        """
        Returns a shallow-copy list of all child edges of this node.

        Note
        ----
        Unless an actual `list` is needed, iterating over the child edges using
        :meth:`Node.child_edge_iter()` is preferable to avoid the overhead of
        list construction.

        Returns
        -------
        children : :py:class:`list` [:class:`Edge`]
           A `list` of :class:`Edge` objects that have `self` as a tail node.
        """
        return list(ch.edge for ch in self._child_nodes)

    def incident_edges(self):
        """
        Return parent and child edges.

        Returns
        -------
        edges : :py:class:`list` [:class:`Edge`]
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
        nodes : :py:class:`list` [:class:`Node`]
            A list with all child nodes and parent node of this node.
        """
        n = [c for c in self._child_nodes]
        if self.parent_node:
            n.append(self.parent_node)
        return n

    def get_adjacent_nodes(self):
        """Legacy synonym for :meth:`Node.adjacent_edges()`"""
        return self.adjacent_nodes()

    def sibling_nodes(self):
        """
        Return all other children of parent, excluding self.

        Returns
        -------
        siblings : :py:class:`list` [:class:`Node`]
            A list of all nodes descended from the same parent as `self`,
            excluding `self`.
        """
        p = self.parent_node
        if not p:
            return []
        sisters = [nd for nd in p.child_nodes() if nd is not self]
        return sisters

    def sister_nodes(self):
        """Legacy synonym for :meth:`Node.sister_nodes()`"""
        return self.sibling_nodes()

    ###########################################################################
    ## Metrics

    def level(self):
        """
        Returns the number of nodes between `self` and the seed node of the tree.

        Returns
        -------
        level : integer
            The number of nodes between `self` and the seed node of the tree,
            or 0 if `self` has no parent.
        """
        if self.parent_node:
            return self.parent_node.level() + 1
        else:
            return 0

    def distance_from_root(self):
        """
        Weighted path length of `self` from root.

        Returns
        -------
        dist : numeric
            Total weight of all edges connecting `self` with the root of the
            tree.
        """
        if self.parent_node and self.edge.length != None:
            if self.parent_node.distance_from_root == None:
                return float(self.edge.length)
            else:
                distance_from_root = float(self.edge.length)
                parent_node = self.parent_node
                # The root is identified when a node with no
                # parent is encountered. If we want to use some
                # other criteria (e.g., where a is_root property
                # is True), we modify it here.
                while parent_node:
                    if parent_node.edge.length != None:
                        distance_from_root = distance_from_root + float(parent_node.edge.length)
                    parent_node = parent_node.parent_node
                return distance_from_root
        elif not self.parent_node and self.edge.length != None:
            return float(self.edge.length)
        elif self.parent_node and self.edge.length == None:
            # what do we do here: parent node exists, but my
            # length does not?
            return float(self.parent_node.edge.length)
        elif not self.parent_node and self.edge.length == None:
            # no parent node, and no edge length
            return 0.0
        else:
            # WTF????
            return 0.0

    def distance_from_tip(self):
        """
        Maximum weighted length of path of `self` to tip.

        If tree is not ultrametric (i.e., descendent edges have different
        lengths), then count the maximum of edge lengths. Note that
        :meth:`Tree.calc_node_ages()` is a more efficient way of doing this
        over the whole tree if this value is need for many or all the nodes on
        the tree.

        Returns
        -------
        dist : numeric
            Maximum weight of edges connecting `self` to tip.
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

    ###########################################################################
    ## Representation

    def description(self, depth=1, indent=0, itemize="", output=None, taxon_namespace=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s)" % self.oid
        else:
            label = " (%s: '%s')" % (self.oid, self.label)
        output_strio.write('%s%sNode object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            leader1 = ' ' * (indent + 4)
            leader2 = ' ' * (indent + 8)
            output_strio.write('\n%s[Edge]' % leader1)
            if self.edge is not None:
                edge_desc = self.edge.description(0)
            else:
                edge_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, edge_desc))

            output_strio.write('\n%s[Taxon]' % leader1)
            if self.taxon is not None:
                taxon_desc = self.taxon.description(0)
            else:
                taxon_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, taxon_desc))

            output_strio.write('\n%s[Parent]' % leader1)
            if self.parent_node is not None:
                parent_node_desc = self.parent_node.description(0)
            else:
                parent_node_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, parent_node_desc))
            if hasattr(self.edge, 'split_bitmask'):
                output_strio.write('\n%s[Clade Mask]' % leader1)
                if taxon_namespace is None:
                    output_strio.write('\n%s%s' % (leader2, self.edge.split_bitmask))
                else:
                    output_strio.write('\n%s%s' % (leader2, taxon_namespace.split_bitmask_string(self.edge.split_bitmask)))
            output_strio.write('\n%s[Children]' % leader1)
            if len(self._child_nodes) == 0:
                output_strio.write('\n%sNone' % leader2)
            else:
                for i, cnd in enumerate(self._child_nodes):
                    output_strio.write('\n%s[%d] %s' % (leader2, i, cnd.description(0)))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    ###########################################################################
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
        edge_lengths = not kwargs.get('suppress_edge_lengths', False)
        edge_lengths = kwargs.get('edge_lengths', edge_lengths)
        child_nodes = self.child_nodes()
        if child_nodes:
            out.write('(')
            f_child = child_nodes[0]
            for child in child_nodes:
                if child is not f_child:
                    out.write(',')
                child._write_newick(out, **kwargs)
            out.write(')')

        out.write(self._get_node_str(**kwargs))
        if edge_lengths:
            e = self.edge
            if e:
                sel = e.length
                if sel is not None:
                    fmt = kwargs.get('edge_length_formatter', None)
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

    def _get_node_str(self, **kwargs):
        """returns a string that is an identifier for the node.  This is called
        by the newick-writing functions, so the kwargs that affect how node
        labels show up in a newick string are the same ones used here:
        `suppress_internal_labels` is a Boolean, and defaults to False.
        """
        is_leaf = (len(self._child_nodes) == 0)
        if not is_leaf:
            if kwargs.get("newick", False):
                return self.as_newick_string()
            if kwargs.get("suppress_internal_labels", False) \
                    or not kwargs.get("include_internal_labels", True):
                return ""
        try:
            t = self.taxon
            rt = kwargs.get("reverse_translate")
            if rt:
                tag = rt(t)
            else:
                tag = t.label
        except AttributeError:
            tag = ""
            try:
                tag = self.label
            except AttributeError:
                if not is_leaf:
                    tag = "n{}".format(id(self))
        preserve_spaces = kwargs.get("preserve_spaces", False)
        raw_labels = kwargs.get("raw_labels", False)
        if raw_labels:
            return tag
        elif " " in tag and "_" in tag:
            if "'" in tag:
                tag.replace("'", "''")
            return "'{}'".format(tag)
        else:
            return tag

    ###########################################################################
    ## alternate representation of tree structure for debugging

    def _get_indented_form(self, **kwargs):
        out = StringIO()
        self._write_indented_form(out, **kwargs)
        return out.getvalue()

    def _write_indented_form(self, out, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        split_bitmasks = kwargs.get("splits", True)
        level = kwargs.get("level", 0)
        ancestors = []
        siblings = []
        n = self
        while n is not None:
            n._write_indented_form_line(out, level, **kwargs)
            n, lev = _preorder_list_manip(n, siblings, ancestors)
            level += lev

    def _get_indented_form_line(self, level, **kwargs):
        out = StringIO()
        self._write_indented_form_line(out, level, **kwargs)
        return out.getvalue()

    def _write_indented_form_line(self, out, level, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        label = _format_node(self, **kwargs)
        if kwargs.get("splits"):
            cm = "%s " % _format_split(self.edge.split_bitmask, **kwargs)
        else:
            cm = ""
        out.write("%s%s%s\n" % ( cm, indentation*level, label))


##############################################################################
## Tree

class Tree(taxon.TaxonNamespaceAssociated, base.Readable, base.Writeable):
    """
    An arborescence, i.e. a fully-connected directed acyclic graph with all
    edges directing away from the root and toward the tips. The "root" of the
    tree is represented by the :attr:`Tree.seed_node` attribute.  In unrooted
    trees, this node is an algorithmic artifact. In rooted trees this node is
    semantically equivalent to the root.
    """

    def _parse_from_stream(cls,
            stream,
            schema,
            collection_offset=None,
            tree_offset=None,
            **kwargs):
        """
        Constructs a new :class:`Tree` object and populates it with data from
        file-like object `stream`.

        If the source defines multiple tree collections (e.g. multiple NEXUS
        "Trees" blocks), then the `collection_offset` argument
        can be used to specify the 0-based index of the tree collection, and
        `tree_offset` argument can be used to specify the 0-based
        index of the tree within the collection, as the source. If
        `collection_offset` is not specified or `None`, then all collections in
        the source are merged before considering `tree_offset`.  If
        `tree_offset` is not specified, then the first tree (offset=0) is
        returned.

        Keyword arguments `**kwargs` are passed directly to
        :meth:`TreeList.read()`, which wraps the actual parsing.

        If no tree is found in the source according to the specified criteria,
        then `None` is returned.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`

        collection_offset : integer or None
            0-based index of tree block or collection in source to be parsed.

        tree_offset : integer or None
            0-based index of tree in source to be parsed. If
            `collection_offset` is `None`, then this is the 0-based index of
            the tree across all collections considered together. Otherwise,
            this is the 0-based index within a particular collection. If
            `tree_offset` is `None` or not specified, then the first tree is
            returned.

        \*\*kwargs : keyword arguments
            Arguments to customize parsing and instantiation this :class:`Tree`
            from the data source, including schema- or format-specific
            handling. The following optional keyword arguments are recognized
            and handled by this constructor:

                `label`
                    The label or description of the new :class:`Tree` object.
                `taxon_namespace`
                    Specifies the :class:`TaxonNamespace` object to be
                    attached to the new :class:`Tree` object.

            All other keyword arguments are passed directly to
            :meth:`TreeList.read()`.  Other keyword arguments may be available,
            depending on the implementation of the reader specialized to handle
            `schema` formats. See documentation for details on keyword
            arguments supported by readers of various schemas.


        Returns
        -------
        tree : :class:`Tree` or `None`
            The :class:`Tree` object corresponding to the tree in the data
            source, or `None` if no valid tree description was found.

        """
        taxon_namespace = taxon.process_kwargs_for_taxon_namespace(kwargs, None)
        label = kwargs.pop("label", None)
        reader = dataio.get_reader(schema, **kwargs)
        if collection_offset is None:
            # coerce all tree products into this list
            tree_list = TreeList(taxon_namespace=taxon_namespace)
            reader.read_trees(
                        stream=stream,
                        taxon_namespace_factory=tree_list._taxon_namespace_pseudofactory,
                        tree_list_factory=tree_list._tree_list_pseudofactory,
                        global_annotations_target=None)
        else:
            tree_lists = reader.read_trees(
                        stream=stream,
                        taxon_namespace_factory=tree_list._taxon_namespace_pseudofactory,
                        tree_list_factory=TreeList,
                        global_annotations_target=None)
            tree_list = tree_lists[collection_offset]
        if not tree_list:
            return None
        if tree_offset is None:
            return tree_list[0]
        else:
            return tree_list[tree_offset]
    _parse_from_stream = classmethod(_parse_from_stream)

    def node_factory(cls, **kwargs):
        """
        Creates and returns a :class:`Node` object.

        Derived classes can override this method to provide support for
        specialized or different types of nodes on the tree.

        Parameters
        ----------

        \*\*kwargs : keyword arguments
            Passed directly to constructor of :class:`Node`.

        Returns
        -------
        node : :class:`Node`
            A new :class:`Node` object.

        """
        return Node(**kwargs)
    node_factory = classmethod(node_factory)

    ###########################################################################
    ## Special/Lifecycle methods

    def __init__(self, *args, **kwargs):
        """
        The constructor can optionally construct a :class:`Tree` object by
        cloning another :class:`Tree` object passed as the first positional
        argument, or out of a data source if `stream` and `schema` keyword
        arguments are passed with a file-like object and a schema-specification
        string object values respectively.

        Parameters
        ----------

        \*args : positional argument, optional
            If given, should be exactly one :class:`Tree` object. The new
            :class:`Tree` will then be a structural clone of this argument.

        \*\*kwargs : keyword arguments, optional
            Arguments to customize parsing and instantiation this :class:`Tree`
            from the data source, including schema- or format-specific
            handling. The following optional keyword arguments are recognized
            and handled by this constructor:

                `stream`
                    The file or file-like object to parse.
                `schema`
                    The :term:`schema` or format of the data in `stream`.
                `label`
                    The label or description of the new :class:`Tree` object.
                `taxon_namespace`
                    Specifies the :class:`TaxonNamespace` object to be
                    attached to the new :class:`Tree` object.

            All other keyword arguments are passed directly to :meth:`TreeList.read()`.
            Other keyword arguments may be available, depending on the implementation
            of the reader specialized to handle `schema` formats.


        If `stream` and `schema` keyword arguments are given, will
        construct this :class:`Tree` object from `schema`-formatted source
        given by file-like object `stream`. `schema` must be a
        recognized and tree file schema, such as `nexus`, `newick`, etc,
        for which a specialized tree list writer is available. If this
        is not implemented for the schema specified, then a
        :class:`UnsupportedSchemaError` is raised. Other keywords will be
        passed to the underlying tree parser.

        Tree objects can thus be instantiated in the following ways::

            # /usr/bin/env python

            try:
                from StringIO import StringIO
            except ImportError:
                from io import StringIO
            from dendropy import Tree, TaxonNamespace

            # empty tree
            t1 = Tree()

            # the canonical way to instantiate a Tree from a data source
            # is the use the `get_from_*` family of static factory methods
            t2 = Tree.get_from_stream(open('treefile.tre', 'rU'), "newick", tree_offset=0)
            t3 = Tree.get_from_path('sometrees.nexus',
                    "nexus",
                    collection_offset=2,
                    tree_offset=1)
            s = "((A,B),(C,D));((A,C),(B,D));"
            t4 = Tree.get_from_string(s, "newick") # tree will be '((A,B),(C,D))'
            t5 = Tree.get_from_string(s, "newick", tree_offset=1) # tree will be '((A,C),(B,D))'

            # tree from stream passed to constructor
            t6 = dendropy.Tree(stream=StringIO("((A,B),(C,D));"), schema="newick")

            # passing keywords to underlying tree parser
            t7 = dendropy.Tree.get_from_string(
                      "((A,B),(C,D));",
                      schema="newick",
                      taxon_namespace=t3.taxon_namespace,
                      encode_splits=True)

            # tree structure deep-copied from another tree
            t8 = dendropy.Tree(t7)
            assert t8 is not t7                             # Trees are distinct
            assert t8.symmetric_difference(t7) == 0         # and structure is identical
            assert t8.taxon_namespace is t7.taxon_namespace             # BUT taxa are not cloned.
            nds3 = [nd for nd in t7.postorder_node_iter()]  # Nodes in the two trees
            nds4 = [nd for nd in t8.postorder_node_iter()]  # are distinct objects,
            for i, n in enumerate(nds3):                    # and can be manipulated
                assert nds3[i] is not nds4[i]               # independentally.
            egs3 = [eg for eg in t7.postorder_edge_iter()]  # Edges in the two trees
            egs4 = [eg for eg in t8.postorder_edge_iter()]  # are also distinct objects,
            for i, e in enumerate(egs3):                    # and can also be manipulated
                assert egs3[i] is not egs4[i]               # independentally.
            lves7 = t7.leaf_nodes()                         # Leaf nodes in the two trees
            lves8 = t8.leaf_nodes()                         # are also distinct objects,
            for i, lf in enumerate(lves3):                  # but order is the same,
                assert lves7[i] is not lves8[i]             # and associated Taxon objects
                assert lves7[i].taxon is lves8[i].taxon     # are the same.

            # to create deep copy of a tree with a different taxon set
            taxa = TaxonNamespace()
            t9 = dendropy.Tree(t7, taxon_namespace=taxa)
            assert t9 is not t7                             # As above, the trees are distinct
            assert t9.symmetric_difference(t7) == 0         # and the structures are identical,
            assert t9.taxon_namespace is not t7.taxon_namespace         # but this time, the taxa *are* different
            assert t9.taxon_namespace is taxa                     # as the given TaxonNamespace is used instead.
            lves3 = t7.leaf_nodes()                         # Leaf nodes (and, for that matter other nodes
            lves5 = t9.leaf_nodes()                         # as well as edges) are also distinct objects
            for i, lf in enumerate(lves3):                  # and the order is the same, as above,
                assert lves7[i] is not lves9[i]             # but this time the associated Taxon
                assert lves7[i].taxon is not lves9[i].taxon # objects are distinct though the taxon
                assert lves7[i].taxon.label == lves9[i].taxon.label # labels are the same.

            # can also call `read()` on a Tree object; each read adds the
            # *replaces* the current tree with the definition specified in the
            # data source
            t10 = Tree()
            t10.read(open('boot1.tre', 'rU'), "newick", tree_offset=0)
            t10.read_from_stream(open('boot2.tre', 'rU'), "newick") # same as above
            t10.read_from_string("((A,B),(C,D));((A,C),(B,D));", "newick", tree_offset=0)
            t10.read_from_path("mle.tre", "newick")

            # to 'switch out' the TaxonNamespace of a tree, replace the reference and
            # reindex the taxa:
            t11 = Tree.get_from_string('((A,B),(C,D));', 'newick')
            taxa = TaxonNamespace()
            t11.taxon_namespace = taxa
            t11.reindex_subcomponent_taxa()

        """
        seed_node = kwargs.pop("seed_node", None)
        super(Tree, self).__init__(*args, **kwargs)
        if seed_node is None:
            self.seed_node = self.node_factory()
        else:
            self.seed_node = seed_node
        self.comments = []
        self._is_rooted = None
        self.weight = None
        self.length_type = None

    # def clone_from(self, other):
    #     """
    #     Clones the structure and properties of :class:`Tree` object `other`.

    #     Parameters
    #     ----------
    #     other : :class:`Tree`
    #         Tree object to clone.

    #     Returns
    #     -------
    #     self : :class:`Tree`
    #         Returns `self`.
    #     """
    #     t = copy.deepcopy(other)
    #     for k, v in t.__dict__.iteritems():
    #         if k not in ["_annotations"]:
    #             self.__dict__[k] = v
    #     self.annotations = t.annotations
    #     return self

    # def __deepcopy__(self, memo):
    #     # we treat the taxa as immutable and copy the reference even in a deepcopy
    #     o = TaxonNamespaceLinked.__deepcopy__(self, memo)
    #     memo[id(self)] = o
    #     for k, v in self.__dict__.iteritems():
    #         if k not in ['taxon_namespace', "_annotations"]:
    #             o.__dict__[k] = copy.deepcopy(v, memo)
    #     o.annotations = copy.deepcopy(self.annotations, memo)
    #     memo[id(self.annotations)] = o.annotations
    #     return o

    ###########################################################################
    ## I/O

    def read(self, stream, schema, **kwargs):
        """
        Redefines this :class:`Tree` object based on data from `source`.

        The current :class:`TaxonNamespace` reference will be retained (and modified
        if new operational taxonomic unit concept definitions are
        encountered in the data source), unless a new object or `None`
        is passed using the `taxon_namespace` argument. Note that any metadata
        associated with the tree specified in the source will be *added* to the
        metadata already associated with the current :class:`Tree`. If the current
        tree has any metadata that should not be associated with the tree
        structure being read, call `tree.annotations.drop()` to clear any
        annotations calling this method.

        If the source defines multiple tree collections (e.g. multiple NEXUS
        "Trees" blocks), then the `collection_offset` argument
        can be used to specify the 0-based index of the tree collection, and
        `tree_offset` argument can be used to specify the 0-based
        index of the tree within the collection, as the source. If
        `collection_offset` is not specified or `None`, then all collections in
        the source are merged before considering `tree_offset`.  If
        `tree_offset` is not specified, then the first tree (offset=0) is
        returned.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`

        collection_offset : integer or None
            0-based index of tree block or collection in source to be parsed.

        tree_offset : integer or None
            0-based index of tree in source to be parsed. If
            `collection_offset` is `None`, then this is the 0-based index of
            the tree across all collections considered together. Otherwise,
            this is the 0-based index within a particular collection. If
            `tree_offset` is `None` or not specified, then the first tree is
            returned.

        \*\*kwargs : keyword arguments
            Arguments to customize parsing and instantiation this :class:`Tree`
            from the data source, including schema- or format-specific
            handling. The following optional keyword arguments are recognized
            and handled by this constructor:

                `label`
                    The label or description of the new :class:`Tree` object.
                `taxon_namespace`
                    Specifies the :class:`TaxonNamespace` object to be
                    attached to the new :class:`Tree` object.

            All other keyword arguments are passed directly to
            :meth:`TreeList.read()`.  Other keyword arguments may be available,
            depending on the implementation of the reader specialized to handle
            `schema` formats. See documentation for details on keyword
            arguments supported by readers of various schemas.

        Returns
        -------
        tree : :class:`Tree`
             Returns `self`.

        Raises
        ------
        ValueError
            If no valid trees matching criteria found in source.

        """
        ignore_metadata = kwargs.pop("ignore_metadata", False)
        if "taxon_namespace" not in kwargs and "taxon_set" not in kwargs:
            kwargs["taxon_namespace"] = self.taxon_namespace
        tree = Tree._parse_from_stream(stream, schema, **kwargs)
        if tree is None:
            raise ValueError("Invalid tree source specification")
        self.seed_node = tree.seed_node
        if not ignore_metadata:
            self.annotations.copy_annotations_from(tree)

    def write(self, stream, schema, **kwargs):
        """
        Writes out `self` in `schema` format to a destination given by
        file-like object `stream`.

        Parameters
        ----------
        `schema` : string
            <ust be a recognized and tree file schema, such as "nexus",
            "newick", etc, for which a specialized tree list writer is
            available. If this is not implemented for the schema specified, then
            a UnsupportedSchemaError is raised.

        \*\*kwargs : keyword arguments, optional
            Keyword arguments will be passed directly to the writer for the
            specified schema. See documentation for details on keyword
            arguments supported by writers of various schemas.
        """
        raise NotImplementedError

    ###########################################################################
    ## Node and Edge Collection Access

    def nodes(self, filter_fn=None):
        """
        Returns list of nodes on tree.

        Parameters
        ----------

        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be included in
            the list, or `False` if not. If `filter_fn` is `None` (default),
            then all nodes visited will be included.

        Returns
        -------
        nodes : :py:class:`list` [:class:`Node`]
            List of :class:`Node` objects in the tree.
        """
        nodes = [node for node in self.preorder_node_iter(filter_fn)]
        return nodes

    def leaf_nodes(self):
        """
        Returns list of leaf nodes on the tree.

        Returns
        -------
        nodes : :py:class:`list` [:class:`Node`]
            List of leaf :class:`Node` objects in `self`.
        """
        return [leaf for leaf in self.leaf_node_iter()]

    def internal_nodes(self, exclude_seed_node=False):
        """
        Returns list of internal nodes in the tree.

        Root or seed node is included unless `exclude_seed_node` is `True`.

        Parameters
        ----------
        exclude_seed_node : boolean, optional
            If `False` (default), then the seed node or root is included. If
            `True`, then the seed node is omitted.

        Returns
        -------
        nodes : :py:class:`list` [:class:`Node`]
            List of internal :class:`Node` objects in `self`.
        """
        return [nd for nd in self.preorder_internal_node_iter(exclude_seed_node=exclude_seed_node)]

    def edges(self, filter_fn=None):
        """
        Returns list of edges on tree.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Edge` object as an argument
            and returns `True` if the :class:`Edge` object is to be included,
            or `False` if not. If `filter_fn` is `None` (default), then all
            edges will be included.

        Returns
        -------
        edges : :py:class:`list` [:class:`Edge`]
            List of :class:`Edge` objects in `self`.
        """
        edges = [edge for edge in self.preorder_edge_iter(filter_fn)]
        return edges

    def leaf_edges(self):
        """
        Returns list of leaf edges on the tree.

        Returns
        -------
        edges : :py:class:`list` [:class:`Edge`]
            List of leaf :class:`Edge` objects in `self`.
        """
        return [leaf.edge for leaf in self.leaf_node_iter()]

    def internal_edges(self, exclude_seed_edge=False):
        """
        Returns list of internal edges on tree.

        Parameters
        ----------
        exclude_seed_node : boolean, optional
            If `False` (default), then the edge subtending the seed node or
            root is included. If `True`, then the seed node is omitted.

        Returns
        -------
        edges : :py:class:`list` [:class:`Edge`]
            List of internal :class:`Edge` objects in `self`.
        """
        return [nd.edge for nd in self.preorder_internal_node_iter(exclude_seed_node=exclude_seed_edge)]

    ###########################################################################
    ## Node Finders

    def find_node(self, filter_fn):
        """
        Finds the first node for which `filter_fn(node) == True`.

        For example, if::

            filter_fn = lambda n: hasattr(n, 'genes') and n.genes is not None

        then::

            t.find_node(filter_fn=filter_fn)

        will return all nodes which have an attribute 'genes' and this value
        is not None.

        Parameters
        ----------
        filter_fn : function object
            Takes a single :class:`Node` object as an argument and returns
            `True` if the node should be returned.

        Returns
        -------
        node : :class:`Node` or `None`
            Returns first :class:`Node` object for which the filter function
            `filter_fn` returns `True`, or `None` if no such node exists on
            this tree.
        """
        for node in self.preorder_node_iter(filter_fn):
            return node
        return None

    def find_node_with_label(self, label):
        """
        Returns first node with `label` attribute matching `label` argument.

        Parameters
        ----------
        label : string
            Value for `label` attribute of :class:`Node` object in this tree.

        Returns
        -------
        node : :class:`Node` or `None`
            Returns first :class:`Node` object with `label` attribute having value
            given in `label`, or `None` if no such node is found.

        """
        for node in self.preorder_node_iter():
            if node.label == label:
                return node
        return None

    def find_node_for_taxon(self, taxon):
        """
        Returns node associated with :class:`Taxon` object `taxon`.

        Parameters
        ----------
        taxon : :class:`Taxon` object
            :class:`Taxon` object that should be associated with the node to be
            returned.

        Returns
        -------
        node : :class:`Node` or `None`
            Returns first :class:`Node` object with `taxon` attribute referencing same
            object as `taxon` argument, or `None` if no such node exists.
        """
        for node in self.postorder_node_iter():
            try:
                if node.taxon is taxon:
                    return node
            except AttributeError:
                pass
        return None

    def find_node_with_taxon(self, taxon_filter_fn=None):
        """
        Returns node associated with :class:`Taxon` object for which `taxon_filter_fn`
        returns `True`.

        Parameters
        ----------
        taxon_filter_fn : function object
            Takes a single :class:`Taxon` object as an argument and returns
            `True` if the node associated with that :class:`Taxon` should be
            returned.

        Returns
        -------
        node : :class:`Node` or `None`
            Returns first :class:`Node` object with `taxon` attribute passing filter
            function `taxon_filter_fn`, or `None` if no such node is found.
        """
        for node in self.preorder_node_iter():
            if hasattr(node, "taxon") and node.taxon is not None:
                if taxon_filter_fn(node.taxon):
                    return node
        return None

    def find_node_with_taxon_label(self, label):
        """
        Returns node associated with :class:`Taxon` object with the specified label.

        Parameters
        ----------
        label : string
            Label of :class:`Taxon` object associated with the node to be returned.

        Returns
        -------
        node : :class:`Node` or `None`
            Returns first :class:`Node` object with `taxon` attribute having label
            `label`, or`None` if no such node is found.

        """
        return self.find_node_with_taxon(lambda x: x.label == label)
        # taxon = self.taxon_namespace.get_taxon(label=label)
        # if taxon is None:
        #     return None
        # return self.find_node_with_taxon(lambda x: x is taxon)

    def mrca(self, **kwargs):
        """
        Returns most-recent common ancestor node of a set of taxa on the tree.

        Returns the shallowest node in the tree (the node nearest the tips)
        that has all of the taxa that:

            * are specified by the split bitmask given by the keyword argument
              `split_bitmask`
            * are in the list of Taxon objects given by the keyword argument
              `taxa`
            * have the labels specified by the list of strings given by the
              keyword argument `taxon_labels`

        Returns `None` if no appropriate node is found.  Assumes that edges on
        tree have been decorated with splits hashes. It is possible that split is
        not compatible with the subtree that is returned! (compatibility tests
        are not fully performed).  This function is used to find the
        "insertion point" for a new split via a root to tip search.

        Parameters
        ----------
        \*\*kwargs : keyword arguments
            Exactly one of the following must be specified:

                `split_bitmask` : integer
                    Node object subtended by the first edge compatible with this
                    split will be returned.
                `taxa` : collections.Iterable [:class:`Taxon`]
                    Shallowest node object with descendent nodes associated with
                    all the :class:`Taxon` objects specified will be returned.
                `taxon_labels` : collections.Iterable [string]
                    Shallowest node object with descendent nodes associated
                    with the minimal set of :class:Taxon objects that
                    collectively have all the labels specified in
                    `taxon_labels` will be returned.

            In addition, the following optional keywords are supported:

                `start_node` : :class:`Node`, optional
                    If given, specifies the node at which to start searching.
                    If not, defaults to the root or `seed_node`.

        Returns
        -------
        node : :class:`Node` or `None`
            The most-recent common ancestor of the nodes specified, or `None`
            if no such node exists.
        """
        start_node = kwargs.get("start_node", self.seed_node)
        split_bitmask = None
        if "split_bitmask" in kwargs:
            split_bitmask = kwargs["split_bitmask"]
        else:
            taxa = kwargs.get("taxa", None)
            if taxa is None:
                if "taxon_labels" in kwargs:
                    taxa = self.taxon_namespace.get_taxa(labels=kwargs["taxon_labels"])
                    if len(taxa) != len(kwargs["taxon_labels"]):
                        raise KeyError("Not all labels matched to taxa")
                else:
                    raise TypeError("Must specify one of: 'split_bitmask', 'taxa' or 'taxon_labels'")
            if taxa is None:
                raise ValueError("No taxa matching criteria found")
            split_bitmask = self.taxon_namespace.get_taxa_bitmask(taxa=taxa)

        if split_bitmask is None or split_bitmask == 0:
            raise ValueError("Null split bitmask (0)")

        if not hasattr(start_node.edge, "split_bitmask"):
            treesplit.encode_splits(self, delete_outdegree_one=False)

        if (start_node.edge.split_bitmask & split_bitmask) != split_bitmask:
            return None

        curr_node = start_node
        last_match = start_node
        nd_source = iter(start_node.child_nodes())
        try:
            while True:
                cm = curr_node.edge.split_bitmask
                cms = (cm & split_bitmask)
                if cms:
                    # for at least one taxon cm has 1 and split has 1
                    if cms == split_bitmask:
                        # curr_node has all of the 1's that split has
                        if cm == split_bitmask:
                            return curr_node
                        last_match = curr_node
                        nd_source = iter(curr_node.child_nodes())
                    else:
                        # we have reached a child that has some, but not all of the
                        #   required taxa as descendants, so we return the last_match
                        return last_match
                curr_node = nd_source.next()
        except StopIteration:
            # we shouldn't reach this if all of the descendants are properly
            #   decorated with split_bitmask attributes, but there may be some hacky
            #   context in which we want to allow the function to be called with
            #   leaves that have not been encoded with split_bitmasks.
            return last_match

    ###########################################################################
    ## Node iterators

    def __iter__(self):
        """
        Iterate over nodes on tree in pre-order.

        Example
        -------

        >>> for nd in tree:
        ...    print(nd.label)
        ...

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding the internal nodes of the subtree rooted at
            this node in post-order sequence.
        """
        return self.preorder_node_iter()

    def preorder_node_iter(self, filter_fn=None):
        """
        Pre-order iterator over nodes in tree.

        Visits nodes in `self`, with each node visited before its children.
        Nodes can optionally be filtered by `filter_fn`: only nodes for which
        `filter_fn` returns `True` when called with the node as an argument are
        yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding nodes in `self` in pre-order sequence.
        """
        return self.seed_node.preorder_iter(filter_fn=filter_fn)

    def preorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes in tree.

        Visits internal nodes in `self`, with each node visited before its
        children. In DendroPy, "internal nodes" are nodes that have at least
        one child node, and thus the root or seed node is typically included
        unless `exclude_seed_node` is `True`. Nodes can optionally be filtered
        by `filter_fn`: only nodes for which `filter_fn` returns `True` when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If `False` (default), then the seed node or root is visited. If
            `True`, then the seed node is skipped.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding the internal nodes of `self`.
        """
        return self.seed_node.preorder_internal_node_iter(filter_fn=filter_fn,
                exclude_seed_node=exclude_seed_node)

    def postorder_node_iter(self, filter_fn=None):
        """
        Post-order iterator over nodes of tree.

        Visits nodes in `self`, with each node visited first followed by its
        children. Nodes can optionally be filtered by `filter_fn`: only nodes
        for which `filter_fn` returns `True` when called with the node as an
        argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding the nodes in `self` in post-order sequence.
        """
        return self.seed_node.postorder_iter(filter_fn=filter_fn)

    def postorder_internal_node_iter(self, filter_fn=None, exclude_seed_node=False):
        """
        Pre-order iterator over internal nodes tree.

        Visits internal nodes in `self`, with each node visited after its
        children. In DendroPy, "internal nodes" are nodes that have at least
        one child node, and thus the root or seed node is typically included
        unless `exclude_seed_node` is `True`. Nodes can optionally be filtered
        by `filter_fn`: only nodes for which `filter_fn` returns `True` when
        passed the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.
        exclude_seed_node : boolean, optional
            If `False` (default), then the seed node or root is visited. If
            `True`, then the seed node is skipped.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding the internal nodes of `self` in post-order
            sequence.
        """
        return self.seed_node.postorder_internal_node_iter(filter_fn=filter_fn,
                exclude_seed_node=exclude_seed_node)

    def levelorder_node_iter(self, filter_fn=None):
        """
        Level-order iteration over nodes of tree.

        Visits nodes in `self`, with each node and other nodes at the same
        level (distance from root) visited before their children.  Nodes can
        optionally be filtered by `filter_fn`: only nodes for which `filter_fn`
        returns `True` when called with the node as an argument are visited.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding nodes of `self` in level-order sequence.
        """
        return self.seed_node.levelorder_iter(filter_fn=filter_fn)

    def level_order_node_iter(self, filter_fn=None):
        """
        Deprecated: use :meth:`Tree.levelorder_node_iter()` instead.
        """
        warnings.warn("Use 'levelorder_node_iter()' instead of 'level_order_node_iter()'",
                FutureWarning, stacklevel=2)
        return self.seed_node.levelorder_iter(filter_fn=filter_fn)

    def inorder_node_iter(self, filter_fn=None):
        """
        In-order iteration over nodes of tree.

        Visits nodes in `self`, with each node visited in-between its children.
        Only valid for strictly-bifurcating trees. Nodes can optionally be
        filtered by `filter_fn`: only nodes for which `filter_fn` returns
        `True` when called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding nodes of `self` in infix or in-order sequence.
        """
        return self.seed_node.inorder_iter(filter_fn=filter_fn)

    def leaf_node_iter(self, filter_fn=None):
        """
        Iterate over all tips or leaves of tree.

        Visits all leaf or tip in `self`. Nodes can optionally be filtered by
        `filter_fn`: only nodes for which `filter_fn` returns `True` when
        called with the node as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding leaf nodes in `self`.
        """
        return self.seed_node.leaf_iter(filter_fn=filter_fn)

    def leaf_iter(self, filter_fn=None):
        """
        Deprecated: use :meth:`Tree.leaf_node_iter()` instead.
        """
        warnings.warn("Use 'leaf_node_iter()' instead of 'leaf_iter()'",
                FutureWarning, stacklevel=2)
        return self.seed_node.leaf_iter(filter_fn=filter_fn)

    def ageorder_node_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Iterator over nodes of tree in order of the age of the node (i.e., the
                time since the present).

        Iterates over nodes in order of age ('age' is as given by the `age`
        attribute, which is usually the sum of edge lengths from tips
        to node, i.e., time since present).
        If `include_leaves` is `True` (default), leaves are included in the
        iteration; if `include_leaves` is `False`, leaves will be skipped.
        If `descending` is `False` (default), younger nodes will be returned
        before older ones; if `True`, older nodes will be returned before
        younger ones.

        Parameters
        ----------
        include_leaves : boolean, optional
            If `True` (default), then leaf nodes are included in the iteration.
            If `False`, then leaf nodes are skipped.
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.
        descending : boolean, optional
            If `False` (default), then younger nodes are visited before older
            ones. If `True`, then older nodes are visited before younger ones.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            Iterator over age-ordered sequence of nodes of `self`.
        """
        if self.seed_node.age is None:
            self.calc_node_ages()
        return self.seed_node.ageorder_iter(include_leaves=include_leaves,
                filter_fn=filter_fn,
                descending=descending)

    def age_order_node_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Deprecated: use :meth:`Tree.ageorder_node_iter()` instead.
        """
        warnings.warn("Use 'ageorder_node_iter()' instead of 'age_order_node_iter()'",
                FutureWarning, stacklevel=2)
        return self.ageorder_node_iter(include_leaves=include_leaves,
                filter_fn=filter_fn,
                descending=descending)

    ###########################################################################
    ## Edge iterators

    def preorder_edge_iter(self, filter_fn=None):
        """
        Pre-order iterator over nodes in tree.

        Visits nodes in `self`, with each node visited before its children.
        Nodes can optionally be filtered by `filter_fn`: only nodes for which
        `filter_fn` returns `True` when called with the node as an argument are
        yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Node` object as an argument
            and returns `True` if the :class:`Node` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all nodes visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Node`]
            An iterator yielding nodes in `self` in pre-order sequence.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.preorder_iter(filter_fn=f):
            yield nd.edge

    def preorder_internal_edge_iter(self, filter_fn=None, exclude_seed_edge=False):
        """
        Pre-order iterator over internal edges in tree.

        Visits internal edges in `self`, with each edge visited before its
        children. In DendroPy, "internal edges" are edges that have at least
        one child edge, and thus the root or seed edge is typically included
        unless `exclude_seed_edge` is `True`. Edges can optionally be filtered
        by `filter_fn`: only edges for which `filter_fn` returns `True` when
        passed the edge as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Edge` object as an argument
            and returns `True` if the :class:`Edge` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all edges visited will be yielded.
        exclude_seed_edge : boolean, optional
            If `False` (default), then the edge subtending the seed node or
            root is visited. If `True`, then this edge is skipped.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Edge`]
            An iterator yielding the internal edges of `self`.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.preorder_internal_node_iter(filter_fn=f,
                exclude_seed_node=exclude_seed_edge):
            yield nd.edge

    def postorder_edge_iter(self, filter_fn=None):
        """
        Post-order iterator over edges of tree.

        Visits edges in `self`, with each edge visited first followed by its
        children. Edges can optionally be filtered by `filter_fn`: only edges
        for which `filter_fn` returns `True` when called with the edge as an
        argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Edge` object as an argument
            and returns `True` if the :class:`Edge` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all edges visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Edge`]
            An iterator yielding the edges in `self` in post-order sequence.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.postorder_iter(filter_fn=f):
            yield nd.edge

    def postorder_internal_edge_iter(self, filter_fn=None, exclude_seed_edge=False):
        """
        Pre-order iterator over internal edges tree.

        Visits internal edges in `self`, with each edge visited after its
        children. In DendroPy, "internal edges" are edges that have at least
        one child edge, and thus the root or seed edge is typically included
        unless `exclude_seed_edge` is `True`. Edges can optionally be filtered
        by `filter_fn`: only edges for which `filter_fn` returns `True` when
        passed the edge as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Edge` object as an argument
            and returns `True` if the :class:`Edge` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all edges visited will be yielded.
        exclude_seed_edge : boolean, optional
            If `False` (default), then the seed edge or root is visited. If
            `True`, then the seed edge is skipped.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Edge`]
            An iterator yielding the internal edges of `self` in post-order
            sequence.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.postorder_internal_node_iter(filter_fn=f,
                exclude_seed_node=exclude_seed_edge):
            yield nd.edge

    def levelorder_edge_iter(self, filter_fn=None):
        """
        Level-order iteration over edges of tree.

        Visits edges in `self`, with each edge and other edges at the same
        level (distance from root) visited before their children.  Edges can
        optionally be filtered by `filter_fn`: only edges for which `filter_fn`
        returns `True` when called with the edge as an argument are visited.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Edge` object as an argument
            and returns `True` if the :class:`Edge` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all edges visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Edge`]
            An iterator yielding edges of `self` in level-order sequence.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.levelorder_iter(filter_fn=f):
            yield nd.edge

    def level_order_edge_iter(self, filter_fn=None):
        """
        Deprecated: use :meth:`Tree.levelorder_edge_iter()` instead.
        """
        warnings.warn("Use 'levelorder_edge_iter()' instead of 'level_order_edge_iter()'",
                FutureWarning, stacklevel=2)
        return self.levelorder_edge_iter(filter_fn=filter_fn)

    def inorder_edge_iter(self, filter_fn=None):
        """
        In-order iteration over edges of tree.

        Visits edges in `self`, with each edge visited in-between its children.
        Only valid for strictly-bifurcating trees. Edges can optionally be
        filtered by `filter_fn`: only edges for which `filter_fn` returns
        `True` when called with the edge as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Edge` object as an argument
            and returns `True` if the :class:`Edge` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all edges visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Edge`]
            An iterator yielding edges of `self` in infix or in-order sequence.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.inorder_iter(filter_fn=f):
            yield nd.edge

    def leaf_edge_iter(self, filter_fn=None):
        """
        Iterate over all tips or leaves of tree.

        Visits all leaf or tip in `self`. Edges can optionally be filtered by
        `filter_fn`: only edges for which `filter_fn` returns `True` when
        called with the edge as an argument are yielded.

        Parameters
        ----------
        filter_fn : function object, optional
            A function object that takes a :class:`Edge` object as an argument
            and returns `True` if the :class:`Edge` object is to be yielded by
            the iterator, or `False` if not. If `filter_fn` is `None`
            (default), then all edges visited will be yielded.

        Returns
        -------
        itor : :py:class:`collections.Iterator` [:class:`Edge`]
            An iterator yielding leaf edges in `self`.
        """
        if filter_fn is not None:
            f = lambda x : filter_fn(x.edge)
        else:
            f = None
        for nd in self.seed_node.leaf_iter(filter_fn=f):
            yield nd.edge

    ###########################################################################
    ## Taxa Management

    def infer_taxa(self):
        """
        Returns a new TaxonNamespace object populated with taxa from this
        tree.
        """
        taxon_namespace = TaxonNamespace()
        for node in self.postorder_node_iter():
            if node.taxon is not None and (node.taxon not in taxon_namespace):
                taxon_namespace.add(node.taxon)
        self.taxon_namespace = taxon_namespace
        return taxon_namespace

    def reindex_subcomponent_taxa(self):
        """
        Reassigns node taxon objects
        """
        for node in self.postorder_node_iter():
            t = node.taxon
            if t:
                node.taxon = self.taxon_namespace.require_taxon(label=t.label)

    def unassign_taxa(self, exclude_leaves=False, exclude_internal=False):
        """
        Strips taxon assignments from tree. If `exclude_leaves` is True,
        then taxa on leaves will be retained. If `exclude_internal` is True,
        then taxa on internal nodes will be retained. The `taxon_namespace` is not
        affected by this operation.
        """
        for nd in self.postorder_node_iter():
            if (len(nd._child_nodes) == 0) and not exclude_leaves:
                nd.taxon = None
            elif (len(nd._child_nodes) > 0) and not exclude_internal:
                nd.taxon = None

    def randomly_assign_taxa(self, create_required_taxa=True, rng=None):
        """
        Randomly assigns taxa to leaf nodes. If the number of taxa defined in
        the taxon set of the tree is more than the number of tips, then a random
        subset of taxa in `taxon_namespace` will be assigned to the tips of tree.
        If the number of tips is more than the number of taxa in the `taxon_namespace`,
        and `add_extra_taxa` is not True [default], then new Taxon
        objects will be created and added to the `taxon_namespace`; if `create_required_taxa`
        is False, then an exception is raised.

        In addition, a Random() object or equivalent can be passed using `rng`;
        otherwise GLOBAL_RNG is used.
        """
        if rng is None:
            rng = GLOBAL_RNG
        if len(self.taxon_namespace) == 0:
            for i, nd in enumerate(self.leaf_nodes()):
                nd.taxon = self.taxon_namespace.require_taxon(label=("T%d" % (i+1)))
        else:
            taxa = [t for t in self.taxon_namespace]
            for i, nd in enumerate(self.leaf_nodes()):
                if len(taxa) > 0:
                    nd.taxon = taxa.pop(rng.randint(0, len(taxa)-1))
                else:
                    if not create_required_taxa:
                        raise ValueError("TaxonNamespace has %d taxa, but tree has %d tips" % (len(self.taxon_namespace), len(self.leaf_nodes())))
                    label = "T%d" % (i+1)
                    k = 0
                    while self.taxon_namespace.has_taxon(label=label):
                        label = "T%d" % (i+1+k)
                        k += 1
                    nd.taxon = self.taxon_namespace.require_taxon(label=label)

    ###########################################################################
    ## Structure

    def _get_rooting_state_is_undefined(self):
        return self._is_rooted is None
    rooting_state_is_undefined = property(_get_rooting_state_is_undefined)

    def _get_is_rooted(self):
        return self._is_rooted
    def _set_is_rooted(self, val):
        self._is_rooted = val
    is_rooted = property(_get_is_rooted, _set_is_rooted)

    def _get_is_unrooted(self):
        return not self._is_rooted
    def _set_is_unrooted(self, val):
        self._is_rooted = not val
    is_unrooted = property(_get_is_unrooted, _set_is_unrooted)

    def deroot(self):
        "Converts a degree-2 node at the root to a degree-3 node."
        seed_node = self.seed_node
        if not seed_node:
            return
        child_nodes = seed_node.child_nodes()
        if len(child_nodes) != 2:
            return

        if len(child_nodes[1].child_nodes()) >= 2:
            to_keep, to_del = child_nodes
        elif len(child_nodes[0].child_nodes()) >= 2:
            to_del, to_keep = child_nodes
        else:
            return
        to_del_edge = to_del.edge
        try:
            to_keep.edge.length += to_del_edge.length
        except:
            pass
        to_del_edge.collapse()
        self.is_rooted = False
        return self.seed_node

    def reseed_at(self, new_seed_node, update_splits=False, delete_outdegree_one=True):
        """
        Takes an internal node, `new_seed_node` that must already be in the tree and
        rotates the tree such that `new_seed_node` is the `seed_node` of the tree.
        This is a 'soft' rerooting -- i.e., changes the tree representation so
        tree traversal behaves as if the tree is rooted at 'new_seed_node', but
        it does not actually change the tree's rooting state.  If
        `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        if new_seed_node.is_leaf():
            raise ValueError('Rooting at a leaf is not supported')

        old_par = new_seed_node.parent_node
        if old_par is None:
            return
        full_encode = False
        if update_splits:
            try:
                taxa_mask = self.seed_node.edge.split_bitmask
            except:
                full_encode = True
                update_splits = False
        to_edge_dict = None
        if update_splits:
            to_edge_dict = getattr(self, "split_edges", None)

        if old_par is self.seed_node:
            root_children = old_par.child_nodes()
            if len(root_children) == 2 and delete_outdegree_one:
                # root (old_par) was of degree 2, thus we need to suppress the
                #   node
                fc = root_children[0]
                if fc is new_seed_node:
                    sister = root_children[1]
                else:
                    assert root_children[1] is new_seed_node
                    sister = fc
                if new_seed_node.edge.length:
                    sister.edge.length += new_seed_node.edge.length
                edge_to_del = new_seed_node.edge
                new_seed_node.edge = old_par.edge
                if update_splits:
                    assert new_seed_node.edge.split_bitmask == taxa_mask
                if to_edge_dict:
                    del to_edge_dict[edge_to_del.split_bitmask]
                new_seed_node.add_child(sister, edge_length=sister.edge.length)
                self.seed_node = new_seed_node
                return
        else:
            self.reseed_at(old_par,
                    update_splits=update_splits,
                    delete_outdegree_one=delete_outdegree_one)
        old_par.edge, new_seed_node.edge = new_seed_node.edge, old_par.edge
        e = old_par.edge
        if update_splits:
            if to_edge_dict:
                del to_edge_dict[e.split_bitmask]
            e.split_bitmask = (~(e.split_bitmask)) & taxa_mask
            if to_edge_dict:
                to_edge_dict[e.split_bitmask] = e
            assert new_seed_node.edge.split_bitmask == taxa_mask
        old_par.remove_child(new_seed_node)
        new_seed_node.add_child(old_par, edge_length=e.length)
        self.seed_node = new_seed_node
        if full_encode:
            treesplit.encode_splits(self, delete_outdegree_one=delete_outdegree_one)
        return self.seed_node

    def to_outgroup_position(self, outgroup_node, update_splits=False, delete_outdegree_one=True):
        """Reroots the tree at the parent of `outgroup_node` and makes `outgroup_node` the first child
        of the new root.  This is just a convenience function to make it easy
        to place a clade as the first child under the root.
        Assumes that `outgroup_node` and `outgroup_node.parent_node` and are in the tree/
        If `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        p = outgroup_node.parent_node
        assert p is not None
        self.reseed_at(p, update_splits=update_splits, delete_outdegree_one=delete_outdegree_one)
        p.remove_child(outgroup_node)
        p.add_child(outgroup_node, edge_length=outgroup_node.edge.length, pos=0)
        return self.seed_node

    def reroot_at_node(self, new_root_node, update_splits=False, delete_outdegree_one=True):
        """
        Takes an internal node, `new_seed_node` that must already be in the tree and
        roots the tree at that node.
        This is a 'hard' rerooting -- i.e., changes the tree
        representation so tree traversal behaves as if the tree is rooted at
        'new_seed_node', *and* changes the tree's rooting state.
        If `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        self.reseed_at(new_seed_node=new_root_node,
                update_splits=False,
                delete_outdegree_one=delete_outdegree_one)
        self.is_rooted = True
        if update_splits:
            self.update_splits(delete_outdegree_one=delete_outdegree_one)
        return self.seed_node

    def reroot_at_edge(self,
            edge,
            length1=None,
            length2=None,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Takes an internal edge, `edge`, adds a new node to it, and then roots
        the tree on the new node.
        `length1` and `length2` will be assigned to the new (sub-)edge leading
        to the old parent of the original edge, while `length2` will be
        assigned to the old child of the original edge.
        If `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        old_tail = edge.tail_node
        old_head = edge.head_node
        new_seed_node = old_tail.new_child(edge_length=length1)
        old_tail.remove_child(old_head)
        new_seed_node.add_child(old_head, edge_length=length2)
        self.reroot_at_node(new_seed_node,
                update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)
        return self.seed_node

    def reroot_at_midpoint(self, update_splits=False, delete_outdegree_one=True):
        """
        Reroots the tree at the the mid-point of the longest distance between
        two taxa in a tree.
        Sets the rooted flag on the tree to True.
        If `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        from dendropy import treecalc
        pdm = treecalc.PatristicDistanceMatrix(self)
        n1, n2 = pdm.max_dist_nodes
        plen = float(pdm.max_dist) / 2
        mrca_node = pdm.mrca(n1.taxon, n2.taxon)
        #assert mrca_node is self.mrca(taxa=[n1.taxon, n2.taxon])
        #mrca_node = self.mrca(taxa=[n1.taxon, n2.taxon])
        cur_node = n1

        break_on_node = None # populated *iff* midpoint is exactly at an existing node
        target_edge = None
        head_node_edge_len = None

        # going up ...
        while cur_node is not mrca_node:
            if cur_node.edge.length > plen:
                target_edge = cur_node.edge
                head_node_edge_len = plen #cur_node.edge.length - plen
                plen = 0
                break
            elif cur_node.edge.length < plen:
                plen -= cur_node.edge.length
                cur_node = cur_node.parent_node
            else:
                break_on_node = cur_node

        assert break_on_node is not None or target_edge is not None

        if break_on_node:
            self.reseed_at(break_on_node, update_splits=False, delete_outdegree_one=delete_outdegree_one)
            new_seed_node = break_on_node
        else:
            tail_node_edge_len = target_edge.length - head_node_edge_len
            old_head_node = target_edge.head_node
            old_tail_node = target_edge.tail_node
            old_tail_node.remove_child(old_head_node)
            new_seed_node = Node()
            new_seed_node.add_child(old_head_node, edge_length=head_node_edge_len)
            old_tail_node.add_child(new_seed_node, edge_length=tail_node_edge_len)
            self.reseed_at(new_seed_node, update_splits=False, delete_outdegree_one=delete_outdegree_one)
        self.is_rooted = True
        if update_splits:
            self.update_splits(delete_outdegree_one=False)
        return self.seed_node

    def delete_outdegree_one_nodes(self):
        for nd in self.postorder_node_iter():
            children = nd._child_nodes
            if len(children) == 1:
                if nd.edge.length is not None:
                    if children[0].edge.length is None:
                        children[0].edge.length = nd.edge.length
                    else:
                        children[0].edge.length += nd.edge.length
                if nd.parent_node is not None:
                    parent = nd.parent_node
                    pos = parent._child_nodes.index(nd)
                    parent.remove_child(nd)
                    parent.add_child(children[0], pos=pos)
                    # assert children[0].parent_node is parent
                    # assert children[0] in parent._child_nodes
                    nd.parent_node = None
                else:
                    # assert nd is self.seed_node
                    self.seed_node = children[0]
                    self.seed_node.parent_node = None

    def collapse_unweighted_edges(self,
            threshold=0.0000001,
            update_splits=False):
        """
        Collapse all edges with edge lengths less than or equal to
        ``threshold``.
        """
        for e in self.postorder_edge_iter():
            if e.length <= threshold:
               e.collapse()
        if update_splits:
            self.update_splits()

    def resolve_polytomies(self, update_splits=False, rng=None):
        """
        Arbitrarily resolve polytomies using 0-length splits.

        If `rng` is an object with a sample() method then the polytomy will be
            resolved by sequentially adding (generating all tree topologies
            equiprobably
            rng.sample() should behave like random.sample()
        If `rng` is not passed in, then polytomy is broken deterministically by
            repeatedly joining pairs of children.
        """
        polytomies = []
        for node in self.postorder_node_iter():
            if len(node.child_nodes()) > 2:
                polytomies.append(node)
        for node in polytomies:
            children = node.child_nodes()
            nc = len(children)
            if nc > 2:
                if rng:
                    to_attach = children[2:]
                    for child in to_attach:
                        node.remove_child(child)
                    attachment_points = children[:2] + [node]
                    while len(to_attach) > 0:
                        next_child = to_attach.pop()
                        next_sib = rng.sample(attachment_points, 1)[0]
                        next_attachment = Node()
                        p = next_sib.parent_node
                        p.add_child(next_attachment)
                        p.remove_child(next_sib)
                        next_attachment.add_child(next_sib)
                        next_attachment.add_child(next_child)
                        attachment_points.append(next_attachment)
                else:
                    while len(children) > 2:
                        nn1 = Node()
                        nn1.edge.length = 0
                        c1 = children[0]
                        c2 = children[1]
                        node.remove_child(c1)
                        node.remove_child(c2)
                        nn1.add_child(c1)
                        nn1.add_child(c2)
                        node.add_child(nn1)
                        children = node.child_nodes()
        if update_splits:
            self.update_splits()

    def prune_subtree(self,
            node,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes subtree starting at `node` from tree.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        if node.parent_node is None:
            raise TypeError('Node has no parent and is implicit root: cannot be pruned')
        node.parent_node.remove_child(node)
        if delete_outdegree_one:
            self.delete_outdegree_one_nodes()
        if update_splits:
            self.update_splits()

    def prune_leaves_without_taxa(self,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes all terminal nodes that have their ``taxon`` attribute set to
        ``None``.
        """
        for nd in self.leaf_iter():
            if nd.taxon is None:
                nd.edge.tail_node.remove_child(nd)
        if delete_outdegree_one:
            self.delete_outdegree_one_nodes()
        if update_splits:
            self.update_splits()

    def xprune_taxa(self, taxa, update_splits=False, delete_outdegree_one=True):
        """
        Removes terminal nodes associated with Taxon objects given by the container
        `taxa` (which can be any iterable, including a TaxonNamespace object) from `self`.
        """
        nodes_to_retain = []
        for nd in self.postorder_node_iter():
            if nd.taxon is not None and nd.taxon not in taxa:
                nodes_to_retain.append(nd)
        parent_nodes = []
        nodes_to_retain.append(self.seed_node)
        for nd in list(nodes_to_retain):
            parent_node = nd.parent_node
            while parent_node is not None and parent_node not in nodes_to_retain:
                nodes_to_retain.append(parent_node)
                parent_node = parent_node.parent_node
        # print ">>"
        # for nd in sorted(nodes_to_retain):
        #     print nd.oid
        # print "--"
        to_process = [self.seed_node]
        while to_process:
            nd = to_process.pop(0)
            children = nd._child_nodes
            for ch in children:
                if ch not in nodes_to_retain:
                    nd._child_nodes.remove(ch)
                    # ch.edge.tail_node.remove_child(ch)
            to_process.extend(nd._child_nodes)
        if delete_outdegree_one:
            self.delete_outdegree_one_nodes()
        # print self.as_string("newick")
        # for nd in sorted(self.postorder_node_iter()):
        #     print nd.oid
        # print "<<\n"

    def prune_taxa(self, taxa, update_splits=False, delete_outdegree_one=True):
        """
        Removes terminal nodes associated with Taxon objects given by the container
        `taxa` (which can be any iterable, including a TaxonNamespace object) from `self`.
        """
        nodes_to_remove = []
        for nd in self.postorder_node_iter():
            if nd.taxon and nd.taxon in taxa:
                nd.edge.tail_node.remove_child(nd)
        self.prune_leaves_without_taxa(update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)

    def prune_nodes(self, nodes, prune_leaves_without_taxa=False, update_splits=False, delete_outdegree_one=True):
        for nd in nodes:
            if nd.edge.tail_node is None:
                raise Exception("Attempting to remove root node or node without parent")
            nd.edge.tail_node.remove_child(nd)
        if prune_leaves_without_taxa:
            self.prune_leaves_without_taxa(update_splits=update_splits,
                    delete_outdegree_one=delete_outdegree_one)

    def prune_taxa_with_labels(self,
            labels,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes terminal nodes that are associated with Taxon objects with
        labels given by `labels`.
        """
        taxa = self.taxon_namespace.get_taxa(labels=labels)
        self.prune_taxa(taxa=taxa,
                update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)

    def retain_taxa(self,
            taxa,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes terminal nodes that are not associated with any
        of the Taxon objects given by ``taxa`` (which can be any iterable, including a
        TaxonNamespace object) from the ``self``.
        """
        to_prune = [t for t in self.taxon_namespace if t not in taxa]
        self.prune_taxa(to_prune,
                update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)

    def retain_taxa_with_labels(self,
            labels,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes terminal nodes that are not associated with Taxon objects with
        labels given by `labels`.
        """
        taxa = self.taxon_namespace.get_taxa(labels=labels)
        self.retain_taxa(taxa=taxa,
                update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)

    def randomly_reorient_tree(self, rng=None, update_splits=False):
        """
        Randomly picks a new rooting position and rotates the branches around all
        internal nodes in the `self`. If `update_splits` is True, the the `split_bitmask`
        and `split_edges` attributes kept valid.
        """
        if rng is None:
            rng = GLOBAL_RNG # use the global rng by default
        nd = rng.sample(self.nodes(), 1)[0]
        if nd.is_leaf():
            self.to_outgroup_position(nd, update_splits=update_splits)
        else:
            self.reseed_at(nd, update_splits=update_splits)
        self.randomly_rotate(rng=rng)

    def randomly_rotate(self, rng=None):
        "Randomly rotates the branches around all internal nodes in `self`"
        if rng is None:
            rng = GLOBAL_RNG # use the global rng by default
        internal_nodes = self.internal_nodes()
        for nd in internal_nodes:
            c = nd.child_nodes()
            rng.shuffle(c)
            nd.set_child_nodes(c)

    def ladderize(self, ascending=True):
        """
        Sorts child nodes in ascending (if ``ascending`` is ``False``) or
        descending (if ``ascending`` is ``False``) order in terms of the number of
        children each child node has.
        """
        node_desc_counts = {}
        for nd in self.postorder_node_iter():
            if len(nd._child_nodes) == 0:
                node_desc_counts[nd] = 0
            else:
                total = 0
                for child in nd._child_nodes:
                    total += node_desc_counts[child]
                total += len(nd._child_nodes)
                node_desc_counts[nd] = total
                nd._child_nodes.sort(key=lambda n: node_desc_counts[n], reverse=not ascending)

    def truncate_from_root(self, distance_from_root):
        self.calc_node_root_distances()
        new_terminals = []
        for nd in self.preorder_node_iter():
            if not nd._parent_node:
                # root node
                # TODO: stringictly speaking, this might be a terminal if distance_from_root == 0
                pass
            else:
                if nd.root_distance == distance_from_root:
                    new_terminals.append(nd)
                elif nd.root_distance > distance_from_root and nd._parent_node.root_distance < distance_from_root:
                    # cut above current node
                    nd.edge.length = distance_from_root - nd._parent_node.root_distance
                    nd.root_distance = distance_from_root
                    new_terminals.append(nd)
        for nd in new_terminals:
            for ch in nd.child_nodes():
                nd.remove_child(ch)

    def update_splits(self, **kwargs):
        """
        Recalculates split hashes for tree.
        """
        treesplit.encode_splits(self, **kwargs)

    ###########################################################################
    ## Ages, depths, branch lengths etc. (mutation)

    def scale_edges(self, edge_len_multiplier):
        """Multiplies every edge length in `self` by `edge_len_multiplier`"""
        for e in self.postorder_edge_iter():
            if e.length is not None:
                e.length *= edge_len_multiplier

    def set_edge_lengths_from_node_ages(self, allow_negative_edges=False):
        """
        Sets the edge lengths of the tree so that the path lengths from the
        tips equal the value of the `age` attribute of the nodes.
        """
        for nd in self.preorder_node_iter():
            if nd.parent_node is not None:
                #if nd.parent_node.age < nd.age:
                #    nd.edge.length = 0.0
                #else:
                #    nd.edge.length = nd.parent_node.age - nd.age
                if not allow_negative_edges and nd.parent_node.age < nd.age:
                    #if nd.parent_node is self.seed_node:
                    #    # special case seed node
                    #    nd.parent_node.age = nd.age + nd.edge_length
                    #else:
                    #    raise ValueError('Parent node age (%s: %s) is younger than descendent (%s: %s)'
                    #            % (nd.parent_node.oid, nd.parent_node.age, nd.oid, nd.age))
                    raise ValueError('Parent node age (%s: %s) is younger than descendent (%s: %s)'
                            % (nd.parent_node.oid, nd.parent_node.age, nd.oid, nd.age))
                nd.edge.length = nd.parent_node.age - nd.age

    ###########################################################################
    ## Ages, depths, branch lengths etc. (calculation)

    def calc_node_ages(self, check_prec=0.0000001):
        """
        Adds an attribute called "age" to  each node, with the value equal to
        the sum of edge lengths from the node to the tips. If the lengths of
        different paths to the node differ by more than `check_prec`, then a
        ValueError exception will be raised indicating deviation from
        ultrametricity. If `check_prec` is negative or False, then this check
        will be skipped.
        """
        ages = []
        for node in self.postorder_node_iter():
            ch = node.child_nodes()
            if len(ch) == 0:
               node.age = 0.0
            else:
                first_child = ch[0]
                node.age = first_child.age + first_child.edge.length
                if not (check_prec is None or check_prec is False or check_prec < 0):
                    for nnd in ch[1:]:
                        ocnd = nnd.age + nnd.edge.length
                        if abs(node.age - ocnd) > check_prec:
                            # raise ValueError("Tree is not ultrametric. Node '{}': expecting {}, but found {}".format(node.label, node.age, ocnd))
                            raise ValueError("Tree is not ultrametric")
            ages.append(node.age)
        return ages

    def calc_node_root_distances(self, return_leaf_distances_only=True):
        """
        Adds attribute "root_distance" to each node, with value set to the
        sum of edge lengths from the node to the root. Returns list of
        distances. If `return_leaf_distances_only` is True, then only
        leaf distances will be true.
        """
        dists = []
        for node in self.preorder_node_iter():
            if node._parent_node is None:
                node.root_distance = 0.0
            else:
                node.root_distance = node.edge.length + node._parent_node.root_distance
            if (not return_leaf_distances_only or node.is_leaf()):
                dists.append(node.root_distance)
        return dists

    def node_ages(self, check_prec=0.0000001):
        """
        Returns list of ages of speciation events / coalescence times on tree.
        """
        try:
            ages = [n.age for n in self.internal_nodes()]
        except AttributeError:
            self.calc_node_ages(check_prec=check_prec)
            ages = [n.age for n in self.internal_nodes()]
        ages.sort()
        return ages

    def length(self):
        """
        Returns sum of edge lengths of self. Edges with no lengths defined
        (None) will be considered to have a length of 0.
        Note that we do not overrride `__len__` as this requires an integer
        return value.
        """
        total = 0
        for edge in self.postorder_edge_iter():
            if edge.length is not None:
                total += edge.length
        return total

    def max_distance_from_root(self):
        """
        Returns distance of node furthest from root.
        """
        dists = self.calc_node_root_distances()
        return max(dists)

    def minmax_leaf_distance_from_root(self):
        """
        Returns pair of values, representing the distance of the leaf closest
        to a furthest from the root.
        """
        dists = self.calc_node_root_distances(return_leaf_distances_only=True)
        return min(dists), max(dists)

    def coalescence_intervals(self):
        """
        Returns list of coalescence intervals of self., i.e., the waiting
        times between successive coalescence events.
        """
        ages = self.node_ages()
        intervals = []
        intervals.append(ages[0])
        for i, d in enumerate(ages[1:]):
            intervals.append(d - ages[i])
        return intervals

    def num_lineages_at(self, distance_from_root):
        """
        Returns the number of lineages on the tree at a particular distance
        from the root.
        """
        self.calc_node_root_distances()
        num_lineages = 0
        for nd in self.preorder_node_iter():
            if not nd._parent_node:
                # root node
                pass
            else:
                if nd.root_distance == distance_from_root:
                    num_lineages += 1
                elif nd.root_distance >= distance_from_root and nd._parent_node.root_distance < distance_from_root:
                    num_lineages += 1
        return num_lineages

    ###########################################################################
    ## Metrics -- Unary

    def B1(self):
        """
        Returns the B1 statistic: the reciprocal of the sum of the maximum
        number of nodes between each interior node and tip over all internal
        nodes excluding root.
        """
        b1 = 0.0
        nd_mi = {}
        for nd in self.postorder_node_iter():
            if nd._parent_node is None:
                continue
            child_nodes = nd._child_nodes
            if len(child_nodes) == 0:
                nd_mi[nd] = 0.0
                continue
            mi = max(nd_mi[ch] for ch in child_nodes)
            mi += 1
            nd_mi[nd] = mi
            b1 += 1.0/mi
        return b1

    def colless_tree_imbalance(self, normalize="max"):
        """
        Returns Colless' tree imbalance or I statistic: the sum of differences
        of numbers of children in left and right subtrees over all internal
        nodes. ``normalize`` specifies the normalization:

            * "max" or True [DEFAULT]
                normalized to maximum value for tree of
                this size
            * "yule"
                normalized to the Yule model
            * "pda"
                normalized to the PDA (Proportional to Distinguishable
                Arrangements) model
            * None or False
                no normalization

        """
        colless = 0.0
        num_leaves = 0
        subtree_leaves = {}
        for nd in self.postorder_node_iter():
            if nd.is_leaf():
                subtree_leaves[nd] = 1
                num_leaves += 1
            else:
                total_leaves = 0
                if len(nd._child_nodes) > 2:
                    raise TypeError("Colless' tree imbalance statistic requires strictly bifurcating trees")
                left = subtree_leaves[nd._child_nodes[0]]
                right = subtree_leaves[nd._child_nodes[1]]
                colless += abs(right-left)
                subtree_leaves[nd] = right + left
        if normalize == "yule":
            colless = float(colless - (num_leaves * math.log(num_leaves)) - (num_leaves * (EULERS_CONSTANT - 1.0 - math.log(2))))/num_leaves
        elif normalize == "pda":
            colless = colless / pow(num_leaves, 3.0/2)
        elif normalize is True or normalize == "max":
            ## note that Mooers 1995 (Evolution 49(2):379-384)
            ## remarks that the correct normalization factor is
            ## 2/((num_leaves - 1) * (num_leaves -2))
            colless = colless * (2.0/(num_leaves * (num_leaves-3) + 2))
        elif normalize is not None and normalize is not False:
            raise TypeError("`normalization` accepts only None, True, False, 'yule' or 'pda' as argument values")
        return colless

    def pybus_harvey_gamma(self, prec=0.00001):
        """Returns the gamma statistic of Pybus and Harvey (2000). This statistic
        is used to test for constancy of birth and death rates over the course of
        a phylogeny.  Under the pure-birth process, the statistic should follow
        a standard Normal distibution: a Normal(mean=0, variance=1).

        If the lengths of different paths to the node differ by more than `prec`,
            then a ValueError exception will be raised indicating deviation from
            ultrametricty.
        Raises a Value Error if the tree is not ultrametric, is non-binary, or has
            only 2 leaves.

        As a side effect a `age` attribute is added to the nodes of the self.

        Pybus and Harvey. 2000. "Testing macro-evolutionary models using incomplete
        molecular phylogenies." Proc. Royal Society Series B: Biological Sciences.
        (267). 2267-2272
        """
        # the equation is given by:
        #   T = \sum_{j=2}^n (jg_j)
        #   C = T \sqrt{\frac{1}{12(n-2)}}
        #   C gamma = \frac{1}{n-2}\sum_{i=2}^{n-1} (\sum_{k=2}^i kg_k) - \frac{T}{2}
        # where n is the number of taxa, and g_2 ... g_n is the vector of waiting
        #   times between consecutive (in time, not along a branch) speciation times.
        node = None
        speciation_ages = []
        n = 0
        if self.seed_node.age is None:
            self.calc_node_ages(check_prec=prec)
        for node in self.postorder_node_iter():
            if len(node.child_nodes()) == 2:
                speciation_ages.append(node.age)
            else:
                n += 1
        if node is None:
            raise ValueError("Empty tree encountered")
        speciation_ages.sort(reverse=True)
        g = []
        older = speciation_ages[0]
        for age in speciation_ages[1:]:
            g.append(older - age)
            older = age
        g.append(older)
        if not g:
            raise ValueError("No internal nodes found (other than the root)")
        assert(len(g) == (n - 1))
        T = 0.0
        accum = 0.0
        for i in xrange(2, n):
            list_index = i - 2
            T += i * float(g[list_index])
            accum += T
        list_index = n - 2
        T += (n) * g[list_index]
        nmt = n - 2.0
        numerator = accum/nmt - T/2.0
        C = T*pow(1/(12*nmt), 0.5)
        return numerator/C

    def N_bar(self):
        """
        Returns the $\bar{N}$ statistic: the average number of nodes above a
        terminal node.
        """
        leaf_count = 0
        nbar = 0
        for leaf_node in self.leaf_iter():
            leaf_count += 1
            for parent in leaf_node.ancestor_iter(inclusive=False):
                nbar += 1
        return float(nbar) / leaf_count

    def sackin_index(self, normalize=True):
        """
        Returns the Sackin's index: the sum of the number of ancestors for each
        tip of the tree. The larger the Sackin's index, the less balanced the
        tree. ``normalize`` specifies the normalization:

            * True [DEFAULT]
                normalized to number of leaves; this results in a value
                equivalent to that given by Tree.N_bar()
            * "yule"
                normalized to the Yule model
            * "pda"
                normalized to the PDA (Proportional to Distinguishable
                Arrangements) model
            * None or False
                no normalization

        """
        leaf_count = 0
        num_anc = 0
        for leaf_node in self.leaf_iter():
            leaf_count += 1
            for parent in leaf_node.ancestor_iter(inclusive=False):
                num_anc += 1
        if normalize == "yule":
            x = sum(1.0/j for j in range(2, leaf_count+1))
            s = float(num_anc - (2 * leaf_count * x))/leaf_count
        elif normalize == "pda":
            s = float(num_anc)/(pow(leaf_count, 3.0/2))
        elif normalize is True:
            s = float(num_anc)/leaf_count
        elif normalize is None or normalize is False:
            s = float(num_anc)
        elif normalize is not None and normalize is not False:
            raise TypeError("`normalization` accepts only None, True, False, 'yule' or 'pda' as argument values")
        return s

    def treeness(self):
        """
        Returns the proportion of total tree length that is taken up by
        internal branches.
        """
        internal = 0.0
        external = 0.0
        for nd in self.postorder_node_iter():
            if not nd.parent_node:
                continue
            if nd.is_leaf():
                external += nd.edge.length
            else:
                internal += nd.edge.length
        return internal/(external + internal)

    ###########################################################################
    ## Metrics -- Comparative

    def find_missing_splits(self, other_tree):
        """
        Returns a list of splits that are in self,  but
        not in `other_tree`.
        """
        missing = []
        if self.taxon_namespace is not other_tree.taxon_namespace:
            raise TypeError("Trees have different TaxonNamespace objects: %s vs. %s" \
                    % (hex(id(self.taxon_namespace)), hex(id(other_tree.taxon_namespace))))
        if not hasattr(self, "split_edges"):
            self.encode_splits()
        if not hasattr(other_tree, "split_edges"):
            other_tree.encode_splits()
        for split in self.split_edges:
            if split in other_tree.split_edges:
                pass
            else:
                missing.append(split)
        return missing

    def symmetric_difference(self, other_tree):
        """
        Returns the symmetric_distance between this tree and the tree given by
        `other`, i.e. the sum of splits found in one but not in both trees.
        """
        t = self.false_positives_and_negatives(other_tree)
        return t[0] + t[1]

    def false_positives_and_negatives(self, other_tree):
        """
        Returns a tuple pair: all splits found in `other` but in self, and all
        splits in self not found in other.
        """
        from dendropy import treecalc
        if other_tree.taxon_namespace is not self.taxon_namespace:
            other_tree = Tree(other_tree, taxon_namespace=self.taxon_namespace)
        return treecalc.false_positives_and_negatives(self, other_tree)

    def robinson_foulds_distance(self, other_tree):
        """
        Returns Robinson-Foulds distance between this tree and `other_tree`.
        """
        from dendropy import treecalc
        if other_tree.taxon_namespace is not self.taxon_namespace:
            other_tree = Tree(other_tree, taxon_namespace=self.taxon_namespace)
        return treecalc.robinson_foulds_distance(self, other_tree)

    def euclidean_distance(self, other_tree):
        """
        Returns Euclidean_distance distance between this tree and `other_tree`.
        """
        from dendropy import treecalc
        if other_tree.taxon_namespace is not self.taxon_namespace:
            other_tree = Tree(other_tree, taxon_namespace=self.taxon_namespace)
        return treecalc.euclidean_distance(self, other_tree)

    def _check_children_for_split_compatibility(self, nd_list, split):
        for nd in nd_list:
            if treesplit.is_compatible(nd.edge.split_bitmask, split, self.taxon_set.all_taxa_bitmask()):
                # see if nd has all of the leaves that are flagged as 1 in the split of interest
                if (nd.edge.split_bitmask & split) == split:
                    return nd
            else:
                return None
        return None

    def is_compatible_with_split(self, split):
        nd = self.seed_node
        while True:
            if nd.edge.split_bitmask == split:
                return True
            nd = self._check_children_for_split_compatibility(nd._child_nodes, split)
            if nd is None:
                return False

    def is_compatible_with_tree(self, other):
        pass

    ###########################################################################
    ## Metadata

    def strip_comments(self):
        """
        Remove comments from tree/nodes.
        """
        self.comments = []
        for nd in self.postorder_node_iter():
            nd.comments = []
            nd.edge.comments = []

    ###########################################################################
    ## Representation

    def __str__(self):
        "Dump Newick string."
        return "%s" % self.as_newick_string()

    def __repr__(self):
        return "<Tree object at %s>" % (hex(id(self)))

    def description(self, depth=1, indent=0, itemize="", output=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s)" % self.oid
        else:
            label = " (%s: '%s')" % (self.oid, self.label)
        output_strio.write('%s%sTree object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            newick_str = self.as_newick_string()
            if not newick_str:
                newick_str = "()"
            if depth == 1:
                output_strio.write(': %s' % newick_str)
            elif depth >= 2:
                num_nodes = len([nd for nd in self.preorder_node_iter()])
                num_edges = len([ed for ed in self.preorder_edge_iter()])
                output_strio.write(': %d Nodes, %d Edges' % (num_nodes, num_edges))
                if self.taxon_namespace is not None:
                    output_strio.write("\n%s[Taxon Set]\n" % (" " * (indent+4)))
                    self.taxon_namespace.description(depth=depth-1, indent=indent+8, itemize="", output=output_strio)
                output_strio.write('\n%s[Tree]' % (" " * (indent+4)))
                output_strio.write('\n%s%s' % (" " * (indent+8), newick_str))
                if depth >= 3:
                    output_strio.write("\n%s[Nodes]" % (" " * (indent+4)))
                    for i, nd in enumerate(self.preorder_node_iter()):
                        output_strio.write('\n')
                        nd.description(depth=depth-3, indent=indent+8, itemize="[%d] " % i, output=output_strio, taxon_namespace=self.taxon_namespace)
                    output_strio.write("\n%s[Edges]" % (" " * (indent+4)))
                    for i, ed in enumerate(self.preorder_edge_iter()):
                        output_strio.write('\n')
                        ed.description(depth=depth-3, indent=indent+8, itemize="[%d] " % i, output=output_strio, taxon_namespace=self.taxon_namespace)

        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    def as_python_source(self, tree_obj_name=None, tree_args=None, oids=False):
        """
        Returns string that will rebuild this tree in Python.
        """
        p = []

        if tree_obj_name is None:
            tree_obj_name = "tree_%s" % id(self)
        if self.label is not None:
            label = "'" + self.label + "'"
        else:
            label = "None"
        if oids:
            oid_str = ', oid="%s"' % self.oid
        else:
            oid_str = ""
        if tree_args is None:
            tree_args = ""
        else:
            tree_args = ", " + tree_args
        p.append("%s = dendropy.Tree(label=%s%s%s)" \
            % (tree_obj_name,
               label,
               oid_str,
               tree_args))
        if oids:
            p.append("%s.seed_node.oid = '%s'" % (tree_obj_name, self.seed_node.oid))

        taxon_obj_namer = lambda x: "tax_%s" % id(x)
        for taxon in self.taxon_namespace:
            tobj_name = taxon_obj_namer(taxon)
            if taxon.label is not None:
                label = "'" + taxon.label + "'"
            else:
                label = "None"
            if oids:
                oid_str = ', oid="%s"' % taxon.oid
            else:
                oid_str = ""
            p.append("%s = %s.taxon_namespace.require_taxon(label=%s%s)" \
                % (tobj_name,
                   tree_obj_name,
                   label,
                   oid_str))

        node_obj_namer = lambda x: "nd_%s" % id(x)
        for node in self.preorder_node_iter():
            for child in node.child_nodes():
                if node is self.seed_node:
                    nn = "%s.seed_node" % tree_obj_name
                else:
                    nn = node_obj_namer(node)
                if child.label is not None:
                    label = "'" + child.label + "'"
                else:
                    label = "None"
                if child.taxon is not None:
                    ct = taxon_obj_namer(child.taxon)
                else:
                    ct = "None"
                if oids:
                    oid_str = ', oid="%s"' % child.oid
                else:
                    oid_str = ""
                p.append("%s = %s.new_child(label=%s, taxon=%s, edge_length=%s%s)" %
                        (node_obj_namer(child),
                         nn,
                         label,
                         ct,
                         child.edge.length,
                         oid_str))
                if oids:
                    p.append('%s.edge.oid = "%s"' % (node_obj_namer(child), child.edge.oid))

        return "\n".join(p)

    ###########################################################################
    ## Representation

    def as_ascii_plot(self, **kwargs):
        """
        Returns a string representation a graphic of this tree using ASCII
        characters.

        Keyword arguments:

            `plot_metric`
                A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            `show_internal_node_labels`
                Boolean: whether or not to write out internal node labels.
            `show_internal_node_ids`
                Boolean: whether or not to write out internal node id's.
            `leaf_spacing_factor`
                Positive integer: number of rows between each leaf.
            `display_width`
                Force a particular display width, in terms of number of columns.

        """
        ap = AsciiTreePlot(**kwargs)
        return ap.compose(self)

    def write_ascii_plot(self, stream, **kwargs):
        """
        Writes an ASCII text graphic of this tree to `stream`.

        Keyword arguments:

            `plot_metric`
                A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            `show_internal_node_labels`
                Boolean: whether or not to write out internal node labels.
            `show_internal_node_ids`
                Boolean: whether or not to write out internal node id's.
            `leaf_spacing_factor`
                Positive integer: number of rows between each leaf.
            `display_width`
                Force a particular display width, in terms of number of columns.

        """
        return stream.write(self.as_ascii_plot(**kwargs))

    def print_plot(self, **kwargs):
        """
        Writes an ASCII text graphic of this tree to standard output.

        Keyword arguments:

            ``plot_metric``
                A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            ``show_internal_node_labels``
                Boolean: whether or not to write out internal node labels.
            ``show_internal_node_ids``
                Boolean: whether or not to write out internal node id's.
            ``leaf_spacing_factor``
                Positive integer: number of rows between each leaf.
            ``display_width``
                Force a particular display width, in terms of number of columns.

        """
        import sys
        self.write_ascii_plot(sys.stdout, **kwargs)
        sys.stdout.write("\n")

    def write_as_dot(self, out, **kwargs):
        """Writes the tree to `out` as a DOT formatted digraph"""
        if not kwargs.get("taxon_namespace"):
            kwargs["taxon_namespace"] = self.taxon_namespace
        out.write("digraph G {\n")

        nd_id_to_dot_nd = {}
        for n, nd in enumerate(self.preorder_node_iter()):
            label = _format_node(nd, **kwargs)
            if nd is self.seed_node:
                label = "root %s" % label
            dot_nd = "n%d" % n
            out.write(' %s  [label="%s"];\n' % (dot_nd, label))
            nd_id_to_dot_nd[nd] = dot_nd
        for nd, dot_nd in nd_id_to_dot_nd.iteritems():
            try:
                e = nd.edge
                par_dot_nd = nd_id_to_dot_nd[e.tail_node]
            except:
                pass
            else:
                label = _format_edge(e, **kwargs)
                s = ' %s -> %s [label="%s"];\n' % (par_dot_nd, dot_nd, label)
                out.write(s)
        out.write("}\n")

    ###########################################################################
    ## Debugging/Testing

    def _assign_node_labels_from_taxon(self):
        for nd in self.postorder_node_iter():
            if nd.label is not None:
                continue
            if nd.taxon is not None:
                nd.label = nd.taxon.label

    def _get_indented_form(self, **kwargs):
        out = StringIO()
        self._write_indented_form(out, **kwargs)
        return out.getvalue()

    def write_indented_form(self, out, **kwargs):
        if kwargs.get("splits"):
            if not kwargs.get("taxon_namespace"):
                kwargs["taxon_namespace"] = self.taxon_namespace
        self.seed_node._write_indented_form(out, **kwargs)

    def _debug_check_tree(self, logger_obj=None, **kwargs):
        import logging, inspect
        if logger_obj and logger_obj.isEnabledFor(logging.DEBUG):
            try:
                assert self._debug_tree_is_valid(logger_obj=logger_obj, **kwargs)
            except:
                calling_frame = inspect.currentframe().f_back
                co = calling_frame.f_code
                emsg = "\nCalled from file %s, line %d, in %s" % (co.co_filename, calling_frame.f_lineno, co.co_name)
                _LOG.debug("%s" % str(self))
                _LOG.debug("%s" % self._get_indented_form(**kwargs))
        assert self._debug_tree_is_valid(logger_obj=logger_obj, **kwargs)

    def _debug_tree_is_valid(self, **kwargs):
        """Performs sanity-checks of the tree data structure.

        kwargs:
            `check_splits` if True specifies that the split_edge and split_bitmask attributes
                are checked.
        """
        check_splits = kwargs.get('check_splits', False)
        taxon_namespace = kwargs.get('taxon_namespace')
        if taxon_namespace is None:
            taxon_namespace = self.taxon_namespace
        if check_splits:
            taxa_mask = self.seed_node.edge.split_bitmask
        nodes = {}
        edges = {}
        curr_node = self.seed_node
        assert curr_node.parent_node is None, \
                "{} is seed node, but has non-`None` parent node: {}".format(curr_node, curr_node.parent_node)
        assert curr_node.edge.tail_node is None, \
                "{} is seed node, but edge has non-`None` tail node: {}".format(curr_node, curr_node.edge.parent_node)
        ancestors = []
        siblings = []
        while curr_node:
            assert curr_node not in nodes, \
                    "Node {} seen multiple times".format(curr_node)
            curr_edge = curr_node.edge
            assert curr_edge not in edges, \
                    "Edge of {}, {}, is also an edge of {}".format(curr_node, curr_node.edge, edges[curr_edge])
            edges[curr_edge] = curr_node
            nodes[curr_node] = curr_edge
            assert curr_edge.head_node is curr_node, \
                    "Head node of edge of {}, {}, is {}, not {}".format(curr_node, curr_edge, curr_edge.head_node, curr_node)
            assert curr_edge.tail_node is curr_node.parent_node, \
                    "Tail node of edge of {}, {}, is {}, but parent node is {}".format(curr_node, curr_edge, curr_edge.tail_node, curr_node.parent_node)
            if check_splits:
                cm = 0
                split_bitmask = curr_edge.split_bitmask
                assert (split_bitmask | taxa_mask) == taxa_mask, \
                        "Split mask error: {} (taxa: {})".format(split_bitmask, taxa_mask)
            c = curr_node._child_nodes
            if c:
                for child in c:
                    assert child.parent_node is curr_node, \
                            "Child of {}, {}, has {} as parent".format(curr_node, child, child.parent_node)
                    if check_splits:
                        cm |= child.edge.split_bitmask
            elif check_splits:
                assert curr_node.taxon is not None, \
                        "Cannot check splits: {} is a leaf node, but its `taxon` attribute is `None`".format(curr_node)
                cm = taxon_namespace.taxon_bitmask(curr_node.taxon)
            if check_splits:
                assert (cm & taxa_mask) == split_bitmask, \
                        "Split mask error: {} (taxa: {}, split: {})".format(cm, taxa_mask, split_bitmask)
                assert self.split_edges[split_bitmask] == curr_edge, \
                        "Expecting edge {} for split {}, but instead found {}".format(curr_edge, split_bitmask, self.split_edges[split_bitmask])
            curr_node, level = _preorder_list_manip(curr_node, siblings, ancestors)
        if check_splits:
            for s, e in self.split_edges.iteritems():
                assert e in edges
        return True

    def _as_newick_string(self, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.
        For production purposes, use the the full-fledged 'as_string()'
        method of the object.
        """
        return self.seed_node._as_newick_string(**kwargs)

    def _print_newick(self, **kwargs):
        """
        Convenience method to newick string representation of this tree
        to the standard output stream.
        """
        import sys
        sys.stdout.write(self._as_newick_string(**kwargs))
        sys.stdout.write("\n")

    def _write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.  For
        production purposes, use the the full-fledged 'write_to_stream()'
        method of the object.
        """
        self.seed_node._write_newick(out, **kwargs)

##############################################################################
## TreeList

class TreeList(taxon.TaxonNamespaceAssociated, base.Annotable, base.Readable, base.Writeable):
    """
    A collection of :class:`Tree` objects, all referencing the same "universe" of
    opeational taxonomic unit concepts through the same :class:`TaxonNamespace`
    object reference.
    """

    def _parse_from_stream(cls,
            stream,
            schema,
            collection_offset=None,
            tree_offset=None,
            **kwargs):
        """
        Constructs a new :class:`TreeList` object and populates it with trees from
        file-like object `stream`.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`

        collection_offset : integer or None
            0-based index indicating collection of trees to parse. If `None`,
            then all tree collections are retrieved, with each distinct
            collection parsed into a separate :class:`TreeList` object. If the
            tree colleciton offset index is equal or greater than the number of
            tree collections in the data source, then IndexError is raised.
            Negative offsets work like negative list indexes; e.g., a
            `collection_offset` of -1 means to read the last collection of
            trees in the data source. For data formats that do not support the
            concept of distinct tree collections (e.g. NEWICK) are considered
            single-collection data source (i.e, the only acceptable
            `collection_offset` values are -1 or 0).

        tree_offset : integer or None
            0-based index indicating particular tree within a particular
            collection of trees at which to begin reading.  If not specified or
            `None` (default), then all trees are parsed.  Otherwise, must be an
            integer value up the length of the collection minus 1.  A positive
            offset indicates the number of trees in the collection to skip;
            e.g. a `tree_offset` of 20 means to skip the first 20 trees in the
            collection.  Negative offsets work like negative list indexes;
            e.g., a `tree_offset` value of -10 means to retrieve the last 10
            trees in the collection.  If the tree offset index is equal or
            greater than the number of trees in the collection, then IndexError
            is raised. Requires that a particular tree collection has been
            identified using the `tree_collection_offset` parameter: if
            `tree_collection_offset` is not specified, a `TypeError` is raised.

        \*\*kwargs : keyword arguments
            Arguments to customize parsing, instantiation, processing, and
            accession of :class:`Tree` objects read from the data source, including
            schema- or format-specific handling.

            The following optional keyword arguments are recognized and handled
            by this function:

                * `label` specifies the label or description of the new
                  :class:`TreeList`.
                * `taxon_namespace` specifies the :class:`TaxonNamespace` object to be
                  attached to the new :class:`TreeList` object.
                * `tree_list` : **SPECIAL** If passed a :class:`TreeList` using
                  this keyword, then this instance is populated and returned
                  (instead of a new instance being created).

            All other keyword arguments are passed directly to `TreeList.read()`.
            Other keyword arguments may be available, depending on the implementation
            of the reader specialized to handle `schema` formats.

        Notes
        -----
        Note that in most cases, even if `collection_offset` and `tree_offset`
        are specified to restrict the trees returned, the *entire* data source
        is still parsed and processed. So this is not more efficient than
        reading all the trees and then manually-extracting them later; just
        more convenient. If you need just a single subset of trees from a data
        source, there is no gain in efficiency. If you need multiple trees or
        subsets of trees from the same data source, it would be much more
        efficient to read the entire data source, and extract trees as needed.

        Returns
        -------
        A :class:`TreeList` object.

        """
        reader = dataio.get_reader(schema, **kwargs)
        taxon_namespace = taxon.process_kwargs_for_taxon_namespace(kwargs, None)
        label = kwargs.pop("label", None)

        # Accommodate an existing TreeList object being passed
        tree_list = kwargs.pop("tree_list", None)
        if tree_list is None:
            tree_list = cls(label=label, taxon_namespace=taxon_namespace)

        if collection_offset is None:
            if tree_offset is not None:
                raise TypeError("Cannot specify `tree_offset` without specifying `collection_offset`")
            # coerce all tree products into this list
            reader.read_trees(
                        stream=stream,
                        taxon_namespace_factory=tree_list._taxon_namespace_pseudofactory,
                        tree_list_factory=tree_list._tree_list_pseudofactory,
                        global_annotations_target=None)
        else:
            tree_lists = reader.read_trees(
                        stream=stream,
                        taxon_namespace_factory=tree_list._taxon_namespace_pseudofactory,
                        tree_list_factory=tree_list.__class__,
                        global_annotations_target=None)
            if collection_offset < 0:
                raise IndexError("Collection offset out of range: {} (minimum offset = 0)".format(collection_offset))
            if collection_offset >= len(tree_lists):
                raise IndexError("Collection offset out of range: {} (number of collecitons = {}, maximum offset = {})".format(collection_offset, len(tree_lists), len(tree_lists)-1))
            target_tree_list = tree_lists[collection_offset]
            tree_list.copy_annotations_from(target_tree_list)
            if tree_offset is not None:
                # if tree_offset < 0:
                #     raise IndexError("Tree offset out of range: {} (minimum offset = 0)".format(tree_offset))
                if tree_offset >= len(target_tree_list):
                    raise IndexError("Tree offset out of range: {} (collection size = {}, maximum offset = {})".format(tree_offset, len(target_tree_list), len(target_tree_list)-1))
                for tree in target_tree_list[tree_offset:]:
                    tree_list._trees.append(tree)
            else:
                for tree in target_tree_list:
                    tree_list._trees.append(tree)
        return tree_list
        # taxon_namespace = taxon.process_kwargs_for_taxon_namespace(kwargs, None)
        # label = kwargs.pop("label", None)
        # tree_list = cls(label=label,
        #         taxon_namespace=taxon_namespace)
        # tree_list.read(
        #         stream=stream,
        #         schema=schema,
        #         collection_offset=collection_offset,
        #         tree_offset=tree_offset,
        #         **kwargs)
        # return tree_list
    _parse_from_stream = classmethod(_parse_from_stream)

    def tree_factory(cls, *args, **kwargs):
        """
        Creates and returns a :class:`Tree` of a type that this list undestands how to
        manage.

        Parameters
        ----------
        \*args : positional arguments
            Passed directly to constructor of :class:`Tree`.

        \*\*kwargs : keyword arguments
            Passed directly to constructor of :class:`Tree`.

        Returns
        -------
        A :class:`Tree` object.

        """
        tree = Tree(*args, **kwargs)
        return tree
    tree_factory = classmethod(tree_factory)

    ###########################################################################
    ## Special/Lifecycle methods

    def __init__(self, *args, **kwargs):
        """
        Constructs a new :class:`TreeList` object, populating it with any iterable
        container with Tree object members passed as unnamed argument, or from
        a data source if `stream` and `schema` are passed.

        If passed an iterable container, the objects in that container must be
        of type :class:`Tree` (or derived). If the container is of type :class:`TreeList`,
        then, because each :class:`Tree` object must have the same :class:`TaxonNamespace`
        reference as the containing :class:`TreeList`, the trees in the container
        passed as an initialization argument will be **deep**-copied (except
        for associated :class:`TaxonNamespace` and :class:`Taxon` objects, which will
        be shallow-copied). If the container is any other type of
        iterable, then the :class:`Tree` objects will be **shallow**-copied.

        :class:`TreeList` objects can directly thus be instantiated in the
        following ways::

            # /usr/bin/env python

            import StringIO
            from dendropy import TaxonNamespace, Tree, TreeList

            # instantiate an empty tree
            tlst1 = TreeList()

            # the canonical way to instantiate a TreeList from a data source
            # is `get_from_*` family of static factory methods
            tlst2 = TreeList.get_from_stream(open('treefile.tre', 'rU'), "newick")
            tlst3 = TreeList.get_from_path('sometrees.nexus', "nexus")
            tlst4 = TreeList.get_from_string("((A,B),(C,D));((A,C),(B,D));", "newick")

            # can also call `read()` on a TreeList object; each read adds
            # (appends) the tree(s) found to the TreeList
            tlst5 = TreeList()
            tlst5.read(open('boot1.tre', 'rU'), "newick")
            tlst5.read_from_stream(open('boot2.tre', 'rU'), "newick") # same as above
            tlst5.read_from_string("((A,B),(C,D));((A,C),(B,D));", "newick")
            tlst5.read_from_path("boot3.tre", "newick")

            # populated from list of Tree objects
            tlist6_1 = Tree(stream=StringIO("((A,B),(C,D))"), schema="newick")
            tlist6_2 = Tree(stream=StringIO("((A,C),(B,D))"), schema="newick")
            tlist6 = TreeList([tlist5_1, tlist5_2])

            # tree from data source specified in constructor
            tlst7 = TreeList(stream=StringIO("((A,B),(C,D));((A,C),(B,D));"), schema="newick") # same

            # passing keywords to underlying tree parser
            tlst8 = TreeList(stream=StringIO("((A,B),(C,D));((A,C),(B,D));"),
                             schema="newick",
                             taxon_namespace=tlst3.taxon_namespace,
                             encode_splits=True)

            # deep-copied (but shallow-copy taxa) from another tree list
            tlst9 = TreeList(t4)

            # same
            tlst10 = TreeList([Tree(t) for t in tlst5])

            # Subsets of trees can be read:
            # (Note that in most cases, the entire data source is parsed, so
            # this is not more efficient than reading all the trees and
            # then manually-extracting them later; just more convenient

            # skip the first 100 trees in the first collection of trees
            trees = TreeList.get_from_path("mcmc.tre", "newick",
                        collection_offset=0, tree_offset=100)

            # get the last 10 trees in the first collection of trees
            trees = TreeList.get_from_path("mcmc.tre", "newick",
                        collection_offset=0, tree_offset=-10)

            # get the last 10 trees in the second-to-last collection of trees
            trees = TreeList.get_from_path("mcmc.xml", "nexml",
                        collection_offset=-2, tree_offset=-10)

        """
        super(TreeList, self).__init__(*args, **kwargs)
        self._trees = []

    ###########################################################################
    ## Representation

    # def __str__(self):
    #     return "[{}]".format(", ".join([str(i) for i in self._taxa]))

    # def __repr__(self):
    #     return "<TaxonNamespace {} '{}': [{}]>".format(hex(id(self)), self.label, ", ".join(repr(i) for i in self._taxa))

    # def __hash__(self):
    #     return id(self)

    ###########################################################################
    ## Data I/O

    def _taxon_namespace_pseudofactory(self, *args, **kwargs):
        """
        Dummy factory to coerce all :class:`TaxonNamespace` objects required when
        parsing a data source to reference `self.taxon_namespace`.
        """
        return self.taxon_namespace

    def _tree_list_pseudofactory(self, *args, **kwargs):
        """
        Dummy factory to coerce all :class:`TreeList` objects required when
        parsing a data source to reference `self`.
        """
        return self

    def read(self,
            stream,
            schema,
            collection_offset=None,
            tree_offset=None,
            **kwargs):
        """
        Parses :class:`Tree` objects from data source and adds to this collection.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`.

        collection_offset : integer or None
            0-based index indicating collection of trees to parse. If `None`,
            then all tree collections are retrieved, with each distinct
            collection parsed into a separate :class:`TreeList` object. If the
            tree colleciton offset index is equal or greater than the number of
            tree collections in the data source, then IndexError is raised.
            Negative offsets work like negative list indexes; e.g., a
            `collection_offset` of -1 means to read the last collection of
            trees in the data source. For data formats that do not support the
            concept of distinct tree collections (e.g. NEWICK) are considered
            single-collection data source (i.e, the only acceptable
            `collection_offset` values are -1 or 0).

        tree_offset : integer or None
            0-based index indicating particular tree within a particular
            collection of trees at which to begin reading.  If not specified or
            `None` (default), then all trees are parsed.  Otherwise, must be an
            integer value up the length of the collection minus 1.  A positive
            offset indicates the number of trees in the collection to skip;
            e.g. a `tree_offset` of 20 means to skip the first 20 trees in the
            collection.  Negative offsets work like negative list indexes;
            e.g., a `tree_offset` value of -10 means to retrieve the last 10
            trees in the collection.  If the tree offset index is equal or
            greater than the number of trees in the collection, then IndexError
            is raised. Requires that a particular tree collection has been
            identified using the `tree_collection_offset` parameter: if
            `tree_collection_offset` is not specified, a `TypeError` is raised.

        \*\*kwargs : keyword arguments

            Arguments to customize parsing, instantiation, processing, and
            accession of :class:`Tree` objects read from the data source, including
            schema- or format-specific handling. These will be passed to the
            underlying schema-specific reader for handling.

            General (schema-agnostic) keyword arguments are:

                * `rooted` specifies the default rooting interpretation of the tree.
                * `edge_len_type` specifies the type of the edge lengths (int or
                  float; defaults to 'float')

            Other keyword arguments are available depending on the schema. See
            specific schema handlers (e.g., :class:`NewickReader`, :class:`NexusReader`,
            :class:`NexmlReader`) for more details.

        Notes
        -----
        Note that in most cases, even if `collection_offset` and `tree_offset`
        are specified to restrict the trees read, the *entire* data source
        is still parsed and processed. So this is not more efficient than
        reading all the trees and then manually-extracting them later; just
        more convenient. If you need just a single subset of trees from a data
        source, there is no gain in efficiency. If you need multiple trees or
        subsets of trees from the same data source, it would be much more
        efficient to read the entire data source, and extract trees as needed.

        Returns
        -------
        n : `int`
            The number of :class:`Tree` objects read.

        """
        if "taxon_set" in kwargs or "taxon_namespace" in kwargs:
            raise TypeError("Cannot change `taxon_namespace` when reading into an existing TreeList")
        kwargs["taxon_namespace"] = self.taxon_namespace
        kwargs["tree_list"] = self
        cur_size = len(self._trees)
        TreeList._parse_from_stream(
                stream=stream,
                schema=schema,
                collection_offset=collection_offset,
                tree_offset=tree_offset,
                **kwargs)
        new_size = len(self._trees)
        return new_size - cur_size

    ###########################################################################
    ## List Interface

    # def __cmp__(self, o):
    #     return list.__cmp__(self._taxa, o._taxa)

    def __add__(self, other):
        """
        Creates and returns new :class:`TreeList` with clones of all :class:`Trees` in `self`
        as well as all :class:`Tree` objects in `other`. Note that if `other` is a
        :class:`TreeList`, then the :class:`Trees` are *cloned*; otherwise, they are copied.

        Parameters
        ----------
        other : iterable of :class:`Tree` objects

        Returns
        -------
        :class:`TreeList` object containing clones of :class:`Tree` objects in `self` and
        `other`.
        """
        raise NotImplementedError

    def __contains__(self, tree):
        return tree in self._trees

    def __delitem__(self, tree):
        del self._trees[tree]

    def __eq__(self, other):
        raise NotImplementedError

    def __iter__(self):
        return iter(self._trees)

    def __reversed__(self):
        return reversed(self._trees)

    def __len__(self):
        return len(self._trees)

    def __getitem__(self, tree):
        return self._trees[tree]

    def append(self, tree):
        raise NotImplementedError

    def clear(self):
        raise NotImplementedError

    def extend(self, other):
        raise NotImplementedError

    def index(self, tree):
        raise NotImplementedError

    def insert(self, tree):
        raise NotImplementedError

    def pop(self, tree):
        raise NotImplementedError

    def remove(self, tree):
        raise NotImplementedError

    def reverse(self, tree):
        raise NotImplementedError

    def sort(self, tree):
        raise NotImplementedError

    def new_tree(self, *args, **kwargs):
        tns = taxon.process_kwargs_for_taxon_namespace(kwargs, self.taxon_namespace)
        if tns is not self.taxon_namespace:
            raise TypeError("Cannot create new Tree with different TaxonNamespace")
        kwargs["taxon_namespace"] = self.taxon_namespace
        tree = self.tree_factory(*args, **kwargs)
        self._trees.append(tree)
        return tree

###############################################################################
## AsciiTreePlot

class AsciiTreePlot(object):

    class NullEdgeLengthError(ValueError):
        def __init__(self, *args, **kwargs):
            ValueError.__init__(self, *args, **kwargs)

    def __init__(self, **kwargs):
        """
        __init__ takes the following kwargs:

            * `plot_metric` A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            * `show_internal_node_labels`
                Boolean: whether or not to write out internal node labels.
            * `show_internal_node_ids`
                Boolean: whether or not to write out internal node id's.
            * `leaf_spacing_factor`
                Positive integer: number of rows between each leaf.
            * `display_width`
                Force a particular display width, in terms of number of columns.

        """
        self.plot_metric = kwargs.get('plot_metric', 'depth')
        self.show_internal_node_labels = kwargs.get('show_internal_node_labels', False)
        self.show_internal_node_ids = kwargs.get('show_internal_node_ids', False)
        self.leaf_spacing_factor = kwargs.get('leaf_spacing_factor', 2)
#        self.null_edge_length = kwargs.get('null_edge_length', 0)
        self.display_width = kwargs.get('display_width', None)
        self.reset()

    def reset(self):
        self.grid = []
        self.node_row = {}
        self.node_col = {}
        self.node_offset = {}
        self.current_leaf_row = 0

    def _calc_node_offsets(self, tree):
        if self.plot_metric == 'age' or self.plot_metric == 'depth':

            ## for verification ...
#            tree.calc_node_ages(check_prec=False)
#            for nd in tree.postorder_node_iter():
#                self.node_offset[nd] = nd.age
#            flipped_origin = max(self.node_offset.values())
#            for nd in self.node_offset:
#                self.node_offset[nd] = flipped_origin - self.node_offset[nd]
#            return

            for nd in tree.postorder_node_iter():
                cnds = nd.child_nodes()
                if self.plot_metric == 'depth': # 'number of branchings from tip'
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        depths = [self.node_offset[v] for v in cnds]
                        curr_node_offset = max(depths) + 1
#                        print curr_node_offset, [self.node_offset[v] for v in cnds]
                elif self.plot_metric == 'age': # 'sum of edge weights from tip'
                    # note: no enforcement of ultrametricity!
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        if cnds[0].edge.length is not None:
                            curr_node_offset = self.node_offset[cnds[0]] + cnds[0].edge.length
#                        if len(elens) == 0:
#                            curr_node_offset = self.node_offset[cnds[0]]
#                        else:
#                            curr_node_offset = max(elens) + self.node_offset[cnds[0]]
                else:
                    raise ValueError("Unrecognized plot metric '%s' (must be one of: 'age', 'depth', 'level', or 'length')" % self.plot_metric)
                self.node_offset[nd] = curr_node_offset
            flipped_origin = max(self.node_offset.values())
            for nd in self.node_offset:
                self.node_offset[nd] = flipped_origin - self.node_offset[nd]
        else:
            for nd in tree.preorder_node_iter():
                if self.plot_metric == 'level': # 'number of branchings from root'
                    curr_edge_len = 1
                elif self.plot_metric == 'length': # 'sum of edge weights from root'
                    if nd.edge.length is not None:
                        curr_edge_len = nd.edge.length
                    else:
                        curr_edge_len = 0
                else:
                    raise ValueError("Unrecognized plot metric '%s' (must be one of: 'age', 'depth', 'level', or 'length')" % self.plot_metric)
                if nd.parent_node is None:
                    self.node_offset[nd] = curr_edge_len
                else:
                    self.node_offset[nd] =  curr_edge_len + self.node_offset[nd.parent_node]
#        print "\n".join([str(k) for k in self.node_offset.values()])

    def draw(self, tree, dest):
        dest.write(self.compose(tree))

    def get_label_for_node(self, nd):
        if nd.taxon and nd.taxon.label:
            return nd.taxon.label
        # @TODO: we should have a separate setting for labeling nodes with an
        # id, but thus far when I want to see this, I want
        # internal_nodes_labels too...
        label = []
        if self.show_internal_node_labels and nd.label:
            label.append(nd.label)
        if self.show_internal_node_ids:
            label.append("@")
            label.append(str(id(nd)))
        if not label:
            return "@"
        return "".join(label)

    def compose(self, tree):
        self.reset()
        if self.display_width is None:
            display_width = terminal.terminal_width() - 1
        else:
            display_width = self.display_width
        max_label_len = max([len(self.get_label_for_node(i)) for i in tree.leaf_node_iter()])
        if max_label_len <= 0:
            max_label_len = 0
        #effective_display_width = display_width - max_label_len - len(tree.internal_nodes) - 1
        effective_display_width = display_width - max_label_len - 1
        self._calc_node_offsets(tree)
        widths = [self.node_offset[i] for i in tree.leaf_node_iter() if self.node_offset[i] is not None]
        max_width = float(max(widths))
        if max_width == 0:
            raise AsciiTreePlot.NullEdgeLengthError("Tree cannot be plotted under metric '%s' due to zero or null edge lengths: '%s'" % (self.plot_metric, tree.as_newick_string()))
        edge_scale_factor = float(effective_display_width) / max_width
        self.calc_plot(tree.seed_node,
                       edge_scale_factor=edge_scale_factor)
        for i in range(len(tree.leaf_nodes())*self.leaf_spacing_factor + 1):
            self.grid.append([' ' for i in range(0, display_width)])
        self.draw_node(tree.seed_node)
        display = '\n'.join([''.join(i) for i in self.grid])
        return display

    def calc_plot(self, node, edge_scale_factor):
        """
        First pass through tree, post-order traversal to calculate
        coordinates of each node.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            for n in child_nodes:
                self.calc_plot(n, edge_scale_factor)
            ys = [self.node_row[n] for n in child_nodes]
            self.node_row[node] = int(float((max(ys)-min(ys)) / 2) + min(ys))
        else:
            self.node_row[node] = self.current_leaf_row
            self.current_leaf_row = self.current_leaf_row + self.leaf_spacing_factor
        if node.edge.length is None:
            self.node_col[node] = 1
        else:
            self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)
        self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)

    def draw_label(self, label, row, start_col):
        if label:
            for i in range(len(label)):
                if start_col + i < len(self.grid[row]):
                    self.grid[row][start_col+i] = label[i]

    def draw_node(self, node):
        """
        Second pass through tree, plotting nodes onto given self.grid.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            for i, child_node in enumerate(child_nodes):
                start_row = min([self.node_row[node], self.node_row[child_node]])
                end_row = max([self.node_row[node], self.node_row[child_node]])
                if i == 0:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '/'
                    start_row = start_row+1
                    edge_row = self.node_row[child_node]
                elif i == len(child_nodes)-1:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '\\'
                    edge_row = self.node_row[child_node]
                else:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '+'
                    edge_row = self.node_row[child_node]
                self.draw_node(child_node)
                for x in range(self.node_col[node]+1, self.node_col[child_node]):
                    self.grid[edge_row][x] = '-'
                for y in range(start_row, end_row):
                    self.grid[y][self.node_col[node]] = '|'
            label = []
            if self.show_internal_node_labels or self.show_internal_node_ids:
                label = self.get_label_for_node(node)
                self.draw_internal_text(label, self.node_row[node], self.node_col[node])
            else:
                self.grid[self.node_row[node]][self.node_col[node]]='+'
        else:
            label = self.get_label_for_node(node)
            self.draw_label(label, self.node_row[node], self.node_col[node]+1)

    def draw_internal_text(self, label, r, c):
        row = self.grid[r]
        try:
            for n, letter in enumerate(label):
                row[c + n] = letter
        except:
            pass

###############################################################################
## Helper Functions

def _preorder_list_manip(n, siblings, ancestors):
    """
    Helper function for recursion free preorder traversal, that does
    not rely on attributes of the node other than child_nodes() (thus it
    is useful for debuggging).

    Returns the next node (or None) and the number of levels toward the
    root the function "moved".
    """
    levels_moved = 0
    c = n.child_nodes()
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

def _format_node(nd, **kwargs):
    nf = kwargs.get('node_formatter', None)
    if nf:
        return nf(nd)
    if nd.is_leaf():
        t = nd.taxon
        if t:
            label = t.label
        else:
            label = "anonymous leaf"
    else:
        label = "* %s" % str(nd.oid)
    return label

def _format_edge(e, **kwargs):
    ef = kwargs.get('edge_formatter', None)
    if ef:
        return ef(e)
    return str(e)

def _format_split(split, width=None, **kwargs):
    from dendropy.treesplit import split_as_string
    if width is None:
        width = len(kwargs.get("taxon_namespace"))
    s = split_as_string(split, width, symbol1=kwargs.get("off_symbol"), symbol2=kwargs.get("on_symbol"))
    return s

def _convert_node_to_root_polytomy(nd):
    """If `nd` has two children and at least on of them is an internal node,
    then it will be converted to an out-degree three node (with the edge length
    added as needed).

    Returns a tuple of child nodes that were detached (or() if the tree was not
    modified). This can be useful for removing the deleted node from the split_edges
    dictionary.
    """
    nd_children = nd.child_nodes()
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
        dest_edge_head = nd
    else:
        right_child = nd_children[1]
        dest_edge_head = right_child
    curr_add = None
    if right_child and right_child.is_internal():
        try:
            left_child.edge.length += right_child.edge.length
        except:
            pass
        nd.remove_child(right_child)
        grand_kids = right_child.child_nodes()
        for gc in grand_kids:
            nd.add_child(gc)
        curr_add = right_child
    elif left_child.is_internal():
        try:
            dest_edge_head.edge.length += left_child.edge.length
        except:
            pass
        nd.remove_child(left_child)
        grand_kids = left_child.child_nodes()
        for gc in grand_kids:
            nd.add_child(gc)
        curr_add = left_child
    if curr_add:
        ndl = [curr_add]
        t = _convert_node_to_root_polytomy(nd)
        ndl.extend(t)
        return tuple(ndl)
    return ()

