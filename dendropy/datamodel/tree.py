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

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.datamodel import base
from dendropy.datamodel import taxon
from dendropy import dataio

##############################################################################
## Edge

class Edge(base.Annotable):
    """
    An edge on a tree.
    """

    def __init__(self,
            tail_node=None,
            head_node=None,
            length=None,
            rootedge=False,
            label=None):
        """
        Constructs an edge of a tree.

        Parameters
        ----------

        tail_node : `Node` object
            `Node` object from which this edge originates, i.e., the parent node
            of `head_node`.

        head_node : `Node` object
            `Node` object from at which this edge ends, i.e., the child cnode of
            `tail_node`.

        length : numerical
            A value representing the weight of the edge.

        rootedge : boolean
            Is the child node of this edge the root or seed node of the tree?

        label : string
            Label for this edge.

        """
        base.Annotable.__init__(self)
        self.tail_node = tail_node
        self.head_node = head_node
        self.rootedge = rootedge
        self.length = length

    # def __deepcopy__(self, memo):
    #     o = self.__class__.__new___(self.__class)
    #     memo[id(self)] = o
    #     for k, v in self.__dict__.iteritems():
    #         if not k in ['tail_node', 'head_node', 'length', 'rootedge', "_oid"]:
    #             o.__dict__[k] = copy.deepcopy(v, memo)
    #     return o

    # def collapse(self):
    #     h = self.head_node
    #     if h.is_leaf():
    #         return
    #     t = self.tail_node
    #     if t is None:
    #         return
    #     c = h.child_nodes()
    #     pc = t.child_nodes()
    #     pos = len(pc)
    #     try:
    #         pos = pc.index(h)
    #     except:
    #         pass
    #     for i, ch in enumerate(c):
    #         t.add_child(ch, pos=pos + i)
    #     t.remove_child(h)

    def new_edge(self, *args, **kwargs):
        "Returns a new edge object of the same class of this edge."
        edge = self.__class__(*args, **kwargs)
        return edge

    def is_terminal(self):
        "Returns True if the head node has no children"
        return self.head_node and self.head_node.is_leaf()

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

class Node(base.Annotable):
    """
    A node on a tree.
    """

    ###########################################################################
    ## Life-cycle

    def __init__(self,
            taxon=None,
            label=None,
            edge_length=None):
        """
        Constructs a node.

        Parameters
        ----------

        taxon : `Taxon` object
            The `Taxon` object representing the operational taxonomic unit
            concept associated with this Node.
        label : string
            A label for this node.
        edge_length : numeric
            Length or weight of the edge subtending this node.

        """
        self.taxon = taxon
        self.label = label
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
    ## Iterators

    def preorder_iter(self, filter_fn=None):
        """
        Pre-order traversal of subtree rooted at this node.

        Visits self and all descendant nodes, with each node visited before its
        children. Filtered by `filter_fn`: node is only returned if no
        `filter_fn` is given or if filter_fn returns `True`.

        Parameters
        ----------

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        Returns
        -------
        Iterator over a sequence of nodes resulting from a pre-order traversal
        of the subtree starting at this node.

        """
        stack = [self]
        while stack:
            node = stack.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = node.child_nodes()
            child_nodes.extend(stack)
            stack = child_nodes

    def preorder_internal_node_iter(self, filter_fn=None):
        """
        Pre-order traversal of internal nodes of subtree rooted at this node.

        Visits self and all internal descendant nodes, with any particular node
        visited before its children. Filtered by `filter_fn`: node is only
        returned if no `filter_fn` is given or if filter_fn returns `True`.

        Parameters
        ----------

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        Returns
        -------
        Iterator over a sequence of internal nodes resulting from a pre-order
        traversal of the subtree starting at this node.

        """
        if filter_fn:
            filter_fn = lambda x: (not x.is_leaf() and filter_fn(x)) or None
        else:
            filter_fn = lambda x: (x and not x.is_leaf()) or None
        for node in self.preorder_iter(filter_fn):
            yield node

    def postorder_iter(self, filter_fn=None):
        """
        Post-order traversal of subtree rooted at this node.

        Visits self and all descendant nodes, with amy particular node visited
        after its children. Filtered by `filter_fn`: node is only returned if
        no `filter_fn` is given or if filter_fn returns `True`.

        Parameters
        ----------

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        Returns
        -------
        Iterator over a sequence of nodes resulting from a post-order traversal
        of the subtree starting at this node.

        """
        stack = [(self, False)]
        while stack:
            node, state = stack.pop(0)
            if state:
                if filter_fn is None or filter_fn(node):
                    yield node
            else:
                stack.insert(0, (node, True))
                child_nodes = [(n, False) for n in node.child_nodes()]
                child_nodes.extend(stack)
                stack = child_nodes

    def postorder_internal_node_iter(self, filter_fn=None):
        """
        Post-order traversal of internal nodes of subtree rooted at this node.

        Visits self and all internal descendant nodes, with any particular node
        visited after its children. Filtered by `filter_fn`: node is only
        returned if no `filter_fn` is given or if filter_fn returns `True`.

        Parameters
        ----------

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        Returns
        -------
        Iterator over a sequence of nodes resulting from a post-order traversal
        of the subtree starting at this node.
        """
        if filter_fn:
            filter_fn = lambda x: (not x.is_leaf() and filter_fn(x)) or None
        else:
            filter_fn = lambda x: (x and not x.is_leaf()) or None
        for node in self.postorder_iter(filter_fn):
            yield node

    def level_order_iter(self, filter_fn=None):
        """
        Level-order traversal of subtree rooted at this node.

        Parameters
        ----------

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        Returns
        -------
        Iterator over sequence of nodes of the subtree rooted at this node in
        level-order.

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

    def leaf_iter(self, filter_fn=None):
        """
        Iterate over all leaves that ultimately descend from this node.

        Parameters
        ----------

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        Returns
        -------
        Iterator over a sequence of leaf nodes that have this node as an
        ancestor.

        """
        if filter_fn:
            ff = lambda x: x.is_leaf() and filter_fn(x) or None
        else:
            ff = lambda x: x.is_leaf() and x or None
        for node in self.postorder_iter(ff):
            yield node

    def child_iter(self, filter_fn=None):
        """
        Iterate over all nodes that are the (immediate) children of this node.

        Parameters
        ----------

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        Returns
        -------
        Iterator over a sequence of nodes that have this node as a parent.

        """
        for node in self._child_nodes:
            if filter_fn is None or filter_fn(node):
                yield node

    def ancestor_iter(self, filter_fn=None, inclusive=True):
        """
        Iterates over all ancestors of self.

        Iterate over all nodes that are the ancestors of this node.  If
        `inclusive` is True, self is returned as the first item of the
        sequence.

        Parameters
        ----------

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        inclusive : boolean
            If `True`, includes this node in the sequence.

        Returns
        -------
        Iterator over all predecessor/ancestor nodes of this node.

        """
        if inclusive:
            yield self
        node = self
        while node is not None:
            node = node.parent_node
            if node is not None \
                   and (filter_fn is None or filter_fn(node)):
                yield node

    def age_order_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """

        Iterates over all nodes in subtree rooted at this node in order of age.

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

        include_leaves : boolean
            If `True` (default), then leaf nodes are included in the iteration.
            If `False`, then leaf nodes are skipped.

        filter_fn : function object
            A function object that takes a `Node` object as an argument and
            returns `True` if this node is to be visited during this traversal
            operation.

        descending : boolean
            If `False` (default), then younger nodes are visited before older
            ones. If `True`, then older nodes are visited before younger ones.

        Returns
        -------
        Iterator over age-ordered sequence of nodes in subtree rooted at this
        node.

        """
        if not descending:
            leaves = [nd for nd in self.leaf_iter()]
            queued_pairs = []
            in_queue = set()
            for leaf in leaves:
                age_nd_tuple = (leaf.age, leaf)
                queued_pairs.insert(bisect.bisect(queued_pairs, age_nd_tuple), age_nd_tuple)
                in_queue.add(leaf)
            while queued_pairs:
                next_el = queued_pairs.pop(0)
                age, nd = next_el
                in_queue.remove(nd)
                p = nd.parent_node
                if p and p not in in_queue:
                    age_nd_tuple = (p.age, p)
                    queued_pairs.insert(bisect.bisect(queued_pairs, age_nd_tuple), age_nd_tuple)
                    in_queue.add(p)
                if include_leaves or nd.is_internal():
                    yield nd
        else:
            nds = [(nd.age, nd) for nd in self.preorder_iter()]
            nds.sort(reverse=True)
            for nd in nds:
                if include_leaves or nd[1].is_internal():
                    yield nd[1]

    ###########################################################################
    ## Child Node Access and Manipulation

    def set_child_nodes(self, child_nodes):
        """
        Assigns the set of child nodes for this node.

        Side effects:

            - sets the parent of each child node to this node
            - sets the tail node of each child to self

        Parameters
        ----------

        child_nodes : iterable (of Node objects)

        """
        self._child_nodes = list(child_nodes)
        for nd in self._child_nodes:
            nd.parent_node = self
            nd.edge.tail_node = self

    def set_children(self, child_nodes):
        """Legacy support: delegates to `set_child_nodes()`"""
        return self.set_child_nodes(child_nodes)

    def add_child(self, node, pos=None):
        """
        Adds a child node to this node.

        Results in the parent_node and containing_tree of the node being
        attached set to this node. Returns node that was just attached.

        Parameters
        ----------

        node : Node
            The node to be added as a child of this node.

        pos : integer
            If not `None`, the position in the the sequence of children that
            this child should occupy.

        Returns
        -------
            Node that was added.

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

        Keyword arguments will be passed to `Node` constructor.
        """
        node = self.__class__(**kwargs)
        return self.add_child(node=node)

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
        node : Node object
            The node to be removed.

        suppress_deg_two : boolean
            If `False` (default), no action is taken. If `True`, then if the
            node removal results in a node with degree of two (i.e., a single
            parent and a single child), then it will be removed from
            the tree and its (sole) child will be added as a child of its
            parent (with edge lengths adjusted accordingly).

        Returns
        -------
        Node
            The node removed.

        """
        if not node:
            raise Exception("Tried to remove an non-existing or null node")
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
            raise Exception("Tried to remove a node that is not listed as a child")
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
            raise Exception("Tried to remove an non-existing or null node")
        children = self._child_nodes
        try:
            pos = children.index(node)
        except:
            raise Exception("Tried to remove a node that is not listed as a child")

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
        "Returns the edge subtending this node."
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
        "Returns the length of the edge subtending this node."
        return self._edge.length
    def _set_edge_length(self, v=None):
        """
        Sets the edge subtending this node, and sets head_node of
        `edge` to point to self.
        """
        self._edge.length = v
    edge_length = property(_get_edge_length, _set_edge_length)

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
        "Returns `True` if the node has no child_nodes"
        return bool(not self._child_nodes)

    def is_internal(self):
        "Returns `True` if the node has child_nodes"
        return bool(self._child_nodes)

    def leaf_nodes(self):
        """
        Returns list of all leaf_nodes descended from this node (or just
        list with self as the only member if self is a leaf).

        Note
        ----
        Usage of  `leaf_iter()` is preferable for efficiency reasons unless
        actual list is required.

        """
        return [node for node in \
                self.postorder_iter(lambda x: bool(len(x.child_nodes())==0))]

    def child_nodes(self):
        """
        Returns the a shallow-copy list of all child nodes.

        Note
        ----
        Usage of  `child_iter()` is preferable for efficiency reasons unless
        actual list is required.

        """
        return list(self._child_nodes)

    def incident_edges(self):
        """
        Return parent and child edges.
        """
        e = [c.edge for c in self._child_nodes]
        e.append(self.edge)
        return e

    def get_incident_edges(self):
        """Legacy synonym for 'incident_edges()'"""
        return self.incident_edges()

    def adjacent_nodes(self):
        """Return parent and child nodes."""
        n = [c for c in self._child_nodes]
        if self.parent_node:
            n.append(self.parent_node)
        return n

    def get_adjacent_nodes(self):
        """Legacy synonym for 'get_incident_edges()'"""
        return self.adjacent_nodes()

    def sister_nodes(self):
        """Return all other children of parent, excluding self."""
        p = self.parent_node
        if not p:
            return []
        sisters = [nd for nd in p.child_nodes() if nd is not self]
        return sisters

    ###########################################################################
    ## Metrics

    def level(self):
        "Number of nodes between self and root."
        if self.parent_node:
            return self.parent_node.level() + 1
        else:
            return 0

    def distance_from_root(self):
        """
        Sum of edge lengths from root. Right now, 'root' is taken to
        be a node with no parent node.
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
        Sum of edge lengths from tip to node. If tree is not ultrametric
        (i.e., descendent edges have different lengths), then count the
        maximum of edge lengths. Note that the 'calc_node_ages()' method
        of dendropy.trees.Tree() is a more efficient way of doing this over
        the whole tree.
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
        label = format_node(self, **kwargs)
        if kwargs.get("splits"):
            cm = "%s " % format_split(self.edge.split_bitmask, **kwargs)
        else:
            cm = ""
        out.write("%s%s%s\n" % ( cm, indentation*level, label))


##############################################################################
## Tree

class Tree(taxon.TaxonNamespaceScoped, base.Readable, base.Writeable):
    """
    An arborescence representing a phylogenetic tree.
    Fundamental class that encapsulates functionality and attributes need for
    working with a fully-connected directed acyclic graph (or, more strictly,
    with a root-to-leaf directionality constraint, an "arborescence").
    A `Tree` contains a `seed_node` attribute (from which the entire tree
    springs), which may or may not be the root node. The distinction is
    not consequential in the current implementation, which identifies the root
    node as a node without `child_node` objects.
    """

    def parse_from_stream(cls,
            stream,
            schema,
            collection_offset=None,
            tree_offset=None,
            **kwargs):
        """
        Constructs a new `Tree` object and populates it with data from
        file-like object `stream`.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`

        collection_offset : integer or None
        tree_offset : integer or None

            If the source defines multiple tree collections (e.g. multiple
            NEXUS "Trees" blocks), then the keyword argument
            `collection_offset` can be used to specify the 0-based index of the
            tree collection, and the keyword argument `tree_offset` can be used
            to specify the 0-based index of the tree within the collection, as
            the source. If `collection_offset` is not specified or `None`, then
            all collections in the source are merged before considering
            `tree_offset`.  If `tree_offset` is not specified, then the first
            tree (offset=0) is returned.

        **kwargs : keyword arguments
            Arguments to customize parsing, instantiation, processing, and
            accession of `Tree` objects read from the data source, including
            schema- or format-specific handling.

            The following optional keyword arguments are recognized and handled
            by this function:

                - `label` specifies the label or description of the new
                  `TreeList`.
                - `taxon_namespace` specifies the `TaxonNamespace` object to be
                   attached to the new `Tree` object.

            All other keyword arguments are passed directly to `TreeList.read()`.
            Other keyword arguments may be available, depending on the implementation
            of the reader specialized to handle `schema` formats.

        Returns
        -------
        A `Tree` object (or `None` if not trees were found in the data source).

        """
        taxon_namespace = taxon.TaxonNamespaceScoped.process_kwargs_for_taxon_namespace(kwargs, None)
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
    parse_from_stream = classmethod(parse_from_stream)

    ###########################################################################
    ## Special/Lifecycle methods

    def __init__(self, *args, **kwargs):
        """
        __init__ creates a new Tree object, optionally constructing it by cloning
        another Tree object if this is passed as the first argument, or
        out of a data source if `stream` and `schema` are keyword arguments are
        passed with a file-like object and a schema-specification string object
        values respectively.

        If `stream` and `schema` keyword arguments are given, will
        construct this `Tree` object from `schema`-formatted source
        given by file-like object `stream`. `schema` must be a
        recognized and tree file schema, such as `nexus`, `newick`, etc,
        for which a specialized tree list writer is available. If this
        is not implemented for the schema specified, then a
        `UnsupportedSchemaError` is raised. Other keywords will be
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
        super(Tree, self).__init__(*args, **kwargs)
        self.seed_node = self.new_node()

    ###########################################################################
    ## Node Management

    def new_node(self, *args, **kwargs):
        """
        Creates and returns a `Node` object.

        Does *not* add the new `Node` object to the tree. Derived classes can
        override this method to provide support for specialized or different
        types of nodes on the tree.

        Parameters
        ----------
        *args : positional arguments
            Passed directly to constructor of `Node`.

        **kwargs : keyword arguments
            Passed directly to constructor of `Node`.

        Returns
        -------
        `Node` object.

        """
        return Node(*args, **kwargs)

##############################################################################
## TreeList

class TreeList(taxon.TaxonNamespaceScoped, base.Readable, base.Writeable):
    """
    A collection of `Tree` objects, all referencing the same "universe" of
    opeational taxonomic unit concepts through the same `TaxonNamespace`
    object reference.
    """

    def parse_from_stream(cls,
            stream,
            schema,
            collection_offset=None,
            tree_offset=None,
            **kwargs):
        """
        Constructs a new `TreeList` object and populates it with trees from
        file-like object `stream`.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`

        collection_offset : integer or None
            0-based index indicating collection of trees to parse. If
            `None`, then all tree collections are retrieved, with each distinct
            collection parsed into a separate `TreeList` object. If the
            tree colleciton offset index is equal or greater than the number of
            tree collections in the data source, then IndexError is raised.
            Only applicable when parsing data schemas that support the concept
            of distinct collections of trees (e.g., NEXUS or NeXML).

        tree_offset : integer or None
            0-based index indicating particular tree within a particular
            collection of trees to parse. If `None`, then all trees are parsed.
            Only makes sense (respected) if a particular tree collection has
            been identified using the `tree_collection_offset` parameter. If
            the tree offset index is equal or greater than the number of trees
            in the collection, then IndexError is raised. Only applicable when
            parsing data schemas that support the concept of distinct
            collections of trees (e.g., NEXUS or NeXML).

        **kwargs : keyword arguments
            Arguments to customize parsing, instantiation, processing, and
            accession of `Tree` objects read from the data source, including
            schema- or format-specific handling.

            The following optional keyword arguments are recognized and handled
            by this function:

                - `label` specifies the label or description of the new
                  `TreeList`.
                - `taxon_namespace` specifies the `TaxonNamespace` object to be
                   attached to the new `TreeList` object.

            All other keyword arguments are passed directly to `TreeList.read()`.
            Other keyword arguments may be available, depending on the implementation
            of the reader specialized to handle `schema` formats.

        Returns
        -------
        A `TreeList` object.

        """
        taxon_namespace = taxon.TaxonNamespaceScoped.process_kwargs_for_taxon_namespace(kwargs, None)
        label = kwargs.pop("label", None)
        tree_list = cls(label=label,
                taxon_namespace=taxon_namespace)
        tree_list.read(
                stream=stream,
                schema=schema,
                collection_offset=collection_offset,
                tree_offset=tree_offset,
                **kwargs)
        return tree_list
    parse_from_stream = classmethod(parse_from_stream)

    def tree_factory(cls, *args, **kwargs):
        """
        Creates and returns a `Tree` of a type that this list undestands how to
        manage.

        Parameters
        ----------
        *args : positional arguments
            Passed directly to constructor of `Tree`.

        **kwargs : keyword arguments
            Passed directly to constructor of `Tree`.

        Returns
        -------
        A `Tree` object.

        """
        tree = Tree(*args, **kwargs)
        return tree
    tree_factory = classmethod(tree_factory)

    ###########################################################################
    ## Special/Lifecycle methods

    def __init__(self, *args, **kwargs):
        """
        Constructs a new `TreeList` object, populating it with any iterable
        container with Tree object members passed as unnamed argument, or from
        a data source if `stream` and `schema` are passed.

        If passed an iterable container, the objects in that container must be
        of type `Tree` (or derived). If the container is of type `TreeList`,
        then, because each `Tree` object must have the same `TaxonNamespace`
        reference as the containing `TreeList`, the trees in the container
        passed as an initialization argument will be **deep**-copied (except
        for associated `TaxonNamespace` and `Taxon` objects, which will
        be shallow-copied). If the container is any other type of
        iterable, then the `Tree` objects will be **shallow**-copied.

        `TreeList` objects can directly thus be instantiated in the
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
                             taxon_set=tlst3.taxon_set,
                             encode_splits=True)

            # deep-copied (but shallow-copy taxa) from another tree list
            tlst9 = TreeList(t4)

            # same
            tlst10 = TreeList([Tree(t) for t in tlst5])

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
    #     return hash( (t for t in self._taxa) )

    ###########################################################################
    ## Data I/O

    def _taxon_namespace_pseudofactory(self, *args, **kwargs):
        """
        Dummy factory to coerce all `TaxonNamespace` objects required when
        parsing a data source to reference `self.taxon_namespace`.
        """
        return self.taxon_namespace

    def _tree_list_pseudofactory(self, *args, **kwargs):
        """
        Dummy factory to coerce all `TreeList` objects required when
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
        Parses `Tree` objects from data source and adds to this collection.

        Parameters
        ----------

        stream : file or file-like object
            Source of data.

        schema : string
            Identifier of format of data in `stream`.

        collection_offset : integer or None
            0-based index indicating collection of trees to parse. If
            `None`, then all tree collections are retrieved, with each distinct
            collection parsed into a separate `TreeList` object. If the
            tree colleciton offset index is equal or greater than the number of
            tree collections in the data source, then IndexError is raised.
            Only applicable when parsing data schemas that support the concept
            of distinct collections of trees (e.g., NEXUS or NeXML).

        tree_offset : integer or None
            0-based index indicating particular tree within a particular
            collection of trees to parse. If `None`, then all trees are parsed.
            Only makes sense (respected) if a particular tree collection has
            been identified using the `tree_collection_offset` parameter. If
            the tree offset index is equal or greater than the number of trees
            in the collection, then IndexError is raised. Only applicable when
            parsing data schemas that support the concept of distinct
            collections of trees (e.g., NEXUS or NeXML).

        **kwargs : keyword arguments

            Arguments to customize parsing, instantiation, processing, and
            accession of `Tree` objects read from the data source, including
            schema- or format-specific handling. These will be passed to the
            underlying schema-specific reader for handling.

            General (schema-agnostic) keyword arguments are:

                - `rooted` specifies the default rooting interpretation of the tree.
                - `edge_len_type` specifies the type of the edge lengths (int or
                  float; defaults to 'float')

            Other keyword arguments are available depending on the schema. See
            specific schema handlers (e.g., `NewickReader`, `NexusReader`,
            `NexmlReader`) for more details.

        """
        reader = dataio.get_reader(schema, **kwargs)
        if "taxon_set" in kwargs or "taxon_namespace" in kwargs:
            raise TypeError("Cannot change `taxon_namespace` when reading into an existing TreeList")
        if collection_offset is None:
            assert tree_offset is None
            # coerce all tree products into this list
            reader.read_trees(
                        stream=stream,
                        taxon_namespace_factory=self._taxon_namespace_pseudofactory,
                        tree_list_factory=self._tree_list_pseudofactory,
                        global_annotations_target=None)
        else:
            tree_lists = reader.read_trees(
                        stream=stream,
                        taxon_namespace_factory=self._taxon_namespace_pseudofactory,
                        tree_list_factory=self.__class__,
                        global_annotations_target=None)
            tree_list = tree_lists[collection_offset]
            self.copy_annotations_from(tree_list)
            if tree_offset is not None:
                self.extend(tree_list[tree_offset:])
            else:
                self.extend(tree_lists)

    ###########################################################################
    ## List Interface

    # def __cmp__(self, o):
    #     return list.__cmp__(self._taxa, o._taxa)

    def __add__(self, other):
        """
        Creates and returns new `TreeList` with clones of all `Trees` in `self`
        as well as all `Tree` objects in `other`. Note that if `other` is a
        `TreeList`, then the `Trees` are *cloned*; otherwise, they are copied.

        Parameters
        ----------
        other : iterable of `Tree` objects

        Returns
        -------
        `TreeList` object containing clones of `Tree` objects in `self` and
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
        tree = self.tree_factory(*args, **kwargs)
        self._trees.append(tree)
        return tree
