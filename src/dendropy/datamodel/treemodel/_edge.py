#! /usr/bin/env python
# -*- coding: utf-8 -*-

from io import StringIO
from dendropy.datamodel import basemodel
from dendropy.datamodel.treemodel import _bipartition

class Edge(basemodel.DataObject, basemodel.Annotable):
    """
    An :term:``edge`` on a :term:``tree``.
    """

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------
        head_node : |Node|, optional
            Node from to which this edge links, i.e., the child node of this
            node ``tail_node``.
        length : numerical, optional
            A value representing the weight of the edge.
        rootedge : boolean, optional
            Is the child node of this edge the root or seed node of the tree?
        label : string, optional
            Label for this edge.

        """
        basemodel.DataObject.__init__(self, label=kwargs.pop("label", None))
        self._head_node = kwargs.pop("head_node", None)
        if "tail_node" in kwargs:
            raise TypeError(
                "Setting the tail node directly is no longer supported: instead, set"
                " the parent node of the head node"
            )
        self.rootedge = kwargs.pop("rootedge", None)
        self.length = kwargs.pop("length", None)
        if kwargs:
            raise TypeError("Unsupported keyword arguments: {}".format(kwargs))

        self._bipartition = None
        self.comments = []

    def __copy__(self, memo=None):
        raise TypeError("Cannot directly copy Edge")

    def taxon_namespace_scoped_copy(self, memo=None):
        raise TypeError("Cannot directly copy Edge")

    def __deepcopy__(self, memo=None):
        # call Annotable.__deepcopy__()
        return basemodel.Annotable.__deepcopy__(self, memo=memo)
        # return super(Edge, self).__deepcopy__(memo=memo)

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def __lt__(self, other):
        return id(self) < id(other)

    def _get_tail_node(self):
        if self._head_node is None:
            return None
        return self._head_node._parent_node

    def _set_tail_node(self, node):
        if self._head_node is None:
            raise ValueError("'_head_node' is 'None': cannot assign 'tail_node'")
        # Go through managed property instead of
        # setting attribute to ensure book-keeping
        self._head_node.parent_node = node
    tail_node = property(_get_tail_node, _set_tail_node)

    def _get_head_node(self):
        return self._head_node

    def _set_head_node(self, node):
        # Go through managed property instead of setting attribute to ensure
        # book-keeping; following should also set ``_head_node`` of ``self``
        node.edge = self
    head_node = property(_get_head_node, _set_head_node)

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
        Returns a list of all edges that "share" a node with ``self``.
        """
        he = [i for i in self.head_node.incident_edges() if i is not self]
        te = [i for i in self.tail_node.incident_edges() if i is not self]
        he.extend(te)
        return he
    adjacent_edges = property(get_adjacent_edges)

    def collapse(self, adjust_collapsed_head_children_edge_lengths=False):
        """
        Inserts all children of the head_node of self as children of the
        tail_node of self in the same place in the child_node list that
        head_node had occupied. The edge length and head_node will no longer be
        part of the tree unless ``adjust_collapsed_head_children_edge_lengths``.
        is True.
        """
        to_del = self.head_node
        parent = self.tail_node
        if not parent:
            return
        children = to_del.child_nodes()
        if not children:
            raise ValueError("collapse_self called with a terminal.")
        pos = parent.child_nodes().index(to_del)
        parent.remove_child(to_del)
        for child in children:
            parent.insert_child(pos, child)
            pos += 1
            if adjust_collapsed_head_children_edge_lengths and self.length is not None:
                # print id(child), child.edge.length, self.length
                if child.edge.length is None:
                    child.edge.length = self.length
                else:
                    child.edge.length += self.length

    def invert(self, update_bipartitions=False):
        """
        Changes polarity of edge.
        """
        # self.head_node, self.tail_node = self.tail_node, self.head_node

        if not self.head_node:
            raise ValueError("Cannot invert edge with 'None' for head node")
        if not self.tail_node:
            raise ValueError("Cannot invert edge with 'None' for tail node")

        old_head_node = self.head_node
        new_tail_node = old_head_node
        old_tail_node = self.tail_node
        new_head_node = old_tail_node
        grandparent = old_tail_node._parent_node
        if grandparent is not None:
            for idx, ch in enumerate(grandparent._child_nodes):
                if ch is old_tail_node:
                    grandparent._child_nodes[idx] = old_head_node
                    break
            else:
                # we did not break loop: force insertion of old_head_node if
                # not already there
                if old_head_node not in grandparent._child_nodes:
                    grandparent._child_nodes.append(old_head_node)
        assert old_head_node in old_tail_node._child_nodes
        old_tail_node.remove_child(old_head_node)
        assert old_head_node not in old_tail_node._child_nodes
        old_head_node.add_child(old_tail_node)
        old_tail_node.edge.length, old_head_node.edge.length = (
            old_head_node.edge.length,
            old_tail_node.edge_length,
        )

    def _get_bipartition(self):
        if self._bipartition is None:
            self._bipartition = _bipartition.Bipartition(
                edge=self,
                is_mutable=True,
            )
        return self._bipartition

    def _set_bipartition(self, v=None):
        self._bipartition = v
    bipartition = property(_get_bipartition, _set_bipartition)

    def _get_split_bitmask(self):
        return self.bipartition._split_bitmask

    def _set_split_bitmask(self, h):
        self.bipartition._split_bitmask = h
    split_bitmask = property(_get_split_bitmask, _set_split_bitmask)

    def _get_leafset_bitmask(self):
        return self.bipartition._leafset_bitmask

    def _set_leafset_bitmask(self, h):
        self.bipartition._leafset_bitmask = h
    leafset_bitmask = property(_get_leafset_bitmask, _set_leafset_bitmask)

    def _get_tree_leafset_bitmask(self):
        return self.bipartition._tree_leafset_bitmask

    def _set_tree_leafset_bitmask(self, h):
        self.bipartition._tree_leafset_bitmask = h
    tree_leafset_bitmask = property(
        _get_tree_leafset_bitmask, _set_tree_leafset_bitmask
    )

    def split_as_bitstring(self):
        return self.bipartition.split_as_bitstring()

    def leafset_as_bitstring(self):
        return self.bipartition.leafset_as_bitstring()

    def description(
        self, depth=1, indent=0, itemize="", output=None, taxon_namespace=None
    ):
        """
        Returns description of object, up to level ``depth``.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s, Length=%s)" % (id(self), str(self.length))
        else:
            label = " (%s: '%s', Length=%s)" % (id(self), self.label, str(self.length))
        output_strio.write(
            "%s%sEdge object at %s%s" % (indent * " ", itemize, hex(id(self)), label)
        )
        if depth >= 1:
            leader1 = " " * (indent + 4)
            leader2 = " " * (indent + 8)
            output_strio.write("\n%s[Length]" % leader1)
            if self.length is not None:
                length = self.length
            else:
                length = "None"
            output_strio.write("\n%s%s" % (leader2, length))
            output_strio.write("\n%s[Tail Node]" % leader1)
            if self.tail_node is not None:
                tn = self.tail_node.description(0)
            else:
                tn = "None"
            output_strio.write("\n%s%s" % (leader2, tn))
            output_strio.write("\n%s[Head Node]" % leader1)
            if self.head_node is not None:
                hn = self.head_node.description(0)
            else:
                hn = "None"
            output_strio.write("\n%s%s" % (leader2, hn))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    def _format_edge(self, **kwargs):
        ef = kwargs.get('edge_formatter', None)
        if ef:
            return ef(self)
        return str(self)


