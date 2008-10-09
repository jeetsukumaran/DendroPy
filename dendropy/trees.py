#! /usr/bin/env python

############################################################################
##  trees.py
##
##  Part of the DendroPy phylogenetic computation library.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
This module handles the core definition of tree data structure class,
as well as all the structural classes that make up a tree.
"""

from dendropy import base
from dendropy import taxa

##############################################################################
## TreesBlocks and TreesBlock

class TreesBlock(list, taxa.TaxaLinked):
    """
    Tree manager.
    """

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `elem_id`, `label` and `taxa_block`.
        """
        list.__init__(self, *args)
        taxa.TaxaLinked.__init__(self, *args, **kwargs)
        
    def normalize_taxa(self, taxa_block=None, clear=True):
        """
        Rebuilds taxa block from scratch, or assigns taxon objects from
        given taxa_block based on labels.
        """
        if taxa_block is None:
            taxa_block = self.taxa_block
        if clear:            
            taxa_block.clear()
        for tree in self:
            for node in tree.postorder_node_iter():
                if node.taxon:
                    node.taxon = taxa_block.find_taxon(label=node.taxon.label, update=True)
        taxa_block.sort()
        self.taxa_block = taxa_block
        return taxa_block        

##############################################################################
## Tree

class Tree(base.IdTagged):
    """
    Fundamental class that encapsulates minimal functionality and
    attributes need for working with trees.  A Tree contains a
    seed_node attribute (from which the entire tree springs), which
    may or may not be the root node. The distinction is not
    consequential in the current implementation, which identifies the
    root node as a node without children.
    """

    ## STATIC METHODS #########################################################
    
    def mrca(node1, node2):
        """
        Returns the most-recent common ancestor node of node1 and
        node2.
        """
        mrca_node = None
        for node1_anc in Node.ancestor_iter(node1, inclusive=True):
            for node2_anc in Node.ancestor_iter(node2, inclusive=True):
                if node1_anc == node2_anc:
                    return node1_anc
        return None

    mrca = staticmethod(mrca)

    ## INSTANCE METHODS #######################################################
    
    def __init__(self, elem_id=None, label=None, seed_node=None):
        """
        Initializes a Tree object by defining a base node which must
        be of type `Node` or derived from `Node`.
        """
        base.IdTagged.__init__(self, elem_id=elem_id, label=label)
        self.seed_node = None
        self.length_type = None
        if seed_node is not None:
            self.seed_node = seed_node
        else:
            self.seed_node = Node(elem_id='n0', edge=Edge())
            
    def __str__(self):
        """
        Dump Newick string.
        """
        return self.compose_newick()

    def new_node(self, elem_id=None, label=None):
        """
        Returns a new node object of the class of this tree's seed
        node.
        """
        node = self.seed_node.__class__(elem_id=elem_id,
                                        edge=self.new_edge())
        node.label = label
        return node

    def new_edge(self, elem_id=None):
        """
        Returns a new edge object of the class of this tree's seed
        node's edge.
        """
        edge = self.seed_node.new_edge()
        edge.elem_id = elem_id
        return edge

    ## Easy access to seed_node edge ##

    def _get_seed_edge(self):
        """
        Returns the edge of the base node.
        """
        return self.seed_node.edge
    seed_edge = property(_get_seed_edge)

    ## Convenience Methods ##

    def nodes(self, cmp_fn=None, filter_fn=None):
        """
        Returns list of nodes on the tree, sorted using cmp_fn.
        """
        nodes = [node for node in self.preorder_node_iter(filter_fn)]
        if cmp_fn:
            nodes.sort(cmp_fn)
        return nodes
    
    def leaves(self):
        """
        Returns list of leaves on the tree.
        """
        return [leaf for leaf in self.leaf_iter()]

    def find_node(self, elem_id):
        """
        Finds the first node with matching id.
        """
        filter_fn = lambda x: x.elem_id == elem_id
        found = [node for node in self.preorder_node_iter(filter_fn)]
        if found and len(found) > 0:
            return found[0]
        else:
            return None

    def find_edge(self, elem_id):
        """
        Finds the first edge with matching id.
        """
        filter_fn = lambda x: x.elem_id == elem_id
        found = [edge for edge in self.preorder_edge_iter(filter_fn)]
        if found and len(found) > 0:
            return found[0]
        else:
            return None
    
    ## Node iterators ##

    def preorder_node_iter(self, filter_fn=None):
        """
        Returns preorder iterator over tree nodes.
        """
        for node in self.seed_node.preorder_iter(self.seed_node, filter_fn):
            yield node

    def postorder_node_iter(self, filter_fn=None):
        """
        Returns postorder iterator over tree nodes.
        """
        for node in self.seed_node.postorder_iter(self.seed_node, filter_fn):
            yield node

    def level_order_node_iter(self, filter_fn=None):
        """
        Returns level-order iterator over tree nodes.
        """
        for node in self.seed_node.level_order_iter(self.seed_node, filter_fn):
            yield node

    def leaf_iter(self, filter_fn=None):
        """
        Returns an iterator over tree leaves (order determined by
        postorder tree-traversal).
        """
        for node in self.seed_node.leaf_iter(self.seed_node, filter_fn):
            yield node

    ## Edge iterators ##

    def preorder_edge_iter(self, filter_fn=None):
        """
        Returns preorder iterator over tree edges.
        """
        for node in self.seed_node.preorder_iter(self.seed_node):
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def postorder_edge_iter(self, filter_fn=None):
        """
        Returns postorder iterator over tree edges.
        """
        for node in self.seed_node.postorder_iter(self.seed_node):
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def level_order_edge_iter(self, filter_fn=None):
        """
        Returns level-order iterator over tree edges.
        """
        for node in self.seed_node.level_order_iter(self.seed_node):
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge
                
    ## Taxa ##
    
    def infer_taxa_block(self):
        """
        Returns a new TaxaBlock object populated with taxa from this
        tree.
        """
        taxa_block = taxa.TaxaBlock()
        for node in self.postorder_node_iter():
            if node.taxon and (node.taxon not in taxa_block):
                taxa_block.append(node.taxon)
        taxa_block.sort()
        return taxa_block
        
    def normalize_taxa(self, taxa_block, update_taxa_block=True):
        """
        Reassigns tree taxa objects to corresponding taxa objects in
        given taxa_block, with identity of taxa objects determined by
        labels.
        """
        for node in self.postorder_node_iter():
            if node.taxon:
                node.taxon = taxa_block.find_taxon(label=node.taxon.label, update=update_taxa_block)                    
                
    ## for debugging ##
    def compose_newick(self, include_internal_labels=True):
        return self.seed_node.compose_newick(include_internal_labels=include_internal_labels)
                
    ## basic tree manipulation ##
    def deroot(self):
        if self.seed_node:
            children = self.seed_node.children()
            if children and len(children) == 2:
                if len(children[0].children()) >= 2:
                    new_child = children[0]
                    new_seed = children[1]
                else: #lif len(children[1].children()) >= 2:
                    new_child = children[1]
                    new_seed = children[0]
                new_edge_length = 0.0
                if new_child.edge.length:
                    new_edge_length += new_child.edge.length
                if new_seed.edge.length:
                    new_edge_length += new_seed.edge.length
                    new_seed.edge = None
                self.seed_node = new_seed
                self.seed_node.add_child(new_child)
        
##############################################################################
## Node

class Node(taxa.TaxonLinked):
    """
    A node on a tree, implementing only fundamental behaviour and
    properties.
    """

    ### ITERATORS ############################################################

    def preorder_iter(node, filter_fn=None):
        """
        Preorder traversal of the node and its children.  Returns node
        and all descendants such that node is returned before node's
        children (and their children). Filtered by filter_fn: node is
        only returned if no filter_fn is given or if filter_fn returns
        True.
        """
        if not node:
            return
        stack = [node]
        while stack:
            node = stack.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            children = node.children()
            children.extend(stack)
            stack = children
    preorder_iter = staticmethod(preorder_iter)

    def postorder_iter(node, filter_fn=None):
        """
        Postorder traversal of the node and its children.  Returns node
        and all descendants such that node's children (and their
        children) are visited before node.  Filtered by filter_fn:
        node is only returned if no filter_fn is given or if filter_fn
        returns True.
        """
        stack = [(node, False)]
        while stack:
            node, state = stack.pop(0)
            if state:
                if filter_fn is None or filter_fn(node):
                    yield node
            else:
                stack.insert(0, (node, True))
                children = [(n, False) for n in node.children()]
                children.extend(stack)
                stack = children
    postorder_iter = staticmethod(postorder_iter)

    def leaf_iter(start_nd, filter_fn=None):
        """
        Returns an iterator over the leaves that are descendants `of start_nd`
        (order determined by postorder tree-traversal).
        """
        if filter_fn:
            filter_fn = lambda x: x.is_leaf() and filter_fn(x) or None
        else:
            filter_fn = lambda x: x.is_leaf() and x or None
        for node in start_nd.postorder_iter(start_nd, filter_fn):
            yield node
            
    leaf_iter = staticmethod(leaf_iter)

    def level_order_iter(node, filter_fn=None):
        """
        Level-order traversal of the node and its children. Filtered
        by filter_fn: node is only returned if no filter_fn is given
        or if filter_fn returns True
        """
        if filter_fn is None or filter_fn(node):
            yield node
        remaining = node.children()
        while len(remaining) > 0:
            node = remaining.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            children = node.children()
            remaining.extend(children)
            
    level_order_iter = staticmethod(level_order_iter)

    def ancestor_iter(node, filter_fn=None, inclusive=True):
        """
        Returns all ancestors of node. If `inclusive` is True, `node`
        is returned as the first item of the sequence.
        """
        if inclusive:
            yield node
        while node is not None:
            node = node.parent_node
            if node is not None \
                   and (filter_fn is None or filter_fn(node)):
                yield node

    ancestor_iter = staticmethod(ancestor_iter)

    ## UTILITIES #############################################################

    def nodeset_hash(nodes, attribute='elem_id'):
        """
        Returns a hash of a set of nodes, based on the given
        attribute.
        """
        tags = []
        for node in nodes:
            if hasattr(node, attribute) and getattr(node, attribute) != None:
                value = getattr(node, attribute)
                tags.append(str(value))
        tags.sort()
        return '+'.join(tags)
    
    nodeset_hash = staticmethod(nodeset_hash)
    
    ## INSTANCE METHODS########################################################

    def __init__(self, elem_id=None, label=None, taxon=None, edge=None):
        """
        Inits. Handles keyword arguments: `elem_id` and `label`.
        """
        taxa.TaxonLinked.__init__(self, elem_id=elem_id, label=label)
        self.__edge = None        
        self.__child_nodes = []        
        self.__parent_node = None        
#         self.__next_sib = None
        if edge is not None:
            self.edge = edge
        else:
            self.edge = Edge(head_node=self)
        self.__edge.head_node = self            

    def __str__(self):
        """
        String representation of the object: it's id.
        """
        return str(self.elem_id)

    def is_leaf(self):
        "Returns True if the node has no children"
        return bool(not self.__child_nodes)
                
    ## Low-level methods for manipulating structure ##

    def _get_edge(self):
        """
        Returns the edge subtending this node.
        """
        return self.__edge

    def _set_edge(self, edge=None):
        """
        Sets the edge subtending this node, and sets head_node of
        `edge` to point to self.
        """
        self.__edge = edge
        if self.edge:
            self.edge.head_node = self

    edge = property(_get_edge, _set_edge)
        
    def children(self):
        """
        Returns the a shallow-copy list of all child nodes.
        """
        return list(self.__child_nodes)
    
    def set_children(self, child_nodes):
        """
        Sets the child_nodes for this node.
        Side effects: 
            - sets the parent of each child node to this node
            - sets the tail node of each child to self
#             - sets the next_sib of each child correctly
        """
        self.__child_nodes = child_nodes
        for nidx in range(len(self.__child_nodes)):
            self.__child_nodes[nidx].parent = self
            self.__child_nodes[nidx].edge.tail_node = self
#             if nidx < len(self.__child_nodes)-1:
#                 self.__child_nodes[nidx].next_sib = self.__child_nodes[nidx+1]
#             else:
#                 self.__child_nodes[nidx].next_sib = None
    
    def _get_parent_node(self):
        """Returns the parent node of this node."""
        return self.__parent_node
    
    def _set_parent_node(self, parent):
        """Sets the parent node of this node."""
        self.__parent_node = parent
        self.edge.tail_node = parent
        
    parent_node = property(_get_parent_node, _set_parent_node)

#     def _get_next_sib(self):
#         """Returns the next sibling of this node."""
#         return self.__next_sib
#     
#     def _set_next_sib(self, next_sib):
#         """Sets the next sibling of this node."""
#         self.__next_sib = next_sib
#         
#     next_sib = property(_get_next_sib, _set_next_sib)

    def add_child(self, node, edge_length=None):
        """
        Adds a child node to this node. Results in the parent_node and
        containing_tree of the node being attached set to this node.
        If `edge_length` is given, then the new child's edge length is
        set to this. Returns node that was just attached.
        """
        node.parent_node = self
        node.edge.tail_node = self
        if edge_length != None:
            node.edge.length = edge_length
#         if len(self.__child_nodes) > 0:
#             self.__child_nodes[-1].next_sib = node
        self.__child_nodes.append(node)
        return node

    def new_child(self, elem_id=None, edge_length=None, node_label=None, node_taxon=None):
        """
        Convenience class to create and add a new child to this node.
        """
        node = self.__class__()
        if elem_id is not None:
            node.elem_id = elem_id
        if node_label is not None:
            node.label = node_label
        if node_taxon is not None:
            node.taxon = node_taxon
        return self.add_child(node, edge_length)

    def remove_child(self, node):
        """
        Removes a node from this nodes child set. Results in the
        parent of the node being removed set to None. Returns node
        that was just removed.
        """
        if node and node in self.__child_nodes:
            node.parent_node = None
            node.edge.tail_node = None
            index = self.__child_nodes.index(node)
#             if index > 0:
#                 self.__child_nodes[index-1].next_sib = None
            self.__child_nodes.remove(node)
        else:
            raise Exception("Tried to remove an non-existing or null node")
        return node
        
    ## Basic node metrics ##

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
                par_node = self.parent_node
                # The root is identified when a node with no
                # parent is encountered. If we want to use some
                # other criteria (e.g., where a is_root property
                # is True), we modify it here.
                while par_node:
                    if par_node.edge.length != None:
                        distance_from_root = distance_from_root + float(par_node.edge.length)
                    par_node = par_node.parent_node
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

    def level(self):
        """
        Number of nodes between self and root.
        """
        if self.parent_node:
            return self.parent_node.level + 1
        else:
            return 0
    
    def distance_from_tip(self):
        """
        Sum of edge lengths from tip to node. If tree is not ultrametric
        (i.e., descendent edges have different lengths), then count the
        maximum of edge lengths.
        """
        if not self.__child_nodes:
            return 0.0
        else:
            distance_from_tips = []
            for ch in self.__child_nodes:
                if ch.edge.length is not None:
                    curr_edge_length = ch.edge.length
                else:
                    curr_edge_length = 0.0
                distance_from_tips.append(ch.distance_from_tip() + curr_edge_length)                    
            return float(max(distance_from_tips))

    def leaf_nodes(self):
        """
        Returns list of all leaves descended from this node (or just
        list with self as the only member if self is a leaf).
        """
        return [node for node in \
                self.postorder_iter(self, \
                                    lambda x: bool(len(node.child_nodes)==0))]

    def sib_nodes(self):
        """
        Returns all children of parent except self.
        """
        if self.parent_node is not None:
            return [node for node in self.parent_node.children() \
                    if node != self]
        else:
            return []

    def supratree_nodes(self):
        """
        Returns all nodes on a tree that do not descend from this edge.
        """
        nodes = []
        node = self.parent_node
        while node is not None:
            for sib in node.sib_nodes():
                nodes.extend(sib.infratree_nodes())
            if node.parent_node:
                nodes.append(node.parent_node)
            node = node.parent_node
        return nodes

    def infratree_nodes(self):
        """
        Returns self and all nodes descended from self.
        """
        return [node for node in self.preorder_iter(self)]

    #### BELOW TO BE MOVED INTO A SPLITS CLASS ####

    def local_split_set(self, attribute='elem_id'):
        """
        The split on the edge subtending this node is represented as
        a set with two members: an infratree hash and a supratree hash.
        """
        if self.parent_node and self.__child_nodes:
            infra = self.nodeset_hash(self.infratree_nodes(), attribute)
            supra = self.nodeset_hash(self.supratree_nodes(), attribute)
            return frozenset([infra, supra])
        else:
            return None

    def subtree_splits(self, attribute='elem_id'):
        """
        Returns list of all sets of splits on the subtree descending
        from this node, using `attribute` to compose the hash.
        """
        if self.__child_nodes:
            subsplits = set()
            if self.local_split_set(attribute) != None:
                subsplits.add(self.local_split_set(attribute))
            if self.__child_nodes:
                for child in self.__child_nodes:
                    child_splits = child.subtree_splits(attribute)
                    if child_splits:
                        subsplits.update(child_splits)
            return frozenset(subsplits)
        else:
            return frozenset()

    #### ABOVE TO BE MOVED INTO A SPLITS CLASS ####

    def new_node(self, elem_id=None, label=None):
        """
        Returns a new node object of the same class as this node.
        """
        node = self.__class__(elem_id)
        node.label = label
        return node

    def new_edge(self, elem_id=None):
        """
        Returns a new edge object of the class of this node's edge.
        """
        edge = self.edge.new_edge()
        edge.elem_id = elem_id
        return edge
     
    ### FOR DEBUGGING ### 
    def compose_newick(self, include_internal_labels=True):
        """
        This returns the Node as a NEWICK
        statement according to the given formatting rules.
        """
        statement = ''
        children = self.children()
        if children:
            subnodes = [child.compose_newick() for child in children]
            statement = '(' + ','.join(subnodes) + ')'
            
        if hasattr(self, 'taxon') and self.taxon:
            tag = self.taxon.label
        elif hasattr(self, 'label') and self.label and (len(children)==0 or include_internal_labels):
            tag = self.label
        elif len(children) == 0:
            tag = self.elem_id
        else:
            tag = ""
        if tag.count(' '):
            if not (tag.startswith("\'") and tag.endswith("\'")) \
               and not (tag.startswith("\"") and tag.endswith("\"")):
                tag = "'" + tag + "'"
        
        statement = statement + tag

        if self.edge and self.edge.length != None:
            try:
                statement =  "%s:%f" \
                            % (statement, float(self.edge.length))
            except ValueError:
                statement =  "%s:%s" \
                            % (statement, self.edge.length)
        return statement       

##############################################################################
## Edge

class Edge(base.IdTagged):
    """
    An edge on a tree. This class implements only the core
    functionality needed for trees.
    """

#     ## ITERATORS #############################################################
    
#     def preorder_iter(node, filter_fn=None):
#         """
#         Returns the edge of node and all descendents such that parents are
#         returned before children.
#         """
#         return Edge.__node_edge_iter(node, node.preorder_iter, filter_fn)

#     preorder_iter = staticmethod(preorder_iter)

#     def postorder_iter(node, filter_fn=None):
#         """
#         Returns the edge of node and all descendents such that parents are
#         returned before children.
#         """
#         return Edge.__node_edge_iter(node, node.postorder_iter, \
#                                           filter_fn)

#     postorder_iter = staticmethod(postorder_iter)

#     def __node_edge_iter(node, node_iter, filter_fn=None):
#         """
#         Uses given node_iter to iterate over nodes, returning edges.
#         """
#         node_filter = lambda x: bool(filter_fn is None or filter_fn(x.edge))
#         edge_caster = lambda x: x.edge
#         return utils.RecastingIterator(node_iter(node, node_filter), \
#                                        edge_caster)

#     __node_edge_iter = staticmethod(__node_edge_iter)

    ## CLASS METHODS  ########################################################
    
    def __init__(self,
                 elem_id=None,
                 head_node=None,
                 tail_node=None,
                 length=None):
        """
        Creates an edge from tail_node to head_node.  Modified from
        arbol.
        """
        base.IdTagged.__init__(self, elem_id=elem_id)
        self.__tail_node = None
        self.__head_node = None
        self.__tail_elem_id = None
        self.__head_elem_id = None
        self.rootedge = False

        self.elem_id = elem_id
        if head_node is not None:
            self._set_head_node(head_node)
        if tail_node is not None:
            self._set_tail_node(tail_node)
        elif self.head_node:
            self._set_tail_node(self.head_node.parent_node)
        self.length = length

    def _get_tail_node(self):
        """
        Returns the source node.
        """
        return self.__tail_node

    def _set_tail_node(self, node):
        """
        Sets the source node. Note: this *does* change the
        tail_elem_id as well!
        """
        self.__tail_node = node
        if self.__tail_node:
            self.__tail_elem_id = node.elem_id
        else:
            self.__tail_elem_id = None

    tail_node = property(_get_tail_node, _set_tail_node)

    def _get_tail_elem_id(self):
        """
        Returns the given id of the source node if defined.
        """
        if self.__tail_elem_id == None and self.__tail_node:
            self.__tail_elem_id = self.__tail_node.elem_id
        return self.__tail_elem_id

    def _set_tail_elem_id(self, elem_id):
        """
        Sets the source node id. Note: does *not* change the
        tail_node itself!
        """
        self.__tail_elem_id = elem_id

    tail_elem_id = property(_get_tail_elem_id, _set_tail_elem_id)
                        
    def _get_head_node(self):
        """
        Returns the target node.
        """
        return self.__head_node

    def _set_head_node(self, node):
        """
        Sets the source node. Note: this *does* change the
        head_elem_id as well!
        """
        self.__head_node = node
        if self.__head_node:
            self.__head_elem_id = node.elem_id
        else:
            self.__head_elem_id = None

    head_node = property(_get_head_node, _set_head_node)

    def _get_head_elem_id(self):
        """
        Returns the given id of the target node if defined.
        """
        if self.__head_elem_id == None and self.__head_node:
            self.__head_elem_id = self.__head_node.elem_id
        return self.__head_elem_id

    def _set_head_elem_id(self, elem_id):
        """
        Sets the target node id. Note: does *not* change the
        head_node itself!
        """
        self.__head_elem_id = elem_id

    head_elem_id = property(_get_head_elem_id, _set_head_elem_id)

    def bisect(self, bisecting_node):
        """
        Adds a new node, `node`, of the same class as the bisecting
        node as the ancestor of `self.head_node`, and attaches `node`
        to this node such that `node` is sister to `self.head_node`.
        !UNTESTED!
        """
        new_node = bisecting_node.__class__()
        new_node.add_child(self.head_node)
        new_node.add_child(bisecting_node)
        self.tail_node.remove_child(self.head_node)
        self.tail_node.add_child(new_node)
        self.head_node.parent_node = new_node
        bisecting_node.parent_node = new_node
        return new_node

    def new_edge(self, elem_id=None):
        """
        Returns a new edge object of the same class of this edge.
        """
        edge = self.__class__()
        edge.elem_id = elem_id
        return edge
