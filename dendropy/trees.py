#! /usr/bin/env python

############################################################################
##  trees.py
##
##  Part of the DendroPy library for phylogenetic computing.
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
import math

##############################################################################
## TreesBlock

class TreesBlock(list, taxa.TaxaLinked):
    """
    Tree manager.
    """

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `oid`, `label` and `taxa_block`.
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
                    node.taxon = taxa_block.get_taxon(label=node.taxon.label, update=True)
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
    root node as a node without child_nodes.
    """

    ###########################################################################
    ## Static methods
    
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

    ###########################################################################
    ## Special/Lifecycle methods
    
    def __init__(self, oid=None, label=None, seed_node=None):
        """
        Initializes a Tree object by defining a base node which must
        be of type `Node` or derived from `Node`.
        """
        base.IdTagged.__init__(self, oid=oid, label=label)
        self.seed_node = None
        self.length_type = None
        if seed_node is not None:
            self.seed_node = seed_node
        else:
            self.seed_node = Node(oid='n0', edge=Edge())
            
    def __str__(self):
        """
        Dump Newick string.
        """
        return self.compose_newick()

    ###########################################################################
    ## Getting/accessing methods

    def nodes(self, cmp_fn=None, filter_fn=None):
        """
        Returns list of nodes on the tree, sorted using cmp_fn.
        """
        nodes = [node for node in self.preorder_node_iter(filter_fn)]
        if cmp_fn:
            nodes.sort(cmp_fn)
        return nodes
    
    def leaf_nodes(self):
        """
        Returns list of leaf_nodes on the tree.
        """
        return [leaf for leaf in self.leaf_iter()]

    def find_node(self, filter_fn):
        """
        Finds the first node for which filter_fn(node) = True.
        For example, if::
        
            filter_fn = lambda n: hasattr(n, 'genes') and n.genes is not None

        then::
            
            t.find_node(filter_fn=filter_fn)
            
        will return all nodes which have an attributed 'genes' and this value
        is not None.
        """
        found = [node for node in self.preorder_node_iter(filter_fn)]
        if found and len(found) > 0:
            return found[0]
        else:
            return None
            
    def find_taxon_node(self, taxon_filter_fn=None):
        """
        Finds the first node for which taxon_filter_fn(node.taxon) == True.
        """
        for node in self.preorder_node_iter():
            if hasattr(node, "taxon") and node.taxon is not None:
                if taxon_filter_fn(node.taxon):
                    return node
        return None                    

    def find_edge(self, oid):
        """
        Finds the first edge with matching id.
        """
        filter_fn = lambda x: x.oid == oid
        found = [edge for edge in self.preorder_edge_iter(filter_fn)]
        if found and len(found) > 0:
            return found[0]
        else:
            return None

    ###########################################################################
    ## Node iterators

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
        Returns an iterator over tree leaf_nodes (order determined by
        postorder tree-traversal).
        """
        for node in self.seed_node.leaf_iter(self.seed_node, filter_fn):
            yield node

    ###########################################################################
    ## Edge iterators

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
    
    ###########################################################################
    ## Taxa 
    
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
                node.taxon = taxa_block.get_taxon(label=node.taxon.label, \
                    update=update_taxa_block)       
                    
    ###########################################################################
    ## Structure                    
                    
    def deroot(self):
        """
        Deroot the tree.
        """
        if self.seed_node:
            child_nodes = self.seed_node.child_nodes()
            if child_nodes and len(child_nodes) == 2:
                if len(child_nodes[0].child_nodes()) >= 2:
                    new_child = child_nodes[0]
                    new_seed = child_nodes[1]
                else: #lif len(child_nodes[1].child_nodes()) >= 2:
                    new_child = child_nodes[1]
                    new_seed = child_nodes[0]
                new_edge_length = 0.0
                if new_child.edge.length:
                    new_edge_length += new_child.edge.length
                if new_seed.edge.length:
                    new_edge_length += new_seed.edge.length
                    new_seed.edge = None
                self.seed_node = new_seed
                self.seed_node.add_child(new_child)                    
                
    ###########################################################################
    ## For debugging
    
    def compose_newick(self, include_internal_labels=True):
        return self.seed_node.compose_newick(include_internal_labels=include_internal_labels)
                


        
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
        Preorder traversal of the node and its child_nodes.  Returns node
        and all descendants such that node is returned before node's
        child_nodes (and their child_nodes). Filtered by filter_fn: node is
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
            child_nodes = node.child_nodes()
            child_nodes.extend(stack)
            stack = child_nodes
    preorder_iter = staticmethod(preorder_iter)

    def postorder_iter(node, filter_fn=None):
        """
        Postorder traversal of the node and its child_nodes.  Returns node
        and all descendants such that node's child_nodes (and their
        child_nodes) are visited before node.  Filtered by filter_fn:
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
                child_nodes = [(n, False) for n in node.child_nodes()]
                child_nodes.extend(stack)
                stack = child_nodes
    postorder_iter = staticmethod(postorder_iter)

    def leaf_iter(start_nd, filter_fn=None):
        """
        Returns an iterator over the leaf_nodes that are descendants `of start_nd`
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
        Level-order traversal of the node and its child_nodes. Filtered
        by filter_fn: node is only returned if no filter_fn is given
        or if filter_fn returns True
        """
        if filter_fn is None or filter_fn(node):
            yield node
        remaining = node.child_nodes()
        while len(remaining) > 0:
            node = remaining.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = node.child_nodes()
            remaining.extend(child_nodes)
            
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

    def nodeset_hash(nodes, attribute='oid'):
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

    def __init__(self, oid=None, label=None, taxon=None, edge=None):
        """
        Inits. Handles keyword arguments: `oid` and `label`.
        """
        taxa.TaxonLinked.__init__(self, oid=oid, label=label)
        self.__edge = None        
        self.__child_nodes = []        
        self.__parent_node = None        
        if edge is not None:
            self.edge = edge
        else:
            self.edge = Edge(head_node=self)
        self.__edge.head_node = self            

    def __str__(self):
        """
        String representation of the object: it's id.
        """
        return str(self.oid)

    def is_leaf(self):
        "Returns True if the node has no child_nodes"
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
        if edge:
            edge.head_node = self

    def _get_edge_length(self):
        """
        Returns the length of the edge  subtending this node.
        """
        return self.__edge.length

    def _set_edge_length(self, v=None):
        """
        Sets the edge subtending this node, and sets head_node of
        `edge` to point to self.
        """
        self.__edge.length = v
        

    edge = property(_get_edge, _set_edge)
    edge_length = property(_get_edge_length, _set_edge_length)
        
    def child_nodes(self):
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
        """
        self.__child_nodes = child_nodes
        for nidx in range(len(self.__child_nodes)):
            self.__child_nodes[nidx].parent = self
            self.__child_nodes[nidx].edge.tail_node = self

    def _get_parent_node(self):
        """Returns the parent node of this node."""
        return self.__parent_node
    
    def _set_parent_node(self, parent):
        """Sets the parent node of this node."""
        self.__parent_node = parent
        self.edge.tail_node = parent
        
    parent_node = property(_get_parent_node, _set_parent_node)

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
            node.edge_length = edge_length
#         if len(self.__child_nodes) > 0:
#             self.__child_nodes[-1].next_sib = node
        self.__child_nodes.append(node)
        return node

    def new_child(self, oid=None, edge_length=None, node_label=None, node_taxon=None):
        """
        Convenience class to create and add a new child to this node.
        """
        node = self.__class__()
        if oid is not None:
            node.oid = oid
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
                    curr_edge_length = ch.edge_length
                else:
                    curr_edge_length = 0.0
                distance_from_tips.append(ch.distance_from_tip() + curr_edge_length)                    
            return float(max(distance_from_tips))

    def leaf_nodes(self):
        """
        Returns list of all leaf_nodes descended from this node (or just
        list with self as the only member if self is a leaf).
        """
        return [node for node in \
                self.postorder_iter(self, \
                                    lambda x: bool(len(node.child_nodes)==0))]
     
    ########################################################################### 
    ## for debugging
    
    def compose_newick(self, include_internal_labels=True):
        """
        This returns the Node as a NEWICK
        statement according to the given formatting rules.
        """
        statement = ''
        child_nodes = self.child_nodes()
        if child_nodes:
            subnodes = [child.compose_newick() for child in child_nodes]
            statement = '(' + ','.join(subnodes) + ')'
            
        if hasattr(self, 'taxon') and self.taxon:
            tag = self.taxon.label
        elif hasattr(self, 'label') and self.label and (len(child_nodes)==0 or include_internal_labels):
            tag = self.label
        elif len(child_nodes) == 0:
            tag = self.oid
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

    ## CLASS METHODS  ########################################################
    
    def __init__(self,
                 oid=None,
                 head_node=None,
                 tail_node=None,
                 length=None):
        """
        Creates an edge from tail_node to head_node.  Modified from
        arbol.
        """
        base.IdTagged.__init__(self, oid=oid)
        self.tail_node = None
        self.head_node = None
        self.rootedge = False
        if head_node is not None:
            self.head_node = head_node
        if tail_node is not None:
            self.tail_node = tail_node
        self.length = length

    def new_edge(self, oid=None):
        """
        Returns a new edge object of the same class of this edge.
        """
        edge = self.__class__()
        edge.oid = oid
        return edge


