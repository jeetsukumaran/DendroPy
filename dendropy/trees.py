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
from cStringIO import StringIO

from dendropy import base
from dendropy import taxa
import math
from dendropy import get_logger
_LOG = get_logger('dendropy.trees')

##############################################################################
## Utilities

def add_depth_to_nodes(tree, prec=0.0000001):
    """Takes an ultrametric `tree` and adds a `depth` attribute to each internal
    node.  The `depth` is the sum of edge lengths from the node to the tips.
    
    If the lengths of different paths to the node differ by more than `prec`, 
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    """
    
    node = None    
    for node in tree.postorder_node_iter():
        ch = node.child_nodes()
        if len(ch) == 0:
            node.depth = 0.0
        else:
            first_child = ch[0]
            node.depth = first_child.depth + first_child.edge.length
            last_child = ch[-1]
            for nnd in ch[1:]:
                ocnd = nnd.depth + nnd.edge.length
                if abs(node.depth - ocnd) > prec:
                    raise ValueError("Tree is not ultrametric")
    if node is None:
        raise ValueError("Empty tree encountered")

##############################################################################
## TreesBlock

class TreesBlock(list, taxa.TaxaLinked):
    "Tree manager."

    def __init__(self, *args, **kwargs):
        "Inits. Handles keyword arguments: `oid`, `label` and `taxa_block`."
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
            tree.normalize_taxa(taxa_block)
        taxa_block.sort()
        self.taxa_block = taxa_block
        return taxa_block        

    def __setitem__(self, key, tree):
        """
        Makes sure tree.taxa_block = self.taxa_block.
        """
        if tree.taxa_block is not self.taxa_block:
            tb_mutable = self.taxa_block._is_mutable
            self.taxa_block._is_mutable = True
            tree.normalize_taxa(self.taxa_block)
            self.taxa_block._is_mutable = tb_mutable
        list.__setitem__(self, key, tree)
        
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
    
    def ancestor(node1, node2):
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

    ancestor = staticmethod(ancestor)
    
    ###########################################################################
    ## Special/Lifecycle methods
    
    def __init__(self, oid=None, label=None, seed_node=None, taxa=None):
        """
        Initializes a Tree object by defining a base node which must
        be of type `Node` or derived from `Node`.
        """
        base.IdTagged.__init__(self, oid=oid, label=label)
        self.seed_node = None
        self.length_type = None
        self.is_rooted = False
        if seed_node is not None:
            self.seed_node = seed_node
        else:
            self.seed_node = Node(oid='n0', edge=Edge())
        self.taxa_block = taxa
    def __str__(self):
        "Dump Newick string."
        return self.compose_newick()

    ###########################################################################
    ## Getting/accessing methods

    def nodes(self, cmp_fn=None, filter_fn=None):
        "Returns list of nodes on the tree, sorted using cmp_fn."
        nodes = [node for node in self.preorder_node_iter(filter_fn)]
        if cmp_fn:
            nodes.sort(cmp_fn)
        return nodes
    
    def leaf_nodes(self):
        "Returns list of leaf_nodes on the tree."
        return [leaf for leaf in self.leaf_iter()]

    def internal_nodes(self):
        "Returns list of internal node in the tree."
        return self.nodes(filter_fn=lambda x : not x.is_leaf())
    
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
        for node in self.preorder_node_iter(filter_fn):
            return node
        return None
            
    def find_taxon_node(self, taxon_filter_fn=None):
        "Finds the first node for which taxon_filter_fn(node.taxon) == True."
        for node in self.preorder_node_iter():
            if hasattr(node, "taxon") and node.taxon is not None:
                if taxon_filter_fn(node.taxon):
                    return node
        return None                    

    def find_edge(self, oid):
        "Finds the first edge with matching id."
        for edge in self.preorder_edge_iter():
            if edge.oid == oid:
                return edge
        return None

    ###########################################################################
    ## Node iterators

    def preorder_node_iter(self, filter_fn=None):
        "Returns preorder iterator over tree nodes."
        for node in self.seed_node.preorder_iter(self.seed_node, filter_fn):
            yield node

    def postorder_node_iter(self, filter_fn=None):
        "Returns postorder iterator over tree nodes."
        for node in self.seed_node.postorder_iter(self.seed_node, filter_fn):
            yield node

    def level_order_node_iter(self, filter_fn=None):
        "Returns level-order iterator over tree nodes."
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
        "Returns preorder iterator over tree edges."
        for node in self.seed_node.preorder_iter(self.seed_node):
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def postorder_edge_iter(self, filter_fn=None):
        "Returns postorder iterator over tree edges."
        for node in self.seed_node.postorder_iter(self.seed_node):
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def level_order_edge_iter(self, filter_fn=None):
        "Returns level-order iterator over tree edges."
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
        self.taxa_block = taxa_block
        return taxa_block
        
    def normalize_taxa(self, taxa_block):
        """
        Reassigns tree taxa objects to corresponding taxa objects in
        given taxa_block, with identity of taxa objects determined by
        labels.
        """
        for node in self.postorder_node_iter():
            t = node.taxon
            if t:
                node.taxon = taxa_block.get_taxon(label=t.label)
        self.taxa_block = taxa_block
                    
    ###########################################################################
    ## Structure                    
                    
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
        from dendropy.treemanip import collapse_edge
        collapse_edge(to_del_edge)
                
    ###########################################################################
    ## For debugging
    
    def compose_newick(self, **kwargs):
        """kwargs["reverse_translate"] can be function that takes a taxon and 
           returns the label to appear in the tree."""
        return self.seed_node.compose_newick(**kwargs)
                
    def reroot_at(self, nd, splits=False, delete_deg_two=True):
        """Takes an internal node, `nd` that must already be in the tree and 
        reroots the tree such that `nd` is the `seed_node` of the tree.
        
        If `splits` is True, then the edges' `clade_mask` and the tree's 
            `split_edges` attributes will be updated."""
        old_par = nd.parent_node
        if old_par is None:
            return
        if splits:
            taxa_mask = self.seed_node.edge.clade_mask
        to_edge_dict = None
        if splits:
            to_edge_dict = getattr(self, "split_edges", None)

        if old_par is self.seed_node:
            root_children = old_par.child_nodes()
            if len(root_children) == 2 and delete_deg_two:
                # root (old_par) was of degree 2, thus we need to suppress the
                #   node
                fc = root_children[0]
                if fc is nd:
                    sister = root_children[1]
                else:
                    assert root_children[1] is nd
                    sister = fc
                if nd.edge.length:
                    sister.edge.length += nd.edge.length
                edge_to_del = nd.edge
                nd.edge = old_par.edge
                if splits:
                    assert nd.edge.clade_mask == taxa_mask
                if to_edge_dict:
                    del to_edge_dict[edge_to_del.clade_mask]
                nd.add_child(sister, edge_length=sister.edge.length)
                self.seed_node = nd
                return
        else:
            self.reroot_at(old_par, splits=splits, delete_deg_two=delete_deg_two)
        old_par.edge, nd.edge = nd.edge, old_par.edge
        e = old_par.edge
        if splits:
            if to_edge_dict:
                del to_edge_dict[e.clade_mask]
            e.clade_mask = (~(e.clade_mask)) & taxa_mask
            if to_edge_dict:
                to_edge_dict[e.clade_mask] = e
            assert nd.edge.clade_mask == taxa_mask
        old_par.remove_child(nd)
        nd.add_child(old_par, edge_length=e.length)
        self.seed_node = nd

    def to_outgroup_position(self, nd, splits=False, delete_deg_two=True):
        """Reroots the tree at the parend of `nd` and makes `nd` the first child
        of the new root.  This is just a convenience function to make it easy
        to place a clade as the first child under the root.
        
        Assumes that `nd` and `nd.parent_node` and are in the tree 
        
        If `splits` is True, then the edges' `clade_mask` and the tree's 
            `split_edges` attributes will be updated.
        If `delete_deg_two` is True and the old root of the tree has an 
            outdegree of 2, then the node will be removed from the tree.
        """
        p = nd.parent_node
        assert p is not None
        self.reroot_at(p, splits=splits)
        p.remove_child(nd)
        p.add_child(nd, edge_length=nd.edge.length, pos=0)
        
    def get_indented_form(self, **kwargs):
        out = StringIO()
        self.write_indented_form(out, **kwargs)
        return out.getvalue()

    def write_indented_form(self, out, **kwargs):
        if kwargs.get("splits"):
            if not kwargs.get("taxa"):
                kwargs["taxa"] = self.taxa_block
        self.seed_node.write_indented_form(out, **kwargs)

    def debug_check_tree(self, logger_obj=None, **kwargs):
        import logging, inspect
        if logger_obj and logger_obj.isEnabledFor(logging.DEBUG):
            try:
                assert self._debug_tree_is_valid(logger_obj=logger_obj, **kwargs)
            except:
                calling_frame = inspect.currentframe().f_back
                co = calling_frame.f_code
                emsg = "\nCalled from file %s, line %d, in %s" % (co.co_filename, calling_frame.f_lineno, co.co_name)
                _LOG.debug("%s" % str(self))
                _LOG.debug("%s" % self.get_indented_form(**kwargs))
        assert self._debug_tree_is_valid(logger_obj=logger_obj, **kwargs)

    def _debug_tree_is_valid(self, **kwargs):
        """Performs sanity-checks of the tree data structure.
        
        kwargs:
            `splits` if True specifies that the split_edge and clade_mask attributes
                are checked.
        """
        check_splits = kwargs.get('splits', False)
        taxa_block = kwargs.get('taxa_block')
        if taxa_block is None:
            taxa_block = self.taxa_block
        if check_splits:
            taxa_mask = self.seed_node.edge.clade_mask
        nodes = set()
        edges = set()
        curr_node = self.seed_node
        assert(curr_node.parent_node is None)
        assert(curr_node.edge.tail_node is None)
        ancestors = []
        siblings = []
        while curr_node:
            curr_edge = curr_node.edge
            assert(curr_edge not in edges)
            edges.add(curr_edge)
            assert(curr_node not in nodes)
            nodes.add(curr_node)
            assert(curr_edge.tail_node is curr_node.parent_node)
            assert(curr_edge.head_node is curr_node)
            if check_splits:
                cm = 0
                clade_mask = curr_edge.clade_mask
                assert((clade_mask | taxa_mask) == taxa_mask)
            c = curr_node.child_nodes()
            if c:
                for child in c:
                    assert child.parent_node is curr_node
                    if check_splits:
                        cm |= child.edge.clade_mask
            elif check_splits:
                cm = taxa_block.taxon_bitmask(curr_node.taxon)
            if check_splits:
                assert((cm & taxa_mask) == clade_mask)
                assert self.split_edges[clade_mask] == curr_edge
            curr_node, level = _preorder_list_manip(curr_node, siblings, ancestors)
        if check_splits:
            for s, e in self.split_edges.iteritems():
                assert(e in edges)
        return True

       
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
        "Inits. Handles keyword arguments: `oid` and `label`."
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
        "String representation of the object: it's id."
        return str(self.oid)

    def is_leaf(self):
        "Returns True if the node has no child_nodes"
        return bool(not self.__child_nodes)
                
    ## Low-level methods for manipulating structure ##

    def _get_edge(self):
        "Returns the edge subtending this node."
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
        "Returns the length of the edge  subtending this node."
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
        "Returns the a shallow-copy list of all child nodes."
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
            self.__child_nodes[nidx].parent_node = self
            self.__child_nodes[nidx].edge.tail_node = self

    def _get_parent_node(self):
        """Returns the parent node of this node."""
        return self.__parent_node
    
    def _set_parent_node(self, parent):
        """Sets the parent node of this node."""
        self.__parent_node = parent
        self.edge.tail_node = parent
        
    parent_node = property(_get_parent_node, _set_parent_node)

    def add_child(self, node, edge_length=None, pos=None):
        """
        Adds a child node to this node. Results in the parent_node and
        containing_tree of the node being attached set to this node.
        If `edge_length` is given, then the new child's edge length is
        set to this. Returns node that was just attached.
        """
        node.parent_node = self
        if edge_length != None:
            node.edge_length = edge_length
#         if len(self.__child_nodes) > 0:
#             self.__child_nodes[-1].next_sib = node
        if pos is None:
            self.__child_nodes.append(node)
        else:
            self.__child_nodes.insert(pos, node)
        return node

    def new_child(self, oid=None, edge_length=None, node_label=None, node_taxon=None):
        "Convenience class to create and add a new child to this node."
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
        if not node:
            raise Exception("Tried to remove an non-existing or null node")
        if node in self.__child_nodes:
            node.parent_node = None
            node.edge.tail_node = None
            index = self.__child_nodes.index(node)
#             if index > 0:
#                 self.__child_nodes[index-1].next_sib = None
            self.__child_nodes.remove(node)
        else:
            raise Exception("Tried to remove a node that is not listed as a child")
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

    def level(self):
        "Number of nodes between self and root."
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


    def get_node_str(self, **kwargs):
        """returns a string that is an identifier for the node.  This is called
        by the newick-writing functions, so the kwargs that affect how node 
        labels show up in a newick string are the same ones used here:
            `include_internal_labels` is a Boolean
        """
        is_leaf = (len(self.__child_nodes) == 0)
        include_internal_labels = kwargs.get("include_internal_labels")
        if (not is_leaf) and (not include_internal_labels):
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
                    tag = self.oid
        if "raw_labels" in kwargs:
            return tag
        from dendropy.nexus import NexusWriter 
        return NexusWriter.escape_token(tag)

    ########################################################################### 
    ## for debugging
    
    def compose_newick(self, **kwargs):
        """
        This returns the Node as a NEWICK
        statement according to the given formatting rules.
        """
        out = StringIO()
        self.write_newick(out, **kwargs)
        return out.getvalue()

    def write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK
        statement according to the given formatting rules.
        """
        edge_lengths = kwargs.get('edge_lengths', True)
        child_nodes = self.child_nodes()
        if child_nodes:
            out.write('(')
            f_child = child_nodes[0]
            for child in child_nodes:
                if child is not f_child:
                    out.write(',')
                child.write_newick(out, **kwargs)
            out.write(')')

        out.write(self.get_node_str(**kwargs))
        e = self.edge
        if e:
            sel = e.length
            if sel is not None:
                s = ""
                try:
                    s = float(sel)
                    s = str(s)
                except ValueError:
                    s = str(sel)
                if s:
                    out.write(":%s" % s)

    def get_indented_form(self, **kwargs):
        out = StringIO()
        self.write_indented_form(out, **kwargs)
        return out.getvalue()

    def write_indented_form(self, out, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        clade_masks = kwargs.get("splits", True)
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
            cm = "%s " % format_split(self.edge.clade_mask, **kwargs)
        else:
            cm = ""
        out.write("%s%s%s\n" % ( cm, indentation*level, label))
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
        "Returns a new edge object of the same class of this edge."
        edge = self.__class__()
        edge.oid = oid
        return edge
    def invert(self):
        self.head_node, self.tail_node = self.tail_node, self.head_node
    def is_terminal(self):
        "Returns True if the head node has no children"
        return bool(self.head_node and self.head_node.is_leaf())
    def is_internal(self):
        "Returns True if the head node has children"
        return bool(self.head_node and not self.head_node.is_leaf())

def _preorder_list_manip(n, siblings, ancestors):
    """Helper function for recursion free preorder traversal, that does not 
    rely on attributes of the node other than child_nodes() (thus it is useful
    for debuggging).
    
    Returns the next node (or None) and the number of levels toward the root 
    the function "moved"
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

def format_node(nd, **kwargs):
    if nd.is_leaf():
        t = nd.taxon
        if t:
            label = t.label
        else:
            label = "anonymous leaf"
    else:
        label = "* %s" % str(nd.oid)
    return label


def format_split(split, width=None, **kwargs):
    from dendropy.splits import split_as_string
    if width is None:
        width = len(kwargs.get("taxa"))
    s = split_as_string(split, width, symbol1=kwargs.get("off_symbol"), symbol2=kwargs.get("on_symbol"))
    return s

