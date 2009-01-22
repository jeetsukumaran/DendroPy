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
        if self.edge:
            self.edge.head_node = self

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

def pybus_harvey_gamma(tree, prec=0.00001):
    """Returns the gamma statistic of Pybus and Harvey (2000). This statistic 
    is used to test for constancy of birth and death rates over the course of
    a phylogeny.  Under the pure-birth process, the statistic should follow
    a standard Normal distibution: a Normal(mean=0, variance=1).
    
    If the lengths of different paths to the node differ by more than `prec`, 
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    Raises a Value Error if the tree is not ultrametric, is non-binary, or has
        only 2 leaves.
    
    As a side effect a `depth` attribute is added to the nodes of the tree.
    
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
    speciation_depths = []
    n = 0
    for node in tree.postorder_node_iter():
        ch = node.child_nodes()
        n_ch = len(ch)
        if n_ch == 0:
            node.depth = 0.0
            n += 1
        elif n_ch > 2:
            raise ValueError("Polytomy encountered")
        else:
            first_child = ch[0]
            node.depth = first_child.depth + first_child.edge.length
            last_child = ch[-1]
            for nnd in ch[1:]:
                ocnd = nnd.depth + nnd.edge.length
                if abs(node.depth - ocnd) > prec:
                    raise ValueError("Tree is not ultrametric")
            if n_ch == 2:
                speciation_depths.append(node.depth)
    if node is None:
        raise ValueError("Empty tree encountered")
    speciation_depths.sort(reverse=True)
    g = []
    older = speciation_depths[0]
    for age in speciation_depths[1:]:
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



def _calc_TKP_rate(starting_rate, duration, roeotroe, rng):
    """Returns a simulated rate for the head node of a tree when:
        * the tail node has rate `starting_rate`
        * the time duration of the edge is `duration`
        * the rate of evolution of the rate of evolution is `roeotroe` (this is
            the parameter nu in Kishino, Thorne, and Bruno 2001)
    `rng` is a random number generator.
    
    The model used to generate the rate is the one described by Thorne, Kishino,
    and Painter 1998.  The descendant rates or lognormally distributed.
    The mean rate returned will have an expectation of `starting_rate`
    The variance of the normal distribution for the logarithm of the ending rate 
        is the product of `duration` and `roeotroe`
    """
    rate_var = duration*roeotroe
    if rate_var > 0.0:
        mu = math.log(starting_rate)
        return rng.lognormvariate(mu, math.sqrt(rate_var))
    return starting_rate
    
def _calc_KTB_rate(starting_rate, duration, roeotroe, rng):
    """Returns a simulated rate for the head node of a tree when:
        * the tail node has rate `starting_rate`
        * the time duration of the edge is `duration`
        * the rate of evolution of the rate of evolution is `roeotroe` (this is
            the parameter nu in Kishino, Thorne, and Bruno 2001)
    `rng` is a random number generator.
    
    The model used to generate the rate is the one described by Kishino, Thorne,
    and Bruno 2001.  The descendant rates or lognormally distributed.
    The mean rate returned will have an expectation of `starting_rate`
    The variance of the normal distribution for the logarithm of the ending rate 
        is the product of `duration` and `roeotroe`
    """
    if starting_rate <= 0.0:
        raise ValueError("starting_rate must be positive in the KTB model")
    rate_var = duration*roeotroe
    if rate_var > 0.0:
        # Kishino, Thorne and Bruno corrected the tendency for the rate to
        #   increase seen in teh TKP, 1998 model
        mu = math.log(starting_rate) - (rate_var/2.0)
        return rng.lognormvariate(mu, math.sqrt(rate_var))
    return starting_rate

def _calc_KTB_rates_crop(starting_rate, duration, roeotroe, rng,  min_rate=None, max_rate=None):
    """Returns a descendant rate and mean rate according to the Kishino, Thorne,
    Bruno model.  Assumes that the min_rate <= starting_rate <= max_rate if a max 
    and min are provided.
    rate is kept within in the [min_rate, max_rate] range by cropping at these 
    values and acting is if the cropping occurred at an appropriate point
    in the duration of the branch (based on a linear change in rate from the
    beginning of the random_variate drawn for the end of the branch).
    """
    if roeotroe*duration <= 0.0:
        if (min_rate and starting_rate < min_rate) or (max_rate and starting_rate > max_rate):
            raise ValueError("Parent rate is out of bounds, but no rate change is possible")
    r = _calc_KTB_rate(starting_rate, duration, roeotroe, rng)
    if max_rate and r > max_rate:
        assert(starting_rate <= max_rate)
        p_changing =  (max_rate - starting_rate)/(r - starting_rate)
        mean_changing = (starting_rate + max_rate)/2.0
        mr = p_changing*mean_changing + (1.0 - p_changing)*max_rate
        return max_rate, mr
    elif min_rate and r < min_rate:
        assert(starting_rate >= min_rate)
        p_changing = (starting_rate - min_rate)/(starting_rate - r)
        mean_changing = (starting_rate + min_rate)/2.0
        mr = p_changing*mean_changing + (1.0 - p_changing)*min_rate
        return min_rate, mr
    return r, (starting_rate + r)/2.0

def _calc_KTB_rates_linear_bounce(starting_rate, duration, roeotroe, rng,  min_rate=0.0, max_rate=None):
    """Returns a descendant rate and mean rate according to the Kishino, Thorne,
    Bruno model.  Assumes that the min_rate <= starting_rate <= max_rate if a max 
    and min are provided.
    The rate is kept within in the [min_rate, max_rate] range by "bouncing" off
    of the barriers, where the "collision" is estimated by assuming a linear
    change in rate from the beginning of the random_variate drawn for the end 
    of the branch).
    """
    if roeotroe*duration <= 0.0:
        if (min_rate and starting_rate < min_rate) or (max_rate and starting_rate > max_rate):
            raise ValueError("Parent rate is out of bounds, but no rate change is possible")
    r = _calc_KTB_rate(starting_rate, duration, roeotroe, rng)
    if min_rate is None:
        min_rate = 0.0
    return bounce_constrain(starting_rate, r, min_rate, max_rate)

def bounce_constrain(start_x, x, min_x=None, max_x=None):
    """Returns the value of variable and its mean value over a path.
    We assume that some variable started at `start_x` and moved toward `x`, but 
    has to bounce of barriers specified by `min_x` and `max_x`.
    
    `x` determines the direction and magnitude of the change.
    
    `start_x` must fall in the legal range (between the min and max). If
    `x` is also legal, then (x, (x + start_x)/2.0) will be returned reflecting
    the fact that the arithmetic mean of the endpoints represents the mean value
    of the variable if it took a direct path (at constant rate).
    """
    
    if max_x is not None and min_x is not None:
        assert(max_x > min_x)
    gt_max = (max_x is not None  and x > max_x)
    lt_min = (min_x is not None and x < min_x)
    prev_x = start_x
    prop_dur_remaining = 1.0
    mx = 0.0
    while gt_max or lt_min:
        if gt_max:
            p_changing = (max_x - prev_x)/(x - prev_x)
            mean_changing = (prev_x + max_x)/2.0
            mx += p_changing*prop_dur_remaining*mean_changing
            prop_dur_remaining *= 1.0 - p_changing
            x = 2*max_x - x
            lt_min = (min_x is not None and x < min_x)
            prev_x =  max_x
            gt_max = False
        if lt_min:
            p_changing = (prev_x - min_x)/(prev_x - x)
            mean_changing = (prev_x + min_x)/2.0
            mx += prop_dur_remaining*p_changing*mean_changing
            prop_dur_remaining *= 1.0 - p_changing
            x = 2*min_x - x
            lt_min = False
            gt_max = (max_x is not None  and x > max_x)
            prev_x =  min_x
    mean_changing = (prev_x + x)/2.0
    mx += mean_changing*prop_dur_remaining
    return x, mx



def simulate_mutation_rates(node, rng, **kwargs):
    """Takes a node and a random number generator object, `rng` This function
    "evolves" a set of rates on the subtree descending from the  `node`.
    
    kwargs keys that are used are:
        `roeotroe` --  the rate of evolution of the rate of evolution. This 
            parameter that controls the degree of deviations away from the 
            molecular clock.  
        `min_rate` is the minimum rate (default None)
        `max_rate' is the maximum rate (default None), 
        `model` is a string specifying the name of the model. Currently only the
            KTB (Kishino, Thorne, Bruno) is supported
        `time_attr` is a string that specifies the name of the attribute 
            that returns the branch length in terms of time for a node. The 
            default is "edge_length"
        `rate_attr` is the string that specifies the name of the attribute 
            used to hold the rate for the nodes.  The root of the subtree is
            assumed to have this field on calling of the function.  On success
            all nodes in the subtree will have the attribute.  The default is
            "mutation_rate"
        `mean_rate_attr` if specified this is string that gives the name of 
            attribute in each node that is mean rate for the branch (default is
            None). This is filled in after time_attr and mut_rate_attr are read,
            so it is permissible to have this attribute match one of thos 
            strings (although it will make the model odd if the mean_rate_attr
            is the same as the rate_attr
         `constrain_rate_mode` controls the behavior when the minimum or maximum 
            rate is simulated. The choices are "crop", and "linear_bounce"
            "crop" means that the rate is set to the most extreme value allowed.
            "linear_bounce" considers the path of evolution of rate to be a
                simple line from the ancestor's rate to the proposed rate. The
                point at which the path crosses the extreme value is determined
                and the rate is "reflected off" the limiting rate at that point.
                This causes the rate to avoid the extreme values more than a 
                simulation of small time slices that simply rejects illegal 
                rates.
    
    Currently the only model supported is the one of Kishino, Thorne, and Bruno. 
    "Performance of a Divergence Time Estimation Method under a Probabilistic 
    Model of Rate Evolution." Molecular Biology and Evolution (2001) vol. 18 
    (3) pp. 352-361. This model is specified by the code "KTB". A node's rate
    is a log-normal variate with variance determined by the product of the 
    duration of the branch and the roeotroe parameter.  The mean of the 
    distribution is chosen such that mean of the log-normal is identical to the 
    rate at the parent. The mean_rate for the branch is the average of the rates
    at the endpoints.
    
    """
    nd_iter = Node.preorder_iter(node)
    # skip the first node -- it should already have a rate
    nd_iter.next()
    if kwargs.get("model", "KTB").upper() != "KTB":
        raise ValueError("Only the Kishino-Thorne-Bruno model is supported at this time")
    rate_attr = kwargs.get("rate_attr", "mutation_rate")
    if not rate_attr:
        raise ValueError("rate_attr cannot be an empty string")
    time_attr = kwargs.get("time_attr", "edge_length")
    mean_rate_attr = kwargs.get("time_attr")
    constrain_rate_mode = kwargs.get("constrain_rate_mode", "crop").lower()
    if constrain_rate_mode not in ["crop", "linear_bounce"]:
        raise ValueError('Only "crop" and "linear_bounce" are supported at this time')
    roeotroe = kwargs.get("roeotroe", 1.0)
    min_rate = kwargs.get("min_rate", 0.0)
    if min_rate < 0.0:
        raise ValueError("min_rate cannot be less than 0")
    max_rate = kwargs.get("max_rate")
    anc_rate = getattr(node, rate_attr)
    if max_rate is not None:
        if min_rate is not None:
            if min_rate > max_rate:
                raise ValueError("max_rate must be greater than the min_rate")
            if min_rate == max_rate:
                for nd in nd_iter:
                    setattr(nd, rate_attr, min_rate)
                    if mean_rate_attr:
                        # here we assume that the rate changed from the 
                        #   ancestral rate to the only allowed rate 
                        #   instantaneously, so the mean rat is min_rate
                        setattr(nd, mean_rate_attr, min_rate) 
                return
        if max_rate <= 0.0:
            raise ValueError("max_rate must be positive")
        if anc_rate > max_rate:
            raise ValueError("rate for the incoming node is > max_rate")
    if (min_rate is not None) and anc_rate < min_rate:
        raise ValueError("rate for the incoming node is > max_rate")
            
    if constrain_rate_mode == "crop":
        rate_func = _calc_KTB_rates_crop
    else:
        rate_func = _calc_KTB_rates_linear_bounce
    for nd in nd_iter:
        starting_rate = getattr(nd.parent_node, rate_attr)
        duration = getattr(nd, time_attr)
        r, mr  = rate_func(starting_rate, duration, roeotroe, rng, min_rate, max_rate)
        setattr(nd, rate_attr, r)
        if mean_rate_attr:
            setattr(nd, mean_rate_attr, 0.5*(min_rate + r))
