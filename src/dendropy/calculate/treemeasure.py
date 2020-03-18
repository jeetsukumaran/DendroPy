#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Statistics, metrics, measurements, and values calculated on (single) trees.
"""

import math
from dendropy.calculate import phylogeneticdistance

EULERS_CONSTANT = 0.5772156649015328606065120900824024310421

## legacy: will soon be deprecated
class PatristicDistanceMatrix(phylogeneticdistance.PhylogeneticDistanceMatrix):

    def __init__(self, tree):
        phylogeneticdistance.PhylogeneticDistanceMatrix.__init__(self)
        self.compile_from_tree(tree=tree)

def patristic_distance(tree, taxon1, taxon2, is_bipartitions_updated=False):
    """
    Given a tree with bipartitions encoded, and two taxa on that tree, returns the
    patristic distance between the two. Much more inefficient than constructing
    a PhylogeneticDistanceMatrix object.
    """
    mrca = tree.mrca(taxa=[taxon1, taxon2], is_bipartitions_updated=is_bipartitions_updated)
    dist = 0
    n = tree.find_node(lambda x: x.taxon == taxon1)
    while n != mrca:
        if n.edge.length is not None:
            dist += n.edge.length
        n = n.parent_node
    n = tree.find_node(lambda x: x.taxon == taxon2)
    while n != mrca:
        if n.edge.length is not None:
            dist += n.edge.length
        n = n.parent_node
    return dist

###########################################################################
### Metrics -- Unary

def B1(tree):
    """
    Returns the B1 statistic: the reciprocal of the sum of the maximum
    number of nodes between each interior node and tip over all internal
    nodes excluding root.
    """
    b1 = 0.0
    nd_mi = {}
    for nd in tree.postorder_node_iter():
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

def colless_tree_imbalance(tree, normalize="max"):
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
    for nd in tree.postorder_node_iter():
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
        raise TypeError("``normalization`` accepts only None, True, False, 'yule' or 'pda' as argument values")
    return colless

def pybus_harvey_gamma(tree, prec=0.00001):
    """Returns the gamma statistic of Pybus and Harvey (2000). This statistic
    is used to test for constancy of birth and death rates over the course of
    a phylogeny.  Under the pure-birth process, the statistic should follow
    a standard Normal distibution: a Normal(mean=0, variance=1).

    If the lengths of different paths to the node differ by more than ``prec``,
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    Raises a Value Error if the tree is not ultrametric, is non-binary, or has
        only 2 leaves.

    As a side effect a ``age`` attribute is added to the nodes of the tree.

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
    if tree.seed_node.age is None:
        tree.calc_node_ages(ultrametricity_precision=prec)
    for node in tree.postorder_node_iter():
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
    for i in range(2, n):
        list_index = i - 2
        T += i * float(g[list_index])
        accum += T
    list_index = n - 2
    T += (n) * g[list_index]
    nmt = n - 2.0
    numerator = accum/nmt - T/2.0
    C = T*pow(1/(12*nmt), 0.5)
    return numerator/C

def N_bar(tree):
    """
    Returns the $\bar{N}$ statistic: the average number of nodes above a
    terminal node.
    """
    leaf_count = 0
    nbar = 0
    for leaf_node in tree.leaf_node_iter():
        leaf_count += 1
        for parent in leaf_node.ancestor_iter(inclusive=False):
            nbar += 1
    return float(nbar) / leaf_count

def sackin_index(tree, normalize=True):
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
    for leaf_node in tree.leaf_node_iter():
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
        raise TypeError("``normalization`` accepts only None, True, False, 'yule' or 'pda' as argument values")
    return s

def treeness(tree):
    """
    Returns the proportion of total tree length that is taken up by
    internal branches.
    """
    internal = 0.0
    external = 0.0
    for nd in tree.postorder_node_iter():
        if not nd._parent_node:
            continue
        if nd.is_leaf():
            external += nd.edge.length
        else:
            internal += nd.edge.length
    return internal/(external + internal)
























import sys
import copy

class BandeltNode:
    """
    This is just an auxiliary class that does Bandelt encoding / decoding.
    Users don't need to access this class. 
    """
    def __init__(self, data):
        self.parent = None
        self.left = None
        self.right = None
        self.data = data
        
    def find_node(self, val):
        if self.data == val:
            return self
        else:
            if (self.left == None) and (self.right == None):
                return None
            if self.left != None:
                find_left = self.left.find_node(val)
                if find_left != None:
                    return find_left
            if self.right != None:
                find_right = self.right.find_node(val)
                if find_right != None:
                    return find_right
            return None
    
    # This is a helper function to check whether the tree is correctly constructed.
    def print_details(self):
        parent_data = self.parent.data if (self.parent != None) else self.parent
        left_data = self.left.data if (self.left != None) else self.left
        right_data = self.right.data if (self.right != None) else self.right
        print("Current Value: ", self.data, "; Parent: ", parent_data, "; left: ", left_data, "; right: ", right_data)
        
    # This is a helper function to check whether the tree is correctly constructed.
    def print_subtree(self, indent_num = 1):
        if indent_num == 1:
            print(self.data)
        if self.left != None:
            print('___'*indent_num, self.left.data)
            self.left.print_subtree(indent_num + 1)
        if self.right != None:
            print('___'*indent_num, self.right.data)
            self.right.print_subtree(indent_num + 1)
            
    # This is a helper function to check encoded & decoded tree are same. 
    def compare_subtree(self, compared_root_node):
        if self.data != compared_root_node.data:
            return False
        self_left_data = self.left.data if (self.left != None) else None
        self_right_data = self.right.data if (self.right != None) else None
        compared_root_node_left_data = compared_root_node.left.data if (compared_root_node.left != None) else None
        compared_root_node_right_data = compared_root_node.right.data if (compared_root_node.right != None) else None
        if (self_left_data in [compared_root_node_left_data, compared_root_node_right_data]) and (self_right_data in [compared_root_node_left_data, compared_root_node_right_data]):
            if (self_left_data == compared_root_node_left_data):
                if self.left != None:
                    compare_ans_left = self.left.compare_subtree(compared_root_node.left)
                    if not compare_ans_left:
                        return False
                if self.right != None:
                    compare_ans_right = self.right.compare_subtree(compared_root_node.right)
                    if not compare_ans_right:
                        return False
                
            elif (self_left_data == compared_root_node_right_data):
                if self.left != None:
                    compare_ans_left = self.left.compare_subtree(compared_root_node.right)
                    if not compare_ans_left:
                        return False
                if self.right != None:
                    compare_ans_right = self.right.compare_subtree(compared_root_node.left)
                    if not compare_ans_right:
                        return False
        else:
            return False
        return True
    

def Bandelt_encode(tree):
    """
    Returns (1) Bandelt encoded list, (2)encoding dictionary and (3)decoding dictionary.
    """
    def create_Bandelt_tree(parent, parent_BN_node):    
        for idx, child in enumerate(parent.child_nodes()):
            if child.taxon is not None:
            ## leaf nodes
                child_BN_node = BandeltNode(encode_dict[child.taxon.label.replace(' ', '_')])
            else:
            ## inner nodes
                child_BN_node = BandeltNode(-100)
            if idx == 0:
                parent_BN_node.left = child_BN_node
                child_BN_node.parent = parent_BN_node
            elif idx == 1:
                parent_BN_node.right = child_BN_node
                child_BN_node.parent = parent_BN_node
            create_Bandelt_tree(child, child_BN_node)
            
    def inner_node_indexation(current_node):
        if current_node.left != None:
            left_data = inner_node_indexation(current_node.left)
        if current_node.right != None:
            right_data = inner_node_indexation(current_node.right)
        if current_node.left != None and current_node.right != None:
            # This is the inner node
            if len(left_data) == 0:
                current_node.data = int(-right_data[0])
                return []
            if len(right_data) == 0:
                current_node.data = int(-left_data[0])
                return []
            if len(left_data) != 0 and len(right_data) != 0:
                if abs(left_data[0]) > abs(right_data[0]):
                    current_node.data = int(-left_data[0])
                    return [right_data[0]]
                elif abs(left_data[0]) < abs(right_data[0]):
                    current_node.data = int(-right_data[0])
                    return [left_data[0]]
        elif current_node.left != None and current_node.right == None:
            # Finish!
            pass
        else:
            return [current_node.data]
        
    def find_Bandelt_encode(target_node, val):
        visited_nodes = []
        queue = []
        visited_nodes.append(target_node.data)
        queue.append(target_node)
        encode_num = None
        while queue:
            current_node = queue.pop(0)
            if current_node != None:
                if current_node.data != '*' and abs(int(current_node.data)) < abs(int(val)):
                    encode_num = current_node.data
            if current_node.left != None:
                if current_node.left.data not in visited_nodes:
                    visited_nodes.append(current_node.left.data)
                    queue.append(current_node.left)
            if current_node.right != None:
                if current_node.right.data not in visited_nodes:
                    visited_nodes.append(current_node.right.data)
                    queue.append(current_node.right)   
            if encode_num != None:
                break
        return encode_num
    
    # Get encoding dictionary & decoding dictionary first
    encode_dict = {}
    decode_dict = {}
    for idx, lves in enumerate(tree.leaf_nodes()):
        if idx == 0:
            encode_dict[lves.taxon.label.replace(' ', '_')] = sys.maxsize
            decode_dict[sys.maxsize] = lves.taxon.label.replace(' ', '_')
        else:
            encode_dict[lves.taxon.label.replace(' ', '_')] = idx-1
            decode_dict[idx-1] = lves.taxon.label.replace(' ', '_')

    #######################################################
    #######################################################
    ## Tree checking conditions should be added HERE !! ##
    #######################################################
    #######################################################
    # 1. Bifurcating tree (IQ-Tree output as the input of this function)
    # 
    # What I want to do here is to choose the first leaf node as the root node but I don't know whether there is a better way to 
    # write the checking conditions. Please help me out here. 
    
    # For example, this is the tree created by dendropy with IQ-Tree output. 
    #
    #  (1) Tree 1
    #   /---------------------------------------------------- A
    #   |
    #   +---------------------------------------------------- B
    #   |
    #   |                /----------------------------------- C
    #   \----------------+
    #                    |                 /----------------- D
    #                    \-----------------+
    #                                      \----------------- E
    #
    #  (1) Tree 2
    # I want to construct the tree like the structure below with my own auxiliary BandeltNode class.
    #
    #                /---------------------------------------------------- B
    #   A -----------+
    #                |                /----------------------------------- C
    #                \----------------+
    #                                 |                 /----------------- D
    #                                 \-----------------+
    #                                                   \----------------- E
    #
    # Now I assume the read-in tree will look like Tree 1 and let the first child of seed_node be the root node of my own 
    # Bandelt tree (BandeltNode). 
    
    # Condition: The number of seed_node's children node must be 3 and the first child must be a leaf node which 
    #            wil be used as the root node.
    #if (len(tree.seed_node.child_nodes()) == 3 && 
    #    encode_dict[tree.seed_node.child_nodes()[0].taxon.label.replace(' ', '_')] == sys.maxsize):
    #    pass
            
    root_node = BandeltNode(encode_dict[tree.seed_node.child_nodes()[0].taxon.label.replace(' ', '_')])
    root_inner_node = BandeltNode(-100)
    root_node.left = root_inner_node
    root_inner_node.parent = root_node
    
    root_child_nodes = tree.seed_node.child_nodes()

    if (root_child_nodes[1].taxon is not None):
        sec_layer_left = BandeltNode(encode_dict[root_child_nodes[1].taxon.label.replace(' ', '_')])
    else:
        sec_layer_left = BandeltNode(-100)

    if (root_child_nodes[2].taxon is not None):
        sec_layer_right = BandeltNode(encode_dict[root_child_nodes[2].taxon.label.replace(' ', '_')])
    else:
        sec_layer_right = BandeltNode(-100)

    root_inner_node.left = sec_layer_left
    sec_layer_left.parent = root_inner_node
    root_inner_node.right = sec_layer_right
    sec_layer_right.parent = root_inner_node
    create_Bandelt_tree(tree.seed_node.child_nodes()[1], root_node.left.left)
    create_Bandelt_tree(tree.seed_node.child_nodes()[2], root_node.left.right)
    inner_node_indexation(root_node)
    encode_num_list = []
    for i in range(0, len(tree.leaf_nodes())-2):
        target_node_val = -(i+1)
        found_node = root_node.find_node(target_node_val)
        encode_num = find_Bandelt_encode(found_node, target_node_val)
        encode_num_list.append(encode_num)
        
    return(encode_num_list, encode_dict, decode_dict)



def Bandelt_decode(Bandelt_Encode_list, decode_dict):
    """
    Returns Bandelt decode newick string.
    """
    def post_order_traversal_num_2_name(current_node, decode_dict):
        if current_node.left != None:
            left_newick_string = post_order_traversal_num_2_name(current_node.left, decode_dict)
        if current_node.right != None:
            right_newick_string = post_order_traversal_num_2_name(current_node.right, decode_dict)
        if current_node.data < 0:
            current_node.data = ''
        elif current_node.data >= 0:
            current_node.data = decode_dict[int(current_node.data)]

    def post_order_traversal_newick_string(current_node):
        if current_node.left != None:
            left_newick_string = post_order_traversal_newick_string(current_node.left)
        if current_node.right != None:
            right_newick_string = post_order_traversal_newick_string(current_node.right)

        if current_node.left == None and current_node.right == None:     
            # This node is the terminal vertex
            return str(current_node.data)

        if current_node.left != None and current_node.right != None:
            newick_string = '(' + left_newick_string + ',' + right_newick_string + ')' + str(current_node.data)

        if current_node.left != None and current_node.right == None:
            newick_string = '(' + str(current_node.data) + ',' + left_newick_string + ')'

        return newick_string
    
    BANDELT_NUM = len(Bandelt_Encode_list)
    # Initial with three nodes
    root_node = BandeltNode(sys.maxsize)
    node_1 = BandeltNode(1)
    node_neg_1 = BandeltNode(-1)
    node_0 = BandeltNode(0)
    # Create links between initial three nodes
    root_node.left = node_neg_1
    node_neg_1.parent = root_node
    node_neg_1.left = node_0
    node_neg_1.right = node_1
    node_1.parent = node_neg_1
    node_0.parent = node_neg_1
    
    for i in range(1, BANDELT_NUM):
        added_node_val = -(i+1)
        added_node_neg = BandeltNode(added_node_val)
        added_node_pos = BandeltNode(-added_node_val)

        target_node = root_node.find_node(Bandelt_Encode_list[i])

        # If target node is the left child 
        if target_node.parent.left == target_node:
            target_node.parent.left = added_node_neg
            added_node_neg.parent = target_node.parent
            added_node_neg.left = target_node
            target_node.parent = added_node_neg

            added_node_neg.right = added_node_pos
            added_node_pos.parent = added_node_neg

        # If target node is the right child
        if target_node.parent.right == target_node:
            target_node.parent.right = added_node_neg
            added_node_neg.parent = target_node.parent
            added_node_neg.right = target_node
            target_node.parent = added_node_neg

            added_node_neg.left = added_node_pos
            added_node_pos.parent = added_node_neg
    copy_root_node = copy.deepcopy(root_node)
    post_order_traversal_num_2_name(copy_root_node, decode_dict)
    decode_tree = post_order_traversal_newick_string(copy_root_node)
    decode_tree_newick = decode_tree + ';'
    return decode_tree_newick