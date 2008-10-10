#! /usr/bin/env python

############################################################################
##  treedists.py
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
Tree distances.
"""

from dendropy import splits

def robinson_foulds_calc(length_diffs):
    """
    Given, `length_diffs`, a list of pairs of corresponding (length/weight) values
    of edges from two trees, this returns the Robinson-Foulds distance (sum of
    absolute differences) between the two trees.
    """
    return sum([abs(length_diffs[i][0]-length_diffs[i][1]) for i in range(len(length_diffs))])

def brlen_scores_calc(length_diffs):
    """
    Given, `length_diffs`, a list of pairs of corresponding (length/weight) values
    of edges from two trees, this returns the branch length score (sum of
    squared differences) between the two trees. This is equivalent to the squared
    Euclidean distance between the two trees.
    """
    return sum([pow(length_diffs[i][0]-length_diffs[i][1], 2) for i in range(len(length_diffs))])

def brlen_dists_calc(length_diffs):
    """
    Given, `length_diffs`, a list of pairs of corresponding (length/weight) values
    of edges from two trees, this returns the branch length distance (square root of the 
    sum of squared differences) between the two trees. This is equivalent to the Euclidean
    branch length distance between the two trees.
    """
    return pow(brlen_scores_calc(length_diffs), 0.5)

def splits_distance(tree1, 
                    tree2, 
                    dist_func=robinson_foulds_calc, 
                    edge_length_attr="length",
                    value_type=float):
    """
    Returns distance between two trees, each represented by a dictionary of
    splits (as split_mask strings) to edges, using `dist_func` to calculate the 
    distance based on `edge_length_attr` of the edges. `dist_func` is a function
    that takes a list of pairs of values, where the values correspond to the edge
    lengths of a given split on tree1 and tree2 respectively.
    """
    length_diffs = []
    split_set = set(tree1.split_edges.keys())
    split_set.update(tree2.split_edges.keys())
    for split in split_set:
    
        if split in tree1.split_edges and getattr(tree1.split_edges[split], edge_length_attr):
            value1 = value_type(getattr(tree1.split_edges[split], edge_length_attr))
        elif split in tree1.complemented_split_edges and getattr(tree1.complemented_split_edges[split], edge_length_attr):
            value1 = value_type(getattr(tree1.complemented_split_edges[split], edge_length_attr))            
        else:
            value1 = value_type(0)
            
        if split in tree2.split_edges and getattr(tree2.split_edges[split], edge_length_attr):
            value2 = value_type(getattr(tree2.split_edges[split], edge_length_attr))
        elif split in tree2.complemented_split_edges and getattr(tree2.complemented_split_edges[split], edge_length_attr):
            value1 = value_type(getattr(tree2.complemented_split_edges[split], edge_length_attr))                   
        else:
            value2 = value_type(0)
            
        length_diffs.append((value1,value2))
    return dist_func(length_diffs)
        
def robinson_foulds_distance(tree1, tree2, edge_length_attr="length"):
    """
    Returns Robinson-Foulds distance between two trees based on `edge_length_attr`. 
    Trees need to have been decorated with the `encode_splits` method of the splits 
    module.
    """
    return splits_distance(tree1, 
                           tree2, 
                           dist_func=robinson_foulds_calc,
                           edge_length_attr=edge_length_attr,
                           value_type=float)

def euclidean_distance(tree1, tree2, edge_length_attr="length", value_type=float):
    """
    Returns Euclidean distance (a.k.a. Felsenstein's 2004 `branch length distance`) 
    between two trees based on `edge_length_attr`. 
    Trees need to have been decorated with the `encode_splits` method of the splits 
    module.
    """
    return splits_distance(tree1, 
                           tree2, 
                           dist_func=brlen_dists_calc,
                           edge_length_attr=edge_length_attr,
                           value_type=value_type)
                           
def symmetric_difference(tree1, tree2):
    """
    False pos = splits in tree2 NOT in tree1
    False neg = splits in tree1 NOT in tree2
    """
    sym_diff = 0
    false_positives = 0
    false_negatives = 0
    
    for split in tree1.splits:
        if (split in tree2.splits) or (split in tree2.complemented_splits):
            pass
        else:
            false_negatives = false_negatives + 1
            sym_diff = sym_diff + 1
    
    for split in tree2.splits:
        if (split in tree1.splits) or (split in tree1.complemented_splits):
            pass
        else:
            false_positives = false_positives + 1
            sym_diff = sym_diff + 1   
            
    return sym_diff, false_positives, false_negatives
    

import sys
import newick
if __name__ == "__main__":
    reader = newick.NewickTreeReader()
    writer = newick.NewickTreeWriter()
    
#     trees = reader.read_trees(text='(3,(4,(1,2))); (2,(1,(3,4))); ')
#     print "Distance: %d" % symmetric_difference(trees[0], trees[1])[0]    
    
   # trees = reader.read_trees(text='(3,(4,(1,2))); (2,(1,(3,4))); (1,(2,(3,4))); ((1,2),(3,4)); ((3,4),(1,2)); ((1,4),(3,2));')      
    trees = reader.read_trees(text='(1,2,(3,(6,(4,5)))); (((1, 2),3),(4,5),6));')
    for tree in trees:
        splits.encode_splits(tree, trees.taxa_block)    
    
    for tree1 in trees:
        for tree2 in trees:
            distance = symmetric_difference(tree1, tree2)[0]
            print writer.compose_tree(tree1), "vs.", writer.compose_tree(tree2), "= %d" % distance
            if distance:
                print tree1.splits
                print tree1.complemented_splits
                print tree2.splits
                print tree2.complemented_splits
                sys.exit(1)