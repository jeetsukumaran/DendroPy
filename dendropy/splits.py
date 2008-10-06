#! /usr/bin/env python

############################################################################
##  splits.py
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
##  with this programm. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Split calculation and management.
"""

from dendropy import taxa
from dendropy import trees
from dendropy import treegen

def is_non_singleton_split(split):
    """
    Returns True if a split is NOT between a leaf and the rest of the taxa.
    """
    # ((split-1) & split) is True (non-zero) only
    # if split is not a power of 2, i.e., if split
    # has more than one bit turned on, i.e., if it
    # is a non-trivial split.    
    return ((split-1) & split)

def number_to_bitstring(num):
    """
    Returns a representation of a number as a bit string.
    """
    return num>0 and number_to_bitstring(num>>1)+str(num&1) or ''
    
def split_as_string(split_mask, taxa_block, symbol1='.', symbol2='*'):
    """
    Returns a 'pretty' split representation.
    """
    s = number_to_bitstring(split_mask).rjust(len(taxa_block), '0')
    return s.replace('0', symbol1).replace('1', symbol2)
    
def split_as_string_rev(split_mask, taxa_block, symbol1='.', symbol2='*'):
    """
    Returns a 'pretty' split representation, reversed, with first taxon
    on the left (as given by PAUP*)
    """
    return split_as_string(split_mask=split_mask, 
                           taxa_block=taxa_block, 
                           symbol1=symbol1,
                           symbol2=symbol2)[::-1]
    
def split_taxa_list(split_mask, taxa_block, index=0):
    """
    Returns list of taxa represented by split.
    """
    taxa = []
    while split_mask:
        if split_mask & 1:
            taxa.append(taxa_block[index])
        split_mask = split_mask >> 1
        index += 1
    return taxa

def encode_splits(tree, 
                  taxa_block=None, 
                  edge_split_mask='split_mask', 
                  tree_split_edges_map="split_edges",
                  tree_split_taxa_map="split_taxa",
                  tree_splits_list="splits",
                  tree_complemented_splits='complemented_splits',
                  tree_complemented_split_edges_map="complemented_split_edges"):
    """
    Processes splits on a tree, encoding them as bitmask on each edge. 
    Returns a dictionary where the keys are splits and the values are edges.
    If `tree_split_edges_map` is given, then the dictionary is embedded as an attribute
    of the tree.
    If `tree_split_taxa_map` is given, then an additional tree attribute is set---
    a dictionary of splits to list of taxa represented by the split.
    If `tree_splits_list` is given, a set of all splits is added to the tree as an attribute
    If `tree_complemented_split_edges_maps` is  gvien, a dictionary of all complements of the splits 
    is also added, with the complemented splits as keys and the corresponding edges as values
    """
    if taxa_block is None:
        taxa_block = tree.infer_taxa_block()
    split_map = {}
    for edge in tree.postorder_edge_iter():
        children = edge.head_node.children()
        if children:
            setattr(edge, edge_split_mask, 0)
            for child in children:
                setattr(edge, edge_split_mask, getattr(edge, edge_split_mask) | getattr(child.edge, edge_split_mask))
        else:
            if edge.head_node.taxon:
                setattr(edge, edge_split_mask, taxa_block.taxon_bitmask(edge.head_node.taxon))
            else:
                setattr(edge, edge_split_mask, 0)
        split_map[getattr(edge, edge_split_mask)] = edge
    if tree_split_edges_map:
        setattr(tree, tree_split_edges_map, split_map)
    if tree_split_taxa_map:
        split_taxa = {}
        for split in split_map:
            split_taxa[split] = split_taxa_list(split, taxa_block)
        setattr(tree, tree_split_taxa_map, split_taxa)
    if tree_splits_list:
        setattr(tree, tree_splits_list, set(split_map.keys()))
    if tree_complemented_splits or tree_complemented_split_edges_map:
        taxa_mask = taxa_block.all_taxa_bitmask()
        csplit_edges = {}
        for split in split_map.keys():
            csplit = split ^ taxa_mask
            csplit_edges[csplit] = split_map[split]
        if tree_complemented_splits:
            setattr(tree, tree_complemented_splits, set(csplit_edges.keys()))
        if tree_complemented_split_edges_map:            
            setattr(tree, tree_complemented_split_edges_map, csplit_edges)        
    return split_map
                
############################################################################        
## SplitDistribution

class SplitDistribution(object):
    """
    Collects information regarding splits over multiple trees.
    """
    
    def __init__(self, taxa_block=None):
        """
        What else?
        """
        self.total_trees_counted = 0
        if taxa_block is not None:
            self.taxa_block = taxa_block            
        else:
            self.taxa_block = taxa.TaxaBlock()
        self.splits = []
        self.complemented_splits = []
        self.split_counts = {}
        self.complemented_split_counts = {}
        self.split_edge_lengths = {}
        self.split_node_ages = {}
        self.ignore_edge_lengths = False
        self.ignore_node_ages = False
        self.__split_freqs = None
        self.__complemented_split_freqs = None
        self.__trees_counted_for_freqs = 0
        
    def splits_considered(self):
        """
        Returns 4 values:
            total number of splits counted
            total number of unique splits counted
            total number of non-trivial splits counted
            total number of unique non-trivial splits counted
        """
        num_splits = 0
        num_unique_splits = 0
        num_nt_splits = 0
        num_nt_unique_splits = 0
        for s in self.split_counts:
            num_unique_splits += 1 
            num_splits += self.split_counts[s]
            if is_non_singleton_split(s):
                num_nt_unique_splits += 1
                num_nt_splits += self.split_counts[s]
        return num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits                                                        
        
    def calc_freqs(self):
        """
        Forces recalculation of frequencies.
        """
        self.__split_freqs = {}
        self.__complemented_split_freqs = {}
        if self.total_trees_counted == 0:
            total = 1
        else:
            total = self.total_trees_counted
        for split in self.split_counts:
            self.__split_freqs[split] = float(self.split_counts[split]) / total
        for split in self.complemented_split_counts:
            self.__complemented_split_freqs[split] = float(self.complemented_split_counts[split]) / total
        self.__trees_counted_for_freqs = self.total_trees_counted            
        return self.__split_freqs, self.__complemented_split_freqs                
        
    def _get_split_frequencies(self):
        """
        Returns dictionary of splits : split frequencies.
        """
        if self.__split_freqs is None or self.__trees_counted_for_freqs != self.total_trees_counted:
            self.calc_freqs()
        return self.__split_freqs   
        
    split_frequencies = property(_get_split_frequencies)     
    
    def _get_complemented_split_frequencies(self):
        """
        Returns dictionary of complemented splits : split frequencies.
        """
        if self.__complemented_split_freqs is None or self.__trees_counted_for_freqs != self.total_trees_counted:
            self.calc_freqs()
        return self.__complemented_split_freqs   
        
    complemented_split_frequencies = property(_get_complemented_split_frequencies)      

    def count_splits_on_tree(self, tree):
        """
        Counts splits in this tree and add to totals.
        """
        self.total_trees_counted += 1
        tree.normalize_taxa(taxa_block=self.taxa_block)
        tree.deroot()
        encode_splits(tree, 
                     self.taxa_block, 
                     tree_split_taxa_map=None)        
        for split in tree.splits:
            if split not in self.split_counts:
                self.splits.append(split)
                self.split_counts[split] = 1
            else:
                self.split_counts[split] += 1     
            if not self.ignore_edge_lengths:                
                if split not in self.split_edge_lengths:
                    self.split_edge_lengths[split] = []    
                if tree.split_edges[split].length is not None: 
                    self.split_edge_lengths[split].append(tree.split_edges[split].length)
                else:
                    self.split_edge_lengths[split].append(0.0)
            if not self.ignore_node_ages:  
                if split not in self.split_node_ages:
                    self.split_node_ages[split] = []
                edge = tree.split_edges[split]
                if edge.head_node is not None:
                    self.split_node_ages[split].append(edge.head_node.distance_from_tip())
        for split in tree.complemented_splits:
            if split not in self.complemented_split_counts:
                self.complemented_splits.append(split)            
                self.complemented_split_counts[split] = 1
            else:
                self.complemented_split_counts[split] += 1    

############################################################################        
## SplitDistribution

class TreeSplitSummarizer(object):
    """
    Summarizes a distribution of splits on trees.
    """
    
    def deepest_compatible_node(start_node, split, taxa_mask):
        """
        Assumes:
            - start_node is a node on a tree
            - edges on tree have been decorated with splits
            - start_node is compatible with split (root if all else fails)
        Returns:
            - deepest node in tree (i.e., furthest from root) that is
              compatible with split (i.e, includes all taxa in split) yet
              not equal to the split.
        """
        for node in start_node.children():
            if (node.edge.split_mask ^ taxa_mask) & split:
                pass
            elif node.edge.split_mask != split:
                return TreeSplitSummarizer.deepest_compatible_node(node, split, taxa_mask)
        return start_node
        
    deepest_compatible_node = staticmethod(deepest_compatible_node)        
                    
    def __init__(self):
        """
        Initializes settings.
        """
        self.burnin = None 
        self.support_as_labels = True 
        self.support_as_percentages = False
        self.support_label_decimals = 2               
        self.verbose=False
        self.write_message = None# sys.stderr.write
        self.progress_message_prefix = None
        self.progress_message_suffix = None
        self.total_trees_read = 0    
        self.total_trees_counted = 0         
        self.total_trees_ignored = 0
        
    def send_progress_message(self, msg):
        """
        Writes progress message.
        """
        if self.verbose and self.write_message:
            if self.progress_message_prefix:
                prefix = self.progress_message_prefix
            else:
                prefix = ""
            if self.progress_message_suffix:
                suffix = self.progress_message_suffix
            else:
                suffix = ""                
            self.write_message("%s%s%s" % (prefix, msg, suffix))
            
    def compose_support_label(self, split_support_freq):
        """
        Returns an appropriately composed and formatted support label.
        """                
        if self.support_as_percentages:
            if self.support_label_decimals == 0:
                support_label = "%s" % (split_support_freq * 100)
            else:                
                support_label_template = "%%0.%df" % self.support_label_decimals
                support_label = support_label_template % (split_support_freq * 100)
        else:                
            if self.support_label_decimals < 0:
                support_label_decimals = 2
            else:
                support_label_decimals = self.support_label_decimals
            support_label_template = "%%0.%df" % self.support_label_decimals
            support_label = support_label_template % (split_support_freq * 100)                
        return support_label            
        
    def map_split_support_to_node(self, node, split_support):
        """
        Appropriately sets up a node.
        """
        if self.support_as_labels:
            node.label = self.compose_support_label(split_support)
        else:
            if self.support_as_percentages:
                node.edge.length = split_support * 100
            else:                        
                node.edge.length = split_support
        return node                
        
    def map_split_support_to_tree(self, tree, split_distribution):
        """
        Maps splits support to the given tree.
        """
        split_frequencies = split_distribution.split_frequencies
        complemented_split_frequencies = split_distribution.complemented_split_frequencies
        tree.normalize_taxa(taxa_block=split_distribution.taxa_block)
        splits.encode_splits(tree, 
                             taxa_block=split_distribution.taxa_block,
                             edge_split_mask='split_mask', 
                             tree_split_edges_map="split_edges",
                             tree_split_taxa_map=None,
                             tree_splits_list=None,
                             tree_complemented_splits=None,
                             tree_complemented_split_edges_map=None)
        for split in tree.split_edges:
            if split in split_frequencies:
                split_support = split_frequencies[split]
            elif split in complemented_split_frequencies:
                split_support = complemented_split_frequencies[split]
            else:
                split_support = 0.0
            self.map_split_support_to_node(tree.split_edges[split].head_node, split_support)                  
        return tree            
            
    def tree_from_splits(self, 
                         split_distribution, 
                         min_freq=0.5,
                         include_edge_lengths=True):
        """
        Returns a consensus tree constructed from splits given in `split_distribution`.
        """
        taxa_block = split_distribution.taxa_block
        con_tree = treegen.star_tree(taxa_block)
        split_freqs = split_distribution.split_frequencies
        taxa_mask = taxa_block.all_taxa_bitmask()
        for split in split_freqs:          
            if splits.is_non_singleton_split(split) \
                and (split ^ taxa_mask) \
                and (split_freqs[split] >= min_freq):            
                splits.encode_splits(con_tree, 
                                     taxa_block, 
                                     tree_split_edges_map=None,
                                     tree_split_taxa_map=None,
                                     tree_complemented_split_edges_map=None)
                parent_node = TreeSplitSummarizer.deepest_compatible_node(con_tree.seed_node, split, taxa_mask)
                new_node = trees.Node()
                self.map_split_support_to_node(node=new_node, split_support=split_freqs[split])                
                if self.support_as_labels and include_edge_lengths:
                        new_node.edge.length = float(sum(split_distribution.split_edge_lengths[split])) / len(split_distribution.split_edge_lengths[split])                
                for node in parent_node.children():
                    if (split ^ taxa_mask) & node.edge.split_mask:
                        pass
                    else:
                        parent_node.remove_child(node)
                        new_node.add_child(node)
                parent_node.add_child(new_node)

        ## here we add the support values and/or edge lengths for the terminal taxa ##
        for node in con_tree.leaves():
            split = node.edge.split_mask
            self.map_split_support_to_node(node, 1.0)
            if self.support_as_labels and include_edge_lengths:
                if split in split_distribution.split_edge_lengths:
                    node.edge.length = float(sum(split_distribution.split_edge_lengths[split])) / len(split_distribution.split_edge_lengths[split])
                else:
                    node.edge.length = 0.0
        return con_tree                
                                                
    def count_splits_on_trees(self, tree_files, tree_iterator, split_distribution=None, taxa_block=None):
        """
        Given a list of trees file, a SplitsDistribution object (a new one, or,
        if passed as an argument) is returned collating the split data in the files.
        """
        if split_distribution is None:
            split_distribution = splits.SplitDistribution()
        if taxa_block is not None:
            split_distribution.taxa_block = taxa_block            
#         if not isinstance(tree_files, list):
#             tree_files = [tree_files]
        total_tree_files = len(tree_files)
        for tree_file_idx, tree_file in enumerate(tree_files):
            file_trees_read = 0
            if isinstance(tree_file, str):
                tree_file_obj = open(tree_file, "r")
            else:
                tree_file_obj = tree_file             
            current_file_note = "Tree file %d of %d: " % (tree_file_idx+1, total_tree_files)
            for tree_idx, tree in enumerate(tree_iterator(src=open(tree_file, "r"))):
                self.total_trees_read += 1
                file_trees_read += 1
                if not self.burnin or file_trees_read > self.burnin:
                    self.total_trees_counted += 1
                    self.send_progress_message("%sCounting splits in tree %d" % (current_file_note, (tree_idx+1)))
                    #tree.normalize_taxa(taxa_block=split_distribution.taxa_block, update_taxa_block=True)
                    split_distribution.count_splits_on_tree(tree)
                else:
                    self.total_trees_ignored += 1
                    self.send_progress_message("%sSkipping tree %d (burn-in=%d)" % (current_file_note, (tree_idx+1), self.burnin))        
        return split_distribution                            