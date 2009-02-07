#! /usr/bin/env python

############################################################################
##  charevolve.py
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
Evolves characters on tree.
"""

import copy

from dendropy import GLOBAL_RNG
from dendropy import datasets
from dendropy import characters
from dendropy import charmodels

############################################################################
## Convenience wrappers

def generate_hky_dataset(seq_len,
                         tree_model,                   
                         mutation_rate=1.0, 
                         kappa=1.0,
                         base_freqs=[0.25, 0.25, 0.25, 0.25],
                         root_states=None,    
                         dataset=None,
                         taxa_block=None,
                         rng=None):
    """
    Convenience class to wrap generation of a dataset based on
    the HKY model.
    `seq_len`       : length of sequence (number of characters)
    `tree_model`    : dendropy.trees.Tree object
    `char_model`    : dendropy.charmodels.CharacterModel object
    `mutation_rate` : mutation *modifier* rate (should be 1.0 if branch lengths
                      on tree reflect true expected number of changes
    `root_states`   : vector of root states (length must equal `seq_len`)
    `dataset`       : a dendropy.datasets.Dataset object.
                      if given, the new 
                      dendropy.characters.CharacterBlock object will be added to
                      this (along with a new taxa_block if required). Otherwise,
                      a new dendropy.datasets.Dataset object will be created.
    `taxa_block`    : if given, this will be the taxa manager used; otherwise
                      new default one will be created
    `rng`           : random number generator; if not given, `GLOBAL_RNG` will be
                      used                      
    Returns: a dendropy.datasets.Dataset object object.
    """
    char_model = charmodels.Hky85CharacterModel(kappa=kappa, 
                                                base_freqs=base_freqs)
    return generate_dataset(seq_len=seq_len,
        tree_model=tree_model,
        char_model=char_model,
        mutation_rate=mutation_rate,       
        root_states=root_states, 
        dataset=dataset,
        taxa_block=taxa_block,
        rng=rng)  
        
def generate_hky_characters(seq_len,
                            tree_model,                   
                            mutation_rate=1.0, 
                            kappa=1.0,
                            base_freqs=[0.25, 0.25, 0.25, 0.25],
                            root_states=None,    
                            char_block=None,
                            taxa_block=None,
                            rng=None):
    """
    Convenience class to wrap generation of characters (as a CharacterBlock
    object) based on the HKY model.
    `seq_len`       : length of sequence (number of characters)
    `tree_model`    : dendropy.trees.Tree object
    `mutation_rate` : mutation *modifier* rate (should be 1.0 if branch lengths
                      on tree reflect true expected number of changes
    `root_states`   : vector of root states (length must equal `seq_len`)
    `char_block`    : dendropy.characters.CharacterBlock object.
                      if given, new sequences for taxa on `tree_model` leaf_nodes 
                      will be appended to existing sequences of corresponding 
                      taxa in char_block; if not, a new 
                      dendropy.characters.CharacterBlock object will be created
    `taxa_block`    : if given, this will be the taxa manager used; otherwise
                      new default one will be created
    `rng`           : random number generator; if not given, `GLOBAL_RNG` will be
                      used                      
    Returns: a dendropy.characters.CharacterBlock object.    
    
    Since characters will be appended to existing sequences, you can simulate a
    sequences under a mixed model by calling this method multiple times with 
    different character model parameter values and/or different mutation
    rates, passing in the same `char_block` object each time.    
    """
    char_model = charmodels.Hky85CharacterModel(kappa=kappa, 
                                                base_freqs=base_freqs)
    return generate_characters(seq_len=seq_len,
                               tree_model=tree_model,
                               char_model=char_model,
                               mutation_rate=mutation_rate,       
                               root_states=root_states, 
                               char_block=char_block,
                               taxa_block=taxa_block,
                               rng=rng)    
                
def generate_dataset(seq_len,
                     tree_model,
                     char_model,
                     mutation_rate=1.0,       
                     root_states=None, 
                     dataset=None,
                     taxa_block=None,
                     rng=None):
    """
    Wrapper to conveniently generate a Dataset simulated under
    the given tree and character model.
    `seq_len`       : length of sequence (number of characters)
    `tree_model`    : dendropy.trees.Tree object
    `char_model`    : dendropy.charmodels.CharacterModel object
    `mutation_rate` : mutation *modifier* rate (should be 1.0 if branch lengths
                      on tree reflect true expected number of changes
    `root_states`   : vector of root states (length must equal `seq_len`)
    `dataset`       : a dendropy.datasets.Dataset object.
                      if given, the new 
                      dendropy.characters.CharacterBlock object will be added to
                      this (along with a new taxa_block if required). Otherwise,
                      a new dendropy.datasets.Dataset object will be created.
    `taxa_block`    : if given, this will be the taxa manager used; otherwise
                      new default one will be created
    `rng`           : random number generator; if not given, `GLOBAL_RNG` will be
                      used                      
    Returns: a dendropy.datasets.Dataset object object.
    """
    if dataset is None:
        dataset = datasets.Dataset()
    if taxa_block is not None and taxa_block not in dataset.taxa_blocks:
        taxa_block = dataset.add_taxa_block(taxa_block=taxa_block)
    char_block = generate_characters(seq_len=seq_len,
        tree_model=tree_model,
        char_model=char_model,
        mutation_rate=mutation_rate,       
        root_states=root_states, 
        char_block=None,
        taxa_block=taxa_block,
        rng=None)
    dataset.add_char_block(char_block=char_block)
    return dataset
    
def generate_characters(seq_len,
                        tree_model,
                        char_model,
                        mutation_rate=1.0,       
                        root_states=None, 
                        char_block=None,
                        taxa_block=None,
                        rng=None):
    """
    Wrapper to conveniently generate a characters simulated under
    the given tree and character model.
    `seq_len`       : length of sequence (number of characters)
    `tree_model`    : dendropy.trees.Tree object
    `char_model`    : dendropy.charmodels.CharacterModel object
    `mutation_rate` : mutation *modifier* rate (should be 1.0 if branch lengths
                      on tree reflect true expected number of changes
    `root_states`   : vector of root states (length must equal `seq_len`)
    `char_block`    : dendropy.characters.CharacterBlock object.
                      if given, new sequences for taxa on `tree_model` leaf_nodes 
                      will be appended to existing sequences of corresponding 
                      taxa in char_block; if not, a new 
                      dendropy.characters.CharacterBlock object will be created
    `taxa_block`    : if given, this will be the taxa manager used; otherwise
                      new default one will be created
    `rng`           : random number generator; if not given, `GLOBAL_RNG` will be
                      used
                      
    Returns: a dendropy.characters.CharacterBlock object.
    
    Since characters will be appended to existing sequences, you can simulate a
    sequences under a mixed model by calling this method multiple times with 
    different character models and/or different mutation rates, passing 
    in the same `char_block` object each time.
    """
    char_evolver = CharEvolver(char_model=char_model,
                               mutation_rate=mutation_rate)
    tree = char_evolver.evolve_states(tree=tree_model,
        seq_len=seq_len,
        root_states=None,
        rng=rng)
    char_matrix = char_evolver.compose_char_matrix(tree, taxa_block)
    if char_block is None:
        char_block = characters.DnaCharactersBlock()
        if taxa_block is not None:
            char_block.taxa_block = taxa_block            
    else:
        if taxa_block is not None:
            char_block.taxa_block = taxa_block  
    char_block.extend_matrix(other_matrix=char_matrix, 
        overwrite_existing=False, 
        append_existing=True)
    return char_block
    
############################################################################
## Workhorse class(es)
                             
class CharEvolver(object):
    "Evolves sequences on a tree."

    def __init__(self,
     char_model=None,
     mutation_rate=None,
     seq_attr='sequences',
     char_model_attr="char_model",                 
     edge_length_attr="length",
     edge_rate_attr="mutation_rate",
     seq_label_attr='taxon'):
        "Sets up meta-data dealing with object nomenclature and semantics."
        self.char_model = char_model
        self.mutation_rate = mutation_rate
        self.seq_attr = seq_attr
        self.char_model_attr = char_model_attr        
        self.edge_length_attr = edge_length_attr
        self.edge_rate_attr = edge_rate_attr
        self.seq_label_attr = seq_label_attr
    
    def evolve_states(self,
     tree,
     seq_len,
     root_states=None,
     generate_root_states=True,
     in_place=False,
     rng=None):
        """
        Appends a new sequence of length `seq_len` to a list at each node
        in `tree`.  The attribute name of this list in each node is given
        by `seq_attr`. If `char_model` is None, `tree.char_model` or
        `char_model` at each node must be specified. If `in_place` is
        False, the tree is copied first, otherwise original tree is modified.
        If `root_states` is given, this will be used as the sequence for the root.
        If not, and if `generate_root_states` is True, then the sequence for the 
        root will be drawn from the stationary distribution of the character model.
        """
        if rng is None:
            rng = GLOBAL_RNG
        if not in_place:            
            tree = copy.deepcopy(tree)

        if self.char_model is None:
            char_model = getattr(tree, self.char_model_attr, None)

        # loop through edges in preorder (root->tips)
        for edge in tree.preorder_edge_iter():
            node = edge.head_node
            if not hasattr(node, self.seq_attr):
                setattr(node, self.seq_attr, [])
            seq_list = getattr(node, self.seq_attr)
            if edge.tail_node:
                par = edge.tail_node                
                if len(seq_list) != n_prev_seq:
                    raise ValueError("'%s' length varies among nodes" % self.seq_attr)
                par_seq = getattr(par, self.seq_attr)[-1]
                char_model  = getattr(edge, self.char_model_attr, None) or self.char_model
                length = getattr(edge, self.edge_length_attr)
                mutation_rate = getattr(edge, self.edge_rate_attr, None) or self.mutation_rate
                seq_list.append(char_model.generate_descendant_states(par_seq, length, mutation_rate))
            else:
                # no tail node: root
                n_prev_seq = len(seq_list)
                if root_states is not None:
                    seq_list.append(root_states)
                elif generate_root_states:
                    char_model  = getattr(node.edge, self.char_model_attr, None) or self.char_model 
                    seq_list.append(char_model.stationary_sample(seq_len, rng=rng))
                else:
                    assert n_prev_seq > 0
                    n_prev_seq -= 1
        return tree
        
    def compose_char_matrix(self, tree, taxa_block=None, include=None, exclude=None):
        """
        Returns a CharacterDataMatrix where the keys are the taxa of the leaf_nodes
        of `source_tree` and the values are sequences, with each sequence being the
        concatentation of all the sequences in the list of sequences associated with
        each tip. Specific sequences to be included/excluded can be fine-tuned using
        the `include` and `exclude` args, where `include`=None means to include all
        by default, and `exclude`=None means to exclude all by default.
        """
        char_matrix = characters.CharacterDataMatrix()
        for leaf in tree.leaf_nodes():
            cvec = characters.CharacterDataVector(taxon=leaf.taxon)
            seq_list = getattr(leaf, self.seq_attr)
            for seq_idx, seq in enumerate(seq_list):
                if ((include is None) or (seq_idx in include))  \
                    and ((exclude is None) or (seq_idx not in exclude)):
                    for state in seq:
                        cvec.append(characters.CharacterDataCell(value=state))
            if taxa_block is not None:
                taxon = taxa_block.get_taxon(label=leaf.taxon.label)
            else:
                taxon = leaf.taxon
            char_matrix[taxon] = cvec                                
        return char_matrix        
                        
#     def generate_char_matrix(self, tree, seq_len=1000):
#         """
#         Evolves sequences on tree using default settings, and returns
#         CharacterDataMatrix where keys are taxon (objects) and values
#         are lists of sequence states.
#         """
#         mtree = self.evolve_states(tree=tree, seq_len=seq_len)
#         char_matrix = characters.CharacterDataMatrix()
#         for leaf in tree.leaf_nodes():
#             cvec = characters.CharacterDataVector(taxon=leaf.taxon)
#             states = getattr(leaf, self.seq_attr)[-1]
#             for state in states:
#                 cvec.append(characters.CharacterDataCell(value=state))
#             if taxa_block is not None:
#                 taxon = taxa_block.get_taxon(label=leaf.taxon.label)
#             else:
#                 taxon = leaf.taxon
#             char_matrix[taxon] = cvec                              
#         return char_matrix


