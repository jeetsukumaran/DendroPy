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
Models and modeling of discrete character evolution.
"""

import copy
import math
import itertools
from dendropy.utility import GLOBAL_RNG
from dendropy.mathlib import probability
import dendropy

############################################################################
## Character Evolution Modeling

class DiscreteCharacterEvolutionModel(object):
    "Base class for discrete character substitution models."

    def __init__(self, state_alphabet, rng=None):
        """
        __init__ initializes the state_alphabet to define the character type on which
        this model acts.  The objects random number generator will be `rng` or `GLOBAL_RNG`
        """
        self.state_alphabet = state_alphabet
        if rng is None:
            self.rng = GLOBAL_RNG
        else:
            self.rng = rng

    def pmatrix(self, tlen, rate=1.0):
        """
        Returns a matrix of nucleotide substitution
        probabilities.
        """
        raise NotImplementedError

    def simulate_descendant_states(self,
        ancestral_states,
        edge_length,
        mutation_rate=1.0,
        rng=None):
        """
        Returns descendent sequence given ancestral sequence.
        """
        if rng is None:
            rng = self.rng
        pmat = self.pmatrix(edge_length, mutation_rate)
        multi = probability.sample_multinomial
        desc_states = []
        for state in ancestral_states:
            anc_state_idx = self.state_alphabet.index(state)
            desc_state_idx = multi(pmat[anc_state_idx], rng)
            desc_states.append(self.state_alphabet[desc_state_idx])
        return desc_states

class DiscreteCharacterEvolver(object):
    "Evolves sequences on a tree."

    def __init__(self,
     seq_model=None,
     mutation_rate=None,
     seq_attr='sequences',
     seq_model_attr="seq_model",
     edge_length_attr="length",
     edge_rate_attr="mutation_rate",
     seq_label_attr='taxon'):
        "__init__ sets up meta-data dealing with object nomenclature and semantics."
        self.seq_model = seq_model
        self.mutation_rate = mutation_rate
        self.seq_attr = seq_attr
        self.seq_model_attr = seq_model_attr
        self.edge_length_attr = edge_length_attr
        self.edge_rate_attr = edge_rate_attr
        self.seq_label_attr = seq_label_attr

    def evolve_states(self,
     tree,
     seq_len,
     root_states=None,
     simulate_root_states=True,
     in_place=False,
     rng=None):
        """
        Appends a new sequence of length `seq_len` to a list at each node
        in `tree`.  The attribute name of this list in each node is given
        by `seq_attr`. If `seq_model` is None, `tree.seq_model` or
        `seq_model` at each node must be specified. If `in_place` is
        False, the tree is copied first, otherwise original tree is modified.
        If `root_states` is given, this will be used as the sequence for the root.
        If not, and if `simulate_root_states` is True, then the sequence for the
        root will be drawn from the stationary distribution of the character model.
        """
        if rng is None:
            rng = GLOBAL_RNG
        if not in_place:
            tree = copy.deepcopy(tree)

        if self.seq_model is None:
            seq_model = getattr(tree, self.seq_model_attr, None)

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
                seq_model  = getattr(edge, self.seq_model_attr, None) or self.seq_model
                length = getattr(edge, self.edge_length_attr)
                mutation_rate = getattr(edge, self.edge_rate_attr, None) or self.mutation_rate
                seq_list.append(seq_model.simulate_descendant_states(par_seq, length, mutation_rate))
            else:
                # no tail node: root
                n_prev_seq = len(seq_list)
                if root_states is not None:
                    seq_list.append(root_states)
                elif simulate_root_states:
                    seq_model  = getattr(node.edge, self.seq_model_attr, None) or self.seq_model
                    seq_list.append(seq_model.stationary_sample(seq_len, rng=rng))
                else:
                    assert n_prev_seq > 0
                    n_prev_seq -= 1
        return tree

    def compose_char_map(self, tree, taxon_set=None, include=None, exclude=None):
        """
        Returns a CharacterDataMap where the keys are the taxa of the leaf_nodes
        of `source_tree` and the values are sequences, with each sequence being the
        concatentation of all the sequences in the list of sequences associated with
        each tip. Specific sequences to be included/excluded can be fine-tuned using
        the `include` and `exclude` args, where `include`=None means to include all
        by default, and `exclude`=None means to exclude all by default.
        """
        char_map = dendropy.CharacterDataMap()
        for leaf in tree.leaf_nodes():
            cvec = dendropy.CharacterDataSequence(taxon=leaf.taxon)
            seq_list = getattr(leaf, self.seq_attr)
            for seq_idx, seq in enumerate(seq_list):
                if ((include is None) or (seq_idx in include))  \
                    and ((exclude is None) or (seq_idx not in exclude)):
                    for state in seq:
                        cvec.append(dendropy.CharacterDataCell(value=state))
            if taxon_set is not None:
                taxon = taxon_set.require_taxon(label=leaf.taxon.label)
            else:
                taxon = leaf.taxon
            char_map[taxon] = cvec
        return char_map

##############################################################################
## Wrappers for Convenience

def simulate_discrete_char_dataset(seq_len,
        tree_model,
        seq_model,
        mutation_rate=1.0,
        root_states=None,
        dataset=None,
        rng=None):
    """
    Wrapper to conveniently generate a DataSet simulated under
    the given tree and character model.
    `seq_len`       : length of sequence (number of characters)
    `tree_model`    : dendropy.Tree object
    `seq_model`     : dendropy.model.discrete.DiscreteCharacterEvolutionModel object
    `mutation_rate` : mutation *modifier* rate (should be 1.0 if branch lengths
                      on tree reflect true expected number of changes
    `root_states`   : vector of root states (length must equal `seq_len`)
    `dataset`       : a dendropy.DataSet object.
                      if given, the new
                      dendropy.CharacterMatrix object will be added to
                      this (along with a new taxon_namespace if required). Otherwise,
                      a new dendropy.DataSet object will be created.
    `rng`           : random number generator; if not given, `GLOBAL_RNG` will be
                      used
    Returns: a dendropy.DataSet object object.
    """
    if dataset is None:
        dataset = dendropy.DataSet()
    if tree_model.taxon_namespace not in dataset.taxon_namespaces:
        taxon_namespace = dataset.add_taxon_namespace(tree_model.taxon_namespace)
    else:
        taxon_namespace = tree_model.taxon_namespace
    char_matrix = simulate_discrete_char_matrix(
        seq_len=seq_len,
        tree_model=tree_model,
        seq_model=seq_model,
        mutation_rate=mutation_rate,
        root_states=root_states,
        char_matrix=None,
        rng=None)
    dataset.add_char_matrix(char_matrix=char_matrix)
    return dataset

def simulate_discrete_char_matrix(
        seq_len,
        tree_model,
        seq_model,
        mutation_rate=1.0,
        root_states=None,
        char_matrix=None,
        rng=None):
    """
    Wrapper to conveniently generate a characters simulated under
    the given tree and character model.
    `seq_len`       : length of sequence (number of characters)
    `tree_model`    : dendropy.Tree object
    `seq_model`    : dendropy.model.discrete.DiscreteCharacterEvolutionModel object
    `mutation_rate` : mutation *modifier* rate (should be 1.0 if branch lengths
                      on tree reflect true expected number of changes
    `root_states`   : vector of root states (length must equal `seq_len`)
    `char_matrix`    : dendropy.CharacterMatrix object.
                      if given, new sequences for taxa on `tree_model` leaf_nodes
                      will be appended to existing sequences of corresponding
                      taxa in char_matrix; if not, a new
                      dendropy.CharacterMatrix object will be created
    `rng`           : random number generator; if not given, `GLOBAL_RNG` will be
                      used

    Returns: a dendropy.CharacterMatrix object.

    Since characters will be appended to existing sequences, you can simulate a
    sequences under a mixed model by calling this method multiple times with
    different character models and/or different mutation rates, passing
    in the same `char_matrix` object each time.
    """
    seq_evolver = DiscreteCharacterEvolver(seq_model=seq_model,
                               mutation_rate=mutation_rate)
    tree = seq_evolver.evolve_states(tree=tree_model,
        seq_len=seq_len,
        root_states=None,
        rng=rng)
    char_map = seq_evolver.compose_char_map(tree, tree.taxon_namespace)
    if char_matrix is None:
        char_matrix = dendropy.DnaCharacterMatrix()
        char_matrix.taxon_namespace = tree_model.taxon_namespace
    if char_matrix.taxon_namespace is None:
        char_matrix.taxon_namespace = tree_model.taxon_namespace
    else:
        assert char_matrix.taxon_namespace is tree_model.taxon_namespace, "conflicting taxon sets"
    char_matrix.extend_map(other_map=char_map,
        overwrite_existing=False,
        extend_existing=True)
    return char_matrix

