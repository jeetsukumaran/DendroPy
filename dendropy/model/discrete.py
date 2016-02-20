#! /usr/bin/env python

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
Models and modeling of discrete character evolution.
"""

import copy
import math
import itertools
from dendropy.utility import GLOBAL_RNG
from dendropy.calculate import probability
import dendropy

############################################################################
## Character Evolution Modeling

class DiscreteCharacterEvolutionModel(object):
    "Base class for discrete character substitution models."

    def __init__(self, state_alphabet, stationary_freqs=None, rng=None):
        """
        __init__ initializes the state_alphabet to define the character type on which
        this model acts.  The objects random number generator will be ``rng`` or 'GLOBAL_RNG'
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
            anc_state_idx = state.index
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
            in_place=True,
            rng=None):
        """
        Appends a new sequence of length ``seq_len`` to a list at each node
        in ``tree``.  The attribute name of this list in each node is given
        by ``seq_attr``. If ``seq_model`` is None, ``tree.seq_model`` or
        ``seq_model`` at each node must be specified. If ``in_place`` is
        False, the tree is copied first, otherwise original tree is modified.
        If ``root_states`` is given, this will be used as the sequence for the root.
        If not, and if ``simulate_root_states`` is True, then the sequence for the
        root will be drawn from the stationary distribution of the character model.
        """
        if rng is None:
            rng = GLOBAL_RNG
        if not in_place:
            tree = tree.clone(1) # ==> taxon_namespace_scoped_copy()

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

    def extend_char_matrix_with_characters_on_tree(self,
            char_matrix,
            tree,
            include=None,
            exclude=None):
        """
        Creates a character matrix with new sequences (or extends sequences of
        an existing character matrix if provided via ``char_matrix``),
        where the the sequence for each taxon corresponds to the concatenation
        of all sequences in the list of sequences associated with tip that
        references the given taxon.
        Specific sequences to be included/excluded can be fine-tuned using the
        ``include`` and ``exclude`` args, where ``include=None`` means to include all
        by default, and ``exclude=None`` means to exclude all by default.
        """
        for leaf in tree.leaf_nodes():
            cvec = char_matrix[leaf.taxon]
            seq_list = getattr(leaf, self.seq_attr)
            for seq_idx, seq in enumerate(seq_list):
                if ((include is None) or (seq_idx in include))  \
                    and ((exclude is None) or (seq_idx not in exclude)):
                    for state in seq:
                        cvec.append(state)
        return char_matrix

    def clean_tree(self, tree):
        for nd in tree:
            # setattr(nd, self.seq_attr, [])
            delattr(nd, self.seq_attr)

############################################################################
## Specialized Models: nucldeotides

class NucleotideCharacterEvolutionModel(DiscreteCharacterEvolutionModel):
    "General nucleotide substitution model."

    def __init__(self, base_freqs=None, state_alphabet=None, rng=None):
        "__init__ calls SeqModel.__init__ and sets the base_freqs field"
        if state_alphabet is None:
            state_alphabet = dendropy.DNA_STATE_ALPHABET
        DiscreteCharacterEvolutionModel.__init__(
                self,
                state_alphabet=state_alphabet,
                rng=rng)
        if base_freqs is None:
            self.base_freqs = [0.25, 0.25, 0.25, 0.25]
        else:
            self.base_freqs = base_freqs

    def stationary_sample(self, seq_len, rng=None):
        """
        Returns a NucleotideSequence() object with length ``length``
        representing a sample of characters drawn from this model's
        stationary distribution.
        """
        probs = self.base_freqs
        char_state_indices = [probability.sample_multinomial(probs, rng) for i in range(seq_len)]
        return [self.state_alphabet[idx] for idx in char_state_indices]

    def is_purine(self, state_index):
        """
        Returns True if state_index represents a purine (A or G) row or column
        index: 0, 2
        """
        return state_index % 2 == 0

    def is_pyrimidine(self, state_index):
        """
        Returns True if state_index represents a pyrimidine (C or T) row or column
        index: 1, 3
        """
        return state_index % 2 == 1

    def is_transversion(self, state1_idx, state2_idx):
        """
        Returns True if the change from state1 to state2, as
        represented by the row or column indices, is a transversional
        change.
        """
        return (self.is_purine(state1_idx) and self.is_pyrimidine(state2_idx)) \
               or (self.is_pyrimidine(state1_idx) and self.is_purine(state2_idx))

    def is_purine_transition(self, state1_idx, state2_idx):
        """
        Returns True if the change from state1 to state2, as
        represented by the row or column indices, is a purine
        transitional change.
        """
        return self.is_purine(state1_idx) and self.is_purine(state2_idx)

    def is_pyrimidine_transition(self, state1_idx, state2_idx):
        """
        Returns True if the change from state1 to state2, as
        represented by the row or column indices, is a pyrimidine
        transitional change.
        """
        return self.is_pyrimidine(state1_idx) \
               and self.is_pyrimidine(state2_idx)

    def is_transition(self, state1_idx, state2_idx):
        """
        Returns True if the change from state1 to state2, as
        represented by the row or column indices, is a
        transitional change.
        """
        return (self.is_purine(state1_idx) and self.is_purine(state2_idx)) \
               or (self.is_pyrimidine(state1_idx) and self.is_pyrimidine(state2_idx))

class Hky85(NucleotideCharacterEvolutionModel):
    """
    Hasegawa et al. 1985 model. Implementation following Swofford et
    al., 1996.
    """

    def __init__(self, kappa=1.0, base_freqs=None, state_alphabet=None, rng=None):
        "__init__: if no arguments given, defaults to JC69."
        if state_alphabet is None:
            state_alphabet = dendropy.DNA_STATE_ALPHABET
        NucleotideCharacterEvolutionModel.__init__(
                self,
                base_freqs=base_freqs,
                state_alphabet=state_alphabet,
                rng=rng)
        self.correct_rate = True
        self.kappa = kappa

    def __repr__(self):
        rep = "kappa=%f bases=%s" % (self.kappa, str(self.base_freqs))
        return rep

    def corrected_substitution_rate(self, rate):
        """Returns the factor that we have to multiply to the branch length
        to make branch lengths proportional to # of substitutions per site."""
        if self.correct_rate:
            pia = self.base_freqs[0]
            pic = self.base_freqs[1]
            pig = self.base_freqs[2]
            pit = self.base_freqs[3]
            f = self.kappa*(pia*pig + pic*pit)
            f += (pia + pig)*(pic + pit)
            return (rate * 0.5/f)  # (rate * 0.5/f)
        else:
            return rate

    def pij(self, state_i, state_j, tlen, rate=1.0):
        """
        Returns probability, p_ij, of going from state i to state j
        over time tlen at given rate. (tlen * rate = nu, expected
        number of substitutions)
        """
        nu = self.corrected_substitution_rate(rate) * tlen
        if self.is_purine(state_j):
            sumfreqs = self.base_freqs[0] + self.base_freqs[2]
        else:
            sumfreqs = self.base_freqs[1] + self.base_freqs[3]
        factorA = 1 + (sumfreqs * (self.kappa - 1.0))
        if state_i == state_j:
            pij = self.base_freqs[state_j] \
                  + self.base_freqs[state_j] \
                      * (1.0/sumfreqs - 1) * math.exp(-1.0 * nu) \
                  + ((sumfreqs - self.base_freqs[state_j])/sumfreqs) \
                      * math.exp(-1.0 * nu * factorA)

        elif self.is_transition(state_i, state_j):
            pij = self.base_freqs[state_j] \
                  + self.base_freqs[state_j] \
                      * (1.0/sumfreqs - 1) * math.exp(-1.0 * nu) \
                  - (self.base_freqs[state_j] / sumfreqs) \
                      * math.exp(-1.0 * nu * factorA)
        else:
            pij = self.base_freqs[state_j] * (1.0 - math.exp(-1.0 * nu))
        return pij

    def qmatrix(self, rate=1.0):
        "Returns the instantaneous rate of change matrix."
        rate = self.corrected_substitution_rate(rate)
        qmatrix = []
        for state_i in range(4):
            qmatrix.append([])
            for state_j in range(4):
                if state_i == state_j:
                    # we cheat here and insert a placeholder till the
                    # other cells are calculated
                    qij = 0.0
                else:
                    if self.is_transition(state_i, state_j):
                        qij = rate * self.kappa * self.base_freqs[state_j]
                    else:
                        qij = rate * self.base_freqs[state_j]
                qmatrix[state_i].append(qij)
        for state in range(4):
            qmatrix[state][state] = -1.0 * sum(qmatrix[state])
        return qmatrix

    def pvector(self, state, tlen, rate=1.0):
        """
        Returns a vector of transition probabilities for a given state
        over time ``tlen`` at rate ``rate`` for ``state``. (tlen * rate =
        nu, expected number of substitutions)
        """
        pvec = []
        # in case later we want to allow characters passed in here
        state_i = state
        for state_j in range(4):
            pvec.append(self.pij(state_i, state_j, tlen=tlen, rate=rate))
        return pvec

    def pmatrix(self, tlen, rate=1.0):
        """
        Returns a matrix of nucleotide substitution
        probabilities. Based on analytical solution by Swofford et
        al., 1996. (tlen * rate = nu, expected number of
        substitutions)
        """
        pmatrix = []
        for state_i in range(4):
            pmatrix.append(self.pvector(state_i, tlen=tlen, rate=rate))
        return pmatrix

class Jc69(Hky85):
    """
    Jukes-Cantor 1969 model. Specializes HKY85 such that
    kappa = 1.0, and base frequencies = [0.25, 0.25, 0.25, 0.25].
    """
    def __init__(self, state_alphabet=None, rng=None):
        "__init__: uses Hky85.__init__"
        Hky85.__init__(self,
                kappa=1.0,
                base_freqs=[0.25, 0.25, 0.25, 0.25],
                state_alphabet=state_alphabet,
                rng=rng,
                )



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

    Parameters
    ----------

    seq_len       : int
        Length of sequence (number of characters).
    tree_model    : |Tree|
        Tree on which to simulate.
    seq_model     : dendropy.model.discrete.DiscreteCharacterEvolutionModel
        The character substitution model under which to to evolve the
        characters.
    mutation_rate : float
        Mutation *modifier* rate (should be 1.0 if branch lengths on tree
        reflect true expected number of changes).
    root_states``   : list
        Vector of root states (length must equal ``seq_len``).
    dataset       : |DataSet|
        If given, the new dendropy.CharacterMatrix object will be
        added to this (along with a new taxon_namespace if
        required). Otherwise, a new dendropy.DataSet
        object will be created.
    rng           : random number generator
        If not given, 'GLOBAL_RNG' will be used.

    Returns
    -------

    d : |DataSet|

    """
    if dataset is None:
        dataset = dendropy.DataSet()
    if tree_model.taxon_namespace not in dataset.taxon_namespaces:
        taxon_namespace = dataset.add_taxon_namespace(tree_model.taxon_namespace)
    else:
        taxon_namespace = tree_model.taxon_namespace
    char_matrix = simulate_discrete_chars(
        seq_len=seq_len,
        tree_model=tree_model,
        seq_model=seq_model,
        mutation_rate=mutation_rate,
        root_states=root_states,
        char_matrix=None,
        rng=None)
    dataset.add_char_matrix(char_matrix=char_matrix)
    return dataset

def simulate_discrete_chars(
        seq_len,
        tree_model,
        seq_model,
        mutation_rate=1.0,
        root_states=None,
        char_matrix=None,
        retain_sequences_on_tree=False,
        rng=None):
    """
    Wrapper to conveniently generate a characters simulated under
    the given tree and character model.

    Since characters will be appended to existing sequences, you can simulate a
    sequences under a mixed model by calling this method multiple times with
    different character models and/or different mutation rates, passing
    in the same ``char_matrix`` object each time.

    Parameters
    ----------

    seq_len       : int
        Length of sequence (number of characters).
    tree_model    : |Tree|
        Tree on which to simulate.
    seq_model     : dendropy.model.discrete.DiscreteCharacterEvolutionModel
        The character substitution model under which to to evolve the
        characters.
    mutation_rate : float
        Mutation *modifier* rate (should be 1.0 if branch lengths on tree
        reflect true expected number of changes).
    root_states``   : list
        Vector of root states (length must equal ``seq_len``).
    char_matrix   : |DnaCharacterMatrix|
        If given, new sequences for taxa on ``tree_model`` leaf_nodes will be
        appended to existing sequences of corresponding taxa in char_matrix; if
        not, a new |DnaCharacterMatrix| object will be created.
    retain_sequences_on_tree : bool
        If |False|, sequence annotations will be cleared from tree after
        simulation. Set to |True| if you want to, e.g., evolve and accumulate
        different sequences on tree, or retain information for other purposes.
    rng           : random number generator
        If not given, 'GLOBAL_RNG' will be used.

    Returns
    -------
    d : a dendropy.datamodel.CharacterMatrix object.

    """
    seq_evolver = DiscreteCharacterEvolver(seq_model=seq_model,
                               mutation_rate=mutation_rate)
    tree = seq_evolver.evolve_states(
        tree=tree_model,
        seq_len=seq_len,
        root_states=None,
        rng=rng)
    if char_matrix is None:
        char_matrix = dendropy.DnaCharacterMatrix(taxon_namespace=tree_model.taxon_namespace)
        char_matrix.taxon_namespace = tree_model.taxon_namespace
    else:
        assert char_matrix.taxon_namespace is tree_model.taxon_namespace, "conflicting taxon sets"
    seq_evolver.extend_char_matrix_with_characters_on_tree(
            char_matrix=char_matrix,
            tree=tree)
    if not retain_sequences_on_tree:
        seq_evolver.clean_tree(tree)
    return char_matrix

def hky85_chars(
        seq_len,
        tree_model,
        mutation_rate=1.0,
        kappa=1.0,
        base_freqs=[0.25, 0.25, 0.25, 0.25],
        root_states=None,
        char_matrix=None,
        retain_sequences_on_tree=False,
        rng=None):
    """
    Convenience class to wrap generation of characters (as a CharacterBlock
    object) based on the HKY model.

    Parameters
    ----------

    seq_len       : int
        Length of sequence (number of characters).
    tree_model    : |Tree|
        Tree on which to simulate.
    mutation_rate : float
        Mutation *modifier* rate (should be 1.0 if branch lengths on tree
        reflect true expected number of changes).
    root_states``   : list
        Vector of root states (length must equal ``seq_len``).
    char_matrix   : |DnaCharacterMatrix|
        If given, new sequences for taxa on ``tree_model`` leaf_nodes will be
        appended to existing sequences of corresponding taxa in char_matrix; if
        not, a new |DnaCharacterMatrix| object will be created.
    retain_sequences_on_tree : bool
        If |False|, sequence annotations will be cleared from tree after
        simulation. Set to |True| if you want to, e.g., evolve and accumulate
        different sequences on tree, or retain information for other purposes.
    rng           : random number generator
        If not given, 'GLOBAL_RNG' will be used.

    Returns
    -------
    d : |DnaCharacterMatrix|
        The simulated alignment.

    Since characters will be appended to existing sequences, you can simulate a
    sequences under a mixed model by calling this method multiple times with
    different character model parameter values and/or different mutation
    rates, passing in the same ``char_matrix`` object each time.
    """
    if char_matrix is None:
        char_matrix = dendropy.DnaCharacterMatrix(taxon_namespace=tree_model.taxon_namespace)
    else:
        assert char_matrix.taxon_namespace is tree_model.taxon_namespace
    state_alphabet = char_matrix.default_state_alphabet
    seq_model = Hky85(
            kappa=kappa,
            base_freqs=base_freqs,
            state_alphabet=state_alphabet)
    return simulate_discrete_chars(seq_len=seq_len,
                               tree_model=tree_model,
                               seq_model=seq_model,
                               mutation_rate=mutation_rate,
                               root_states=root_states,
                               char_matrix=char_matrix,
                               rng=rng)

