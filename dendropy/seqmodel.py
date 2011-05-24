#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
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
Models of molecular character evolution.
"""

import math
import itertools

from dendropy.utility import GLOBAL_RNG
from dendropy.utility import probability
import dendropy

class SeqModel(object):
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

    def generate_descendant_states(self,
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

class NucleotideSeqModel(SeqModel):
    "General nucleotide substitution model."

    def __init__(self, base_freqs=None, state_alphabet=None):
        "__init__ calls SeqModel.__init__ and sets the base_freqs field"
        if state_alphabet is None:
            state_alphabet = dendropy.DNA_STATE_ALPHABET
        SeqModel.__init__(self, state_alphabet)
        if base_freqs is None:
            self.base_freqs = [0.25, 0.25, 0.25, 0.25]
        else:
            self.base_freqs = base_freqs

    def stationary_sample(self, seq_len, rng=None):
        """
        Returns a NucleotideSequence() object with length `length`
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

class Hky85SeqModel(NucleotideSeqModel):
    """
    Hasegawa et al. 1985 model. Implementation following Swofford et
    al., 1996.
    """

    def __init__(self, kappa=1.0, base_freqs=None):
        "__init__: if no arguments given, defaults to JC69."
        NucleotideSeqModel.__init__(self, base_freqs=base_freqs)
        self.correct_rate = True
        self.kappa = kappa
        if base_freqs is None:
            self.base_freqs = [0.25, 0.25, 0.25, 0.25]
        else:
            self.base_freqs = base_freqs

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
        over time `tlen` at rate `rate` for `state`. (tlen * rate =
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

class Jc69SeqModel(Hky85SeqModel):
    """
    Jukes-Cantor 1969 model. Specializes HKY85 such that
    kappa = 1.0, and base frequencies = [0.25, 0.25, 0.25, 0.25].
    """
    def __init__(self):
        "__init__: uses Hky85SeqModel.__init__"
        Hky85SeqModel.__init__(self,
                                     kappa=1.0,
                                     base_freqs=[0.25, 0.25, 0.25, 0.25])

