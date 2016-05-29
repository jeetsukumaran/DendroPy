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
Population genetic statistics.
"""

import math
import dendropy
from dendropy.calculate import probability
from dendropy.calculate import combinatorics

###############################################################################
## internal functions: generally taking lower-level data, such as sequences etc.
###############################################################################

def _count_differences(char_sequences, state_alphabet, ignore_uncertain=True):
    """
    Returns pair of values: total number of pairwise differences observed between
    all sequences, and mean number of pairwise differences pair base.
    """
    sum_diff = 0.0
    mean_diff = 0.0
    sq_diff = 0.0
    comps = 0

    #Check that all sequences are the same length
    if len(set([len(seq) for seq in char_sequences])) != 1:
        raise Exception("sequences of unequal length")

    if ignore_uncertain:
        attr = "fundamental_indexes_with_gaps_as_missing"
        _states_to_ignore = [state_alphabet.gap_state, state_alphabet.no_data_state]
        states_to_ignore = set([getattr(char, attr) for char in _states_to_ignore])
    else:
        attr = "fundamental_indexes"
        states_to_ignore = set()

    reduced_char_sequences = []
    for sequence in char_sequences:
        seq = [getattr(char, attr) for char in sequence]
        reduced_char_sequences.append(seq)

    for vidx, i in enumerate(reduced_char_sequences[:-1]):
        for j in reduced_char_sequences[vidx+1:]:
            diff = 0
            counted = 0
            comps += 1
            for cidx, c in enumerate(i):
                c1 = c
                c2 = j[cidx]
                if c1 in states_to_ignore or c2 in states_to_ignore:
                    continue
                counted += 1
                if c1 is not c2:
                    diff += 1
            sum_diff += float(diff)
            # If counted < 0, this means that there is sites between these sequences
            # in which both are not ignored: i.e., one or the other has a gap
            # or an uncertain character. We consider this to mean (maybe
            # somewhat paradoxically) that there are no sites that are
            # different between the sequences. Put less paradoxically: there
            # are no non-ignored sites that are different between the
            # sequences.
            mean_diff += (float(diff) / counted) if counted > 0 else float(diff)
            sq_diff += (diff ** 2)
    return sum_diff, mean_diff / comps, sq_diff

def _nucleotide_diversity(char_sequences, state_alphabet, ignore_uncertain=True):
    """
    Returns $\pi$, the proportional nucleotide diversity, calculated for a
    list of character sequences.
    """
    return _count_differences(char_sequences, state_alphabet, ignore_uncertain)[1]

def _average_number_of_pairwise_differences(char_sequences, state_alphabet, ignore_uncertain=True):
    """
    Returns $k$ (Tajima 1983; Wakely 1996), calculated for a set of sequences:

    k = \frac{\right(\sum \sum \k_{ij}\left)}{n \choose 2}

    where $k_{ij}$ is the number of pairwise differences between the
    $i$th and $j$th sequence, and $n$ is the number of DNA sequences
    sampled.
    """
    sum_diff, mean_diff, sq_diff = _count_differences(char_sequences, state_alphabet, ignore_uncertain)
    return sum_diff / combinatorics.choose(len(char_sequences), 2)

def _num_segregating_sites(char_sequences, state_alphabet, ignore_uncertain=True):
    """
    Returns the raw number of segregating sites (polymorphic sites).
    """
    s = 0
    if ignore_uncertain:
        attr = "fundamental_indexes_with_gaps_as_missing"
        _states_to_ignore = [state_alphabet.gap_state, state_alphabet.no_data_state]
        states_to_ignore= set([getattr(char, attr) for char in _states_to_ignore])
    else:
        attr = "fundamental_indexes"
        states_to_ignore = set()
    for i, c1 in enumerate(char_sequences[0]):
        for v in char_sequences[1:]:
            c2 = v[i]
            f1 = getattr(c1, attr)
            f2 = getattr(c2, attr)
            if f1 in states_to_ignore or f2 in states_to_ignore:
                continue
            if f1 is not f2:
                s += 1
                break
    return s

def _tajimas_d(num_sequences, avg_num_pairwise_differences, num_segregating_sites):

    ### VERIFICATION ###
    ###
    ### Given: num_sequences = 10, num_pairwise_differences = 3.888889, num_segregating_sites = 16
    ###  i.e.: tajimas_d(10, 3.888889, 16)  == -1.44617198561
    ###  Then:    a1 == 2.82896825397
    ###           a2 == 1.53976773117
    ###           b1 == 0.407407407407
    ###           b2 == 0.279012345679
    ###           c1 == 0.0539216450284
    ###           c2 == 0.0472267720013
    ###           e1 == 0.0190605338016
    ###           e2 == 0.0049489277699
    ###           D ==  -1.44617198561

    a1 = sum([1.0/i for i in range(1, num_sequences)])
    a2 = sum([1.0/(i**2) for i in range(1, num_sequences)])
    b1 = float(num_sequences+1)/(3*(num_sequences-1))
    b2 = float(2 * ( (num_sequences**2) + num_sequences + 3 )) / (9*num_sequences*(num_sequences-1))
    c1 = b1 - 1.0/a1
    c2 = b2 - float(num_sequences+2)/(a1 * num_sequences) + float(a2)/(a1 ** 2)
    e1 = float(c1) / a1
    e2 = float(c2) / ( (a1**2) + a2 )
    D = (
        float(avg_num_pairwise_differences - (float(num_segregating_sites)/a1))
        / math.sqrt(
            (e1 * num_segregating_sites )
          + ((e2 * num_segregating_sites) * (num_segregating_sites - 1) ))
        )
    return D

###############################################################################
## friendlier-functions, generally taking a CharacterMatrix
###############################################################################

def num_segregating_sites(char_matrix, ignore_uncertain=True):
    """
    Returns the raw number of segregating sites (polymorphic sites).
    """
    return _num_segregating_sites(
            char_matrix.sequences(),
            char_matrix.default_state_alphabet,
            ignore_uncertain)

def average_number_of_pairwise_differences(char_matrix, ignore_uncertain=True):
    """
    Returns $k$, calculated for a character block.
    """
    return _average_number_of_pairwise_differences(char_matrix.sequences(), char_matrix.default_state_alphabet, ignore_uncertain)

def nucleotide_diversity(char_matrix, ignore_uncertain=True):
    """
    Returns $\pi$, calculated for a character block.
    """
    return _nucleotide_diversity(char_matrix.sequences(), char_matrix.default_state_alphabet, ignore_uncertain)

def tajimas_d(char_matrix, ignore_uncertain=True):
    """
    Returns Tajima's D.
    """
    sequences = char_matrix.sequences()
    num_sequences = len(sequences)
    avg_num_pairwise_differences = _average_number_of_pairwise_differences(sequences, char_matrix.default_state_alphabet, ignore_uncertain=ignore_uncertain)
    num_segregating_sites = _num_segregating_sites(
            sequences,
            char_matrix.default_state_alphabet,
            ignore_uncertain=ignore_uncertain)
    return _tajimas_d(num_sequences, avg_num_pairwise_differences, num_segregating_sites)

def wattersons_theta(char_matrix, ignore_uncertain=True):
    """
    Returns Watterson's Theta (per sequence)
    """
    sequences = char_matrix.sequences()
    num_segregating_sites = _num_segregating_sites(
            sequences,
            char_matrix.default_state_alphabet,
            ignore_uncertain=ignore_uncertain)
    a1 = sum([1.0/i for i in range(1, len(sequences))])
    return float(num_segregating_sites) / a1

###############################################################################
## Classes
###############################################################################

class PopulationPairSummaryStatistics(object):

    def __init__(self, pop1_seqs, pop2_seqs, ignore_uncertain=True):
        self.pop1_seqs = pop1_seqs
        self.pop2_seqs = pop2_seqs
        self.combined_seqs = pop1_seqs + pop2_seqs
        self.ignore_uncertain = ignore_uncertain
        self.state_alphabet = dendropy.DNA_STATE_ALPHABET

        self.average_number_of_pairwise_differences = 0
        self.average_number_of_pairwise_differences_between = 0
        self.average_number_of_pairwise_differences_within = 0
        self.average_number_of_pairwise_differences_net = 0
        self.num_segregating_sites = 0
        self.wattersons_theta = 0.0
        self.wakeleys_psi = 0.0
        self.tajimas_d = 0.0
        if self.ignore_uncertain:
            self.state_attr = "fundamental_indexes_with_gaps_as_missing"
            self.states_to_ignore = set([self.state_alphabet.gap_state, self.state_alphabet.no_data_state])
        else:
            self.state_attr = "fundamental_indexes"
            self.states_to_ignore = set()
        self.calc()

    def calc(self):
        """
        Returns a summary of a set of sequences that can be partitioned into
        the list of lists of taxa given by ``taxon_groups``.
        """
        diffs_x, mean_diffs_x, sq_diff_x = _count_differences(self.pop1_seqs, self.state_alphabet, self.ignore_uncertain)
        diffs_y, mean_diffs_y, sq_diff_y = _count_differences(self.pop2_seqs, self.state_alphabet, self.ignore_uncertain)
        d_x = diffs_x / combinatorics.choose(len(self.pop1_seqs), 2)
        d_y = diffs_y / combinatorics.choose(len(self.pop2_seqs), 2)
        d_xy = self._average_number_of_pairwise_differences_between_populations()
        s2_x = (float(sq_diff_x) / combinatorics.choose(len(self.pop1_seqs), 2) ) - (d_x ** 2)
        s2_y = (float(sq_diff_y) / combinatorics.choose(len(self.pop2_seqs), 2) ) - (d_y ** 2)
        s2_xy = self._variance_of_pairwise_differences_between_populations(d_xy)
        n = len(self.combined_seqs)
        n_x = float(len(self.pop1_seqs))
        n_y = float(len(self.pop2_seqs))
        a = float(n * (n-1))
        ax = float(n_x * (n_x - 1))
        ay = float(n_y * (n_y - 1))
        k = _average_number_of_pairwise_differences(self.combined_seqs, self.state_alphabet, self.ignore_uncertain)
        n = len(self.combined_seqs)

        # Hickerson 2006: pi #
        self.average_number_of_pairwise_differences = k

        # Hickerson 2006: pi_b #
        self.average_number_of_pairwise_differences_between = d_xy

        # Hickerson 2006: pi_w #
        self.average_number_of_pairwise_differences_within = d_x + d_y

        # Hickerson 2006: pi_net #
        self.average_number_of_pairwise_differences_net = d_xy - (d_x + d_y)

        # Hickerson 2006: S #
        self.num_segregating_sites = _num_segregating_sites(
                self.combined_seqs,
                self.state_alphabet,
                self.ignore_uncertain)

        # Hickerson 2006: theta #
        a1 = sum([1.0/i for i in range(1, n)])
        self.wattersons_theta = float(self.num_segregating_sites) / a1

        # Wakeley 1996 #
        self.wakeleys_psi = (float(1)/(a)) * ( ax * (math.sqrt(s2_x)/d_x) + ay * (math.sqrt(s2_y)/d_y) + (2 * n_x * n_y * math.sqrt(s2_xy)/k))

        # Tajima's D #
        self.tajimas_d = _tajimas_d(n, self.average_number_of_pairwise_differences, self.num_segregating_sites)

    def _average_number_of_pairwise_differences_between_populations(self):
        """
        Implements Eq (3) of:

        Wakeley, J. 1996. Distinguishing migration from isolation using the
        variance of pairwise differences. Theoretical Population Biology 49:
        369-386.
        """
        diffs = 0
        for sx in self.pop1_seqs:
            for sy in self.pop2_seqs:
                for cidx, c in enumerate(sx):
                    c1 = c
                    c2 = sy[cidx]
                    if c1 in self.states_to_ignore or c2 in self.states_to_ignore:
                        continue
                    f1 = getattr(c1, self.state_attr)
                    f2 = getattr(c2, self.state_attr)
                    if f1 != f2:
                        diffs += 1
        dxy = float(1)/(len(self.pop1_seqs) * len(self.pop2_seqs)) * float(diffs)
        return dxy

    def _variance_of_pairwise_differences_between_populations(self, mean_diff):
        """
        Implements Eq (10) of:

        Wakeley, J. 1996. Distinguishing migration from isolation using the
        variance of pairwise differences. Theoretical Population Biology 49:
        369-386.
        """
        ss_diffs = 0
        for sx in self.pop1_seqs:
            for sy in self.pop2_seqs:
                diffs = 0
                for cidx, c in enumerate(sx):
                    c1 = c
                    c2 = sy[cidx]
                    if c1 in self.states_to_ignore or c2 in self.states_to_ignore:
                        continue
                    f1 = getattr(c1, self.state_attr)
                    f2 = getattr(c2, self.state_attr)
                    if f1 != f2:
                        diffs += 1
                ss_diffs += (float(diffs - mean_diff) ** 2)
        return float(ss_diffs)/(len(self.pop1_seqs)*len(self.pop2_seqs))

def derived_state_matrix(
        char_matrix,
        ancestral_sequence=None,
        derived_state_alphabet=None,
        ignore_uncertain=True,
        ):
    """
    Given a list of CharDataSequence objects, and a reference ancestral sequence,
    this returns a list of strings corresponding to the list of CharDataSequence
    objects, where a '0' indicates the ancestral state and '1' a derived state.

    e.g.

        Given:
                GGCTAATCTGA
                GCTTTTTCTGA
                GCTCTCTCTTC

        with ancestral sequence:
                GGTTAATCTGA

        this returns:
                0010000000
                0000110000
                0001110011
    """
    if derived_state_alphabet is None:
        derived_state_alphabet = dendropy.StateAlphabet(
                fundamental_states="01",
                polymorphic_states=None,
                ambiguous_states=None,
                no_data_symbol="?",
                gap_symbol="-")
    derived_matrix = dendropy.StandardCharacterMatrix(
            taxon_namespace=char_matrix.taxon_namespace,
            default_state_alphabet=derived_state_alphabet)
    if ignore_uncertain:
        attr = "fundamental_indexes_with_gaps_as_missing"
        states_to_ignore = set([char_matrix.default_state_alphabet.gap_state, char_matrix.default_state_alphabet.no_data_state])
    else:
        attr = "fundamental_indexes"
        states_to_ignore = set()
    if ancestral_sequence is None:
        ancestral_sequence = char_matrix[0]
    ancestral_fundamental_ids = []
    for idx, c1 in enumerate(ancestral_sequence):
        if c1 in states_to_ignore:
            ancestral_fundamental_ids.append(None)
        else:
            ancestral_fundamental_ids.append(getattr(c1, attr))
    for taxon in char_matrix:
        s1 =  char_matrix[taxon]
        for idx, c2 in enumerate(s1):
            if ancestral_fundamental_ids[idx] is None or c2 in states_to_ignore:
                derived_matrix[taxon].append(derived_matrix.default_state_alphabet["?"])
                continue
            f2 = getattr(c2, attr)
            if f2 == ancestral_fundamental_ids[idx]:
                derived_matrix[taxon].append(derived_matrix.default_state_alphabet["0"])
            else:
                derived_matrix[taxon].append(derived_matrix.default_state_alphabet["1"])
    return derived_matrix

def unfolded_site_frequency_spectrum(
        char_matrix,
        ancestral_sequence=None,
        ignore_uncertain=False,
        pad=True):
    """
    Returns the site frequency spectrum of list of CharDataSequence objects given by char_sequences,
    with reference to the ancestral sequence given by ancestral_seq. If ancestral_seq
    is None, then the first sequence in char_sequences is taken to be the ancestral
    sequence.
    """
    dsm = derived_state_matrix(
            char_matrix=char_matrix,
            ancestral_sequence=ancestral_sequence,
            derived_state_alphabet=None,
            ignore_uncertain=ignore_uncertain,
            )
    sites = zip(*dsm.sequences()) # transpose
    freqs = {}
    if pad:
        for i in range(len(char_matrix)+1):
            freqs[i] = 0
    for s in sites:
        p = sum(1 for i in s if i.symbol == "1")
        if p not in freqs:
            freqs[p] = 1
        else:
            freqs[p] += 1
    return freqs
