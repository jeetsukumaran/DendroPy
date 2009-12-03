#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
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
###############################################################################

"""
Population genetic statistics.
"""

import math
from dendropy.utility import probability

def _count_differences(char_vectors, state_alphabet, ignore_uncertain=True):
    """
    Returns pair of values: total number of pairwise differences observed between
    all sequences, and mean number of pairwise differences pair base.
    """
    sum_diff = 0.0
    mean_diff = 0.0
    sq_diff = 0.0
    total_counted = 0
    comps = 0
    for vidx, i in enumerate(char_vectors[:-1]):
        for j in char_vectors[vidx+1:]:
            if len(i) != len(j):
                raise Exception("sequences of unequal length")
            diff = 0
            counted = 0
            comps += 1
            for cidx, c in enumerate(i):
                c1 = c
                c2 = j[cidx]
                if (not ignore_uncertain) \
                    or (c1.value is not state_alphabet.gap \
                        and c2.value is not state_alphabet.gap \
                        and len(c1.value.fundamental_ids) == 1 \
                        and len(c2.value.fundamental_ids) == 1):
                    counted += 1
                    total_counted += 1
                    if c1.value is not c2.value:
                        diff += 1
            sum_diff += float(diff)
            mean_diff += float(diff) / counted
            sq_diff += (diff ** 2)
    return sum_diff, mean_diff / comps, sq_diff

def _nucleotide_diversity(char_vectors, state_alphabet, ignore_uncertain=True):
    """
    Returns $\pi$, the proportional nucleotide diversity, calculated for a
    list of character vectors.
    """
    return _count_differences(char_vectors, state_alphabet, ignore_uncertain)[1]

def _average_number_of_pairwise_differences(char_vectors, state_alphabet, ignore_uncertain=True):
    """
    Returns $k$ (Tajima 1983; Wakely 1996), calculated for a set of sequences:

    k = \frac{\right(\sum \sum \k_{ij}\left)}{n \choose 2}

    where $k_{ij}$ is the number of pairwise differences between the
    $i$th and $j$th sequence, and $n$ is the number of DNA sequences
    sampled.
    """
    sum_diff, mean_diff, sq_diff = _count_differences(char_vectors, state_alphabet, ignore_uncertain)
    return sum_diff / probability.binomial_coefficient(len(char_vectors), 2)

def _num_segregating_sites(char_vectors, state_alphabet, ignore_uncertain=True):
    """
    Returns the raw number of segregating sites (polymorphic sites).
    """
    s = 0
    for vidx, i in enumerate(char_vectors[:-1]):
        for j in char_vectors[vidx+1:]:
            if len(i) != len(j):
                raise Exception("sequences of unequal length")
            for cidx, c in enumerate(i):
                c1 = c
                c2 = j[cidx]
                if (not ignore_uncertain) \
                    or (c1.value is not state_alphabet.gap \
                        and c2.value is not state_alphabet.gap \
                        and len(c1.value.fundamental_ids) == 1 \
                        and len(c2.value.fundamental_ids) == 1):
                    if c1.value is not c2.value:
                        s += 1
                        continue
    return s

def num_segregating_sites(char_matrix, ignore_uncertain=True):
    """
    Returns the raw number of segregating sites (polymorphic sites).
    """
    return _num_segregating_sites(char_matrix.vectors(), char_matrix.default_state_alphabet, ignore_uncertain)

def average_number_of_pairwise_differences(char_matrix, ignore_uncertain=True):
    """
    Returns $k$, calculated for a character block.
    """
    return _average_number_of_pairwise_differences(char_matrix.vectors(), char_matrix.default_state_alphabet, ignore_uncertain)

def nucleotide_diversity(char_matrix, ignore_uncertain=True):
    """
    Returns $\pi$, calculated for a character block.
    """
    return _nucleotide_diversity(char_matrix.vectors(), char_matrix.default_state_alphabet, ignore_uncertain)

def _average_number_of_pairwise_differences_between_populations(char_x, char_y, state_alphabet, ignore_uncertain=True):
    """
    Implements Eq (3) of:

    Wakeley, J. 1996. Distinguishing migration from isolation using the
    variance of pairwise differences. Theoretical Population Biology 49:
    369-386.
    """
    diffs = 0
    for sx in char_x:
        for sy in char_y:
            for cidx, c in enumerate(sx):
                c1 = c
                c2 = sy[cidx]
                if (not ignore_uncertain) \
                    or (c1.value is not state_alphabet.gap \
                        and c2.value is not state_alphabet.gap \
                        and len(c1.value.fundamental_ids) == 1 \
                        and len(c2.value.fundamental_ids) == 1):
                    if c1.value is not c2.value:
                        diffs += 1
    dxy = float(1)/(len(char_x) * len(char_y)) * float(diffs)
    return dxy

def _variance_of_pairwise_differences_between_populations(char_x, char_y, mean_diff, state_alphabet, ignore_uncertain=True):
    """
    Implements Eq (10) of:

    Wakeley, J. 1996. Distinguishing migration from isolation using the
    variance of pairwise differences. Theoretical Population Biology 49:
    369-386.
    """
    ss_diffs = 0
    for sx in char_x:
        for sy in char_y:
            diffs = 0
            for cidx, c in enumerate(sx):
                c1 = c
                c2 = sy[cidx]
                if (not ignore_uncertain) \
                    or (c1.value is not state_alphabet.gap \
                        and c2.value is not state_alphabet.gap \
                        and len(c1.value.fundamental_ids) == 1 \
                        and len(c2.value.fundamental_ids) == 1):
                    if c1.value is not c2.value:
                        diffs += 1
            ss_diffs += (float(diffs - mean_diff) ** 2)
    return float(ss_diffs)/(len(char_x)*len(char_y))

def wakeleys_Psi(char_matrix, taxon_groups, ignore_uncertain=True):
    """
    Returns Wakeley's psi, as described in:

    Wakeley, J. 1996. Distinguishing migration from isolation using the
    variance of pairwise differences. Theoretical Population Biology 49:
    369-386
    """
    char_x = []
    char_y = []
    state_alphabet = char_matrix.default_state_alphabet
    for t in char_matrix.taxon_set:
        if t in taxon_groups[0]:
            char_x.append(char_matrix[t])
        else:
            char_y.append(char_matrix[t])

    diffs_x, mean_diffs_x, sq_diff_x = _count_differences(char_x, state_alphabet, ignore_uncertain)
    diffs_y, mean_diffs_y, sq_diff_y = _count_differences(char_y, state_alphabet, ignore_uncertain)
    d_x = diffs_x / probability.binomial_coefficient(len(char_x), 2)
    d_y = diffs_y / probability.binomial_coefficient(len(char_y), 2)
    d_xy = _average_number_of_pairwise_differences_between_populations(char_x, char_y, state_alphabet, ignore_uncertain)
    s2_x = (float(sq_diff_x) / probability.binomial_coefficient(len(char_x), 2) ) - (d_x ** 2)
    s2_y = (float(sq_diff_y) / probability.binomial_coefficient(len(char_y), 2) ) - (d_y ** 2)
    s2_xy = _variance_of_pairwise_differences_between_populations(char_x, char_y, d_xy, state_alphabet, ignore_uncertain)

    n = len(char_matrix)
    n_x = float(len(char_x))
    n_y = float(len(char_y))
    a = float(n * (n-1))
    ax = float(n_x * (n_x - 1))
    ay = float(n_y * (n_y - 1))
    k = average_number_of_pairwise_differences(char_matrix, ignore_uncertain)
    psi = (float(1)/(a)) * ( ax * (math.sqrt(s2_x)/d_x) + ay * (math.sqrt(s2_y)/d_y) + (2 * n_x * n_y * math.sqrt(s2_xy)/k))
    return psi

def summarize(char_matrix, taxon_groups, ignore_uncertain=True):
    """
    Returns a summary of a set of sequences that can be partitioned into
    the list of lists of taxa given by `taxon_groups`.
    """
    char_x = []
    char_y = []
    state_alphabet = char_matrix.default_state_alphabet
    for t in char_matrix.taxon_set:
        if t in taxon_groups[0]:
            char_x.append(char_matrix[t])
        else:
            char_y.append(char_matrix[t])

    diffs_x, mean_diffs_x, sq_diff_x = _count_differences(char_x, state_alphabet, ignore_uncertain)
    diffs_y, mean_diffs_y, sq_diff_y = _count_differences(char_y, state_alphabet, ignore_uncertain)
    d_x = diffs_x / probability.binomial_coefficient(len(char_x), 2)
    d_y = diffs_y / probability.binomial_coefficient(len(char_y), 2)
    d_xy = _average_number_of_pairwise_differences_between_populations(char_x, char_y, state_alphabet, ignore_uncertain)
    s2_x = (float(sq_diff_x) / probability.binomial_coefficient(len(char_x), 2) ) - (d_x ** 2)
    s2_y = (float(sq_diff_y) / probability.binomial_coefficient(len(char_y), 2) ) - (d_y ** 2)
    s2_xy = _variance_of_pairwise_differences_between_populations(char_x, char_y, d_xy, state_alphabet, ignore_uncertain)

    n = len(char_matrix)
    n_x = float(len(char_x))
    n_y = float(len(char_y))
    a = float(n * (n-1))
    ax = float(n_x * (n_x - 1))
    ay = float(n_y * (n_y - 1))
    k = average_number_of_pairwise_differences(char_matrix, ignore_uncertain)

    summary = {}

    # Hickerson 2006: pi #
    summary["k"] = k

    # Hickerson 2006: pi_b #
    summary["k_b"] = d_xy

    # Hickerson 2006: pi_w #
    summary["k_w"] = d_x + d_y

    # Hickerson 2006: pi_net #
    summary["k_net"] = d_xy - (d_x + d_y)

    # Hickerson 2006: S #
    summary["S"] = num_segregating_sites(char_matrix, ignore_uncertain)

    # Hickerson 2006: theta #
    a1 = sum([1.0/i for i in range(1, len(char_matrix))])
    summary["theta"] = float(summary["S"]) / a1

    # Wakeley 1996 #
    summary["psi"] = (float(1)/(a)) * ( ax * (math.sqrt(s2_x)/d_x) + ay * (math.sqrt(s2_y)/d_y) + (2 * n_x * n_y * math.sqrt(s2_xy)/k))

    # Tajima's D #
    n = len(char_matrix)
    a1 = sum([1.0/i for i in range(1,n)])
    a2 = sum([1.0/(i**2) for i in range(1,n)])
    b1 = float(n+1)/(3*(n-1))
    b2 = float(2 * ( (n**2) + n + 3 )) / (9*n*(n-1))
    c1 = b1 - 1.0/a1
    c2 = b2 - float(n+2)/(a1 * n) + float(a2)/(a1 ** 2)
    e1 = float(c1) / a1
    e2 = float(c2) / ( (a1**2) + a2 )
    S = summary["S"]
    D = float(k - S) / math.sqrt( (e1 * S ) + (e2 * S) * (S -1) )
    summary["D"] = D

    return summary

def derived_state_matrix(char_vectors, ancestral_seq=None):
    """
    Given a list of CharVector objects, and a reference ancestral sequence,
    this returns a list of strings corresponding to the list of CharVector
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
    m = []
    for cv in char_vectors:
        m.append([])
        for i, s in enumerate(cv):
            if cv[i].value is ancestral_seq[i].value:
                m[-1].append(0)
            else:
                m[-1].append(1)
    return m

def unfolded_site_frequency_spectrum(char_vectors, ancestral_seq=None, pad=True):
    """
    Returns the site frequency spectrum of list of CharVector objects given by char_vectors,
    with reference to the ancestral sequence given by ancestral_seq. If ancestral_seq
    is None, then the first sequence in char_vectors is taken to be the ancestral
    sequence.
    """
    if ancestral_seq is None:
        ancestral_seq = char_vectors[0]
    dsm = derived_state_matrix(char_vectors, ancestral_seq)
    sites = zip(*dsm) # transpose
    freqs = {}
    if pad:
        for i in range(len(char_vectors)+1):
            freqs[i] = 0
    for s in sites:
        p = sum(s)
        if p not in freqs:
            freqs[p] = 1
        else:
            freqs[p] += 1
    return freqs


