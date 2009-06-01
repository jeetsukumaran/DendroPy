#! /usr/bin/env python

############################################################################
##  popgenstats.py
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
Population genetic statistics.
"""

import math
from dendropy import distributions

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
    return sum_diff / distributions.binomial_coefficient(len(char_vectors), 2)
    
def average_number_of_pairwise_differences(char_block, ignore_uncertain=True):
    """
    Returns $k$, calculated for a character block.
    """
    return _average_number_of_pairwise_differences(char_block.vectors(), char_block.default_state_alphabet, ignore_uncertain)

def nucleotide_diversity(char_block, ignore_uncertain=True):
    """
    Returns $\pi$, calculated for a character block.
    """
    return _nucleotide_diversity(char_block.vectors(), char_block.default_state_alphabet, ignore_uncertain)

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

def wakeleys_Psi(char_block, taxon_groups, ignore_uncertain=True):
    """         
    Returns Wakeley's Sigma, as described in:
    
    Wakeley, J. 1996. Distinguishing migration from isolation using the 
    variance of pairwise differences. Theoretical Population Biology 49: 
    369-386   
    """
    char_x = []
    char_y = []
    state_alphabet = char_block.default_state_alphabet
    for t in char_block.taxa_block:
        if t in taxon_groups[0]:
            char_x.append(char_block[t])
        else:
            char_y.append(char_block[t])
            
    diffs_x, mean_diffs_x, sq_diff_x = _count_differences(char_x, state_alphabet, ignore_uncertain)            
    diffs_y, mean_diffs_y, sq_diff_y = _count_differences(char_y, state_alphabet, ignore_uncertain)
    d_x = diffs_x / distributions.binomial_coefficient(len(char_x), 2)
    d_y = diffs_y / distributions.binomial_coefficient(len(char_y), 2)
    d_xy = _average_number_of_pairwise_differences_between_populations(char_x, char_y, state_alphabet, ignore_uncertain)
    s2_x = (float(sq_diff_x) / distributions.binomial_coefficient(len(char_x), 2) ) - (d_x ** 2)
    s2_y = (float(sq_diff_y) / distributions.binomial_coefficient(len(char_y), 2) ) - (d_y ** 2)
    s2_xy = _variance_of_pairwise_differences_between_populations(char_x, char_y, d_xy, state_alphabet, ignore_uncertain)
    
    n = len(char_block)
    n_x = float(len(char_x))
    n_y = float(len(char_y))
    a = float(n * (n-1))
    ax = float(n_x * (n_x - 1))
    ay = float(n_y * (n_y - 1))
    k = average_number_of_pairwise_differences(char_block, ignore_uncertain)
    sigma = (float(1)/(a)) * ( ax * (math.sqrt(s2_x)/d_x) + ay * (math.sqrt(s2_y)/d_y) + (2 * n_x * n_y * math.sqrt(s2_xy)/k))
    return sigma