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

from dendropy import distributions

def _count_differences(char_vectors, state_alphabet, ignore_uncertain=True):
    """
    Returns pair of values: total number of pairwise differences observed between
    all sequences, and mean number of pairwise differences pair base.
    """
    sum_pi = 0.0
    mean_pi = 0.0
    total_counted = 0
    comps = 0
    for vidx, i in enumerate(char_vectors):
        for j in char_vectors[vidx+1:]:
            if len(i) != len(j):
                raise Exception("sequences of unequal length")
            diff = 0            
            counted = 0
            comps += 1
            for cidx, c in enumerate(i):
                c1 = c
                c2 = j[cidx]
                if (not ignore_uncertain) or (c1.value is not state_alphabet.gap and c2.value is not state_alphabet.gap and len(c1.value.fundamental_ids) == 1 and len(c2.value.fundamental_ids) == 1):                 
                    counted += 1
                    total_counted += 1
                    if c1.value is not c2.value:
                        diff += 1
            sum_pi += float(diff)    
            mean_pi += float(diff) / counted
    return sum_pi, mean_pi / comps
    
def _nucleotide_diversity(char_vectors, state_alphabet, ignore_uncertain=True):
    """
    Returns $k$ (Tajima 1983; Wakely 1996), calculated for a set of sequences:
    
    k = \frac{\right(\sum \sum \k_{ij}\left)}{n \choose 2}

    where $k_{ij}$ is the number of pairwise differences between the
    $i$th and $j$th sequence, and $n$ is the number of DNA sequences
    sampled.       
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
    sum_pi, mean_pi = _count_differences(char_vectors, state_alphabet, ignore_uncertain)    
    print sum_pi, mean_pi
    return sum_pi / distributions.binomial_coefficient(len(char_vectors), 2)
    
def average_number_of_pairwise_differences(char_block, ignore_uncertain=True):
    """
    Returns $k$, calculated for a character block.
    """
    return _average_number_of_pairwise_differences(char_block.vectors(), char_block.default_state_alphabet, ignore_uncertain)

def nucelotide_diversity(char_block, ignore_uncertain=True):
    """
    Returns $k$, calculated for a character block.
    """
    return _nucleotide_diversity(char_block.vectors(), char_block.default_state_alphabet, ignore_uncertain)


    
        