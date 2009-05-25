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

def _pi(char_vectors, ignore_uncertain=True):
    """
    Returns $\pi$, calculated for a set of sequences:
    
    \pi = \frac{\right(\sum \sum \pi_{ij}\left)}{n \choose 2}

    where $\pi_{ij}$ is the number of pairwise differences between the
    $i$th and $j$th sequence, and $n$ is the number of DNA sequences
    sampled
    """
    sum_pi = 0.0
    for vidx, i in enumerate(char_vectors):
        for j in char_vectors[vidx+1:]:
            if len(i) != len(j):
                raise Exception("sequences of unequal length")
            count = 0            
            for cidx, c in enumerate(i):
                c1 = c
                c2 = j[cidx]
                if not ignore_uncertain or (len(c1.value.fundamental_ids) == 1 and len(c2.value.fundamental_ids) == 1):
                    if c1.value != c2.value:
                        count += 1
            sum_pi += float(count) / len(i)                        
    return sum_pi / distributions.binomial_coefficient(len(char_vectors), 2)
    

def pi(char_block, ignore_uncertain=True):
    """
    Returns $\pi$, calculated for a character block.
    """
    return _pi(char_block.vectors(), ignore_uncertain)


if __name__ == "__main__":
    import StringIO
    from dendropy import datasets   
    d = datasets.Dataset()
    s = StringIO.StringIO("""
#NEXUS 

Begin data;
	Dimensions ntax=4 nchar=55;
	Format datatype=dna gap=-;
	Matrix
seq_1 ATATACGGGGTTA---TTAGA----AAAATGTGTGTGTGTTTTTTTTTTCATGTG
seq_2 ATATAC--GGATA---TTACA----AGAATCTATGTCTGCTTTCTTTTTCATGTG
seq_3 ATATACGGGGATA---TTATA----AGAATGTGTGTGTGTTTTTTTTTTCATGTG
seq_4 ATATACGGGGATA---GTAGT----AAAATGTGTGTGTGTTTTTTTTTTCATGTG
	;
End;
""")
    d.read(s, "NEXUS")
    print pi(d.char_blocks[0])
    
        