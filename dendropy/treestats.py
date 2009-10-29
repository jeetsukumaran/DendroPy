############################################################################
##  treestats.py
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
Various statistics calculated on trees.
"""

import math

def pybus_harvey_gamma(tree, prec=0.00001):
    """Returns the gamma statistic of Pybus and Harvey (2000). This statistic 
    is used to test for constancy of birth and death rates over the course of
    a phylogeny.  Under the pure-birth process, the statistic should follow
    a standard Normal distibution: a Normal(mean=0, variance=1).
    
    If the lengths of different paths to the node differ by more than `prec`, 
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    Raises a Value Error if the tree is not ultrametric, is non-binary, or has
        only 2 leaves.
    
    As a side effect a `depth` attribute is added to the nodes of the tree.
    
    Pybus and Harvey. 2000. "Testing macro-evolutionary models using incomplete
    molecular phylogenies." Proc. Royal Society Series B: Biological Sciences.
    (267). 2267-2272
    """
    # the equation is given by:
    #   T = \sum_{j=2}^n (jg_j)
    #   C = T \sqrt{\frac{1}{12(n-2)}}
    #   C gamma = \frac{1}{n-2}\sum_{i=2}^{n-1} (\sum_{k=2}^i kg_k) - \frac{T}{2}
    # where n is the number of taxa, and g_2 ... g_n is the vector of waiting
    #   times between consecutive (in time, not along a branch) speciation times.
    node = None
    speciation_depths = []
    n = 0
    for node in tree.postorder_node_iter():
        ch = node.child_nodes()
        n_ch = len(ch)
        if n_ch == 0:
            node.depth = 0.0
            n += 1
        elif n_ch > 2:
            raise ValueError("Polytomy encountered")
        else:
            first_child = ch[0]
            node.depth = first_child.depth + first_child.edge.length
            last_child = ch[-1]
            for nnd in ch[1:]:
                ocnd = nnd.depth + nnd.edge.length
                if abs(node.depth - ocnd) > prec:
                    raise ValueError("Tree is not ultrametric")
            if n_ch == 2:
                speciation_depths.append(node.depth)
    if node is None:
        raise ValueError("Empty tree encountered")
    speciation_depths.sort(reverse=True)
    g = []
    older = speciation_depths[0]
    for age in speciation_depths[1:]:
        g.append(older - age)
        older = age
    g.append(older)
    if not g:
        raise ValueError("No internal nodes found (other than the root)")
    assert(len(g) == (n - 1))
    T = 0.0
    accum = 0.0
    for i in xrange(2, n):
        list_index = i - 2
        T += i * float(g[list_index])
        accum += T
    list_index = n - 2
    T += (n) * g[list_index]
    nmt = n - 2.0
    numerator = accum/nmt - T/2.0
    C = T*pow(1/(12*nmt), 0.5)
    return numerator/C

