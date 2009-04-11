#! /usr/bin/env python

############################################################################
##  coalecsent.py
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
Methods for working with Kingman's n-coalescent framework.
"""
import math
from dendropy import GLOBAL_RNG
from dendropy import distributions
from dendropy import treecalc

def discrete_time_to_coalescence(n_genes, 
                                 pop_size=None, 
                                 haploid=True,
                                 rng=None):
    """
    A random draw from the "Kingman distribution" (discrete time version):
    Time to go from n genes to n-1 genes; i.e. waiting time until two
    lineages coalesce. **`pop_size` = HAPLOID population size in the default
    formulation!**
    """
    if pop_size is None or pop_size <= n_genes:
        raise Exception("Population size must be >> num genes")
    if haploid:
        N = pop_size
    else:
        N = pop_size * 2
    if rng is None:
        rng = GLOBAL_RNG    
    p = float(distributions.binomial_coefficient(n_genes, 2)) / N
    tmrca = distributions.geometric_rv(p)
    if pop_size is not None and pop_size >= 0:
        return tmrca * pop_size
    else:
        return tmrca
            
def time_to_coalescence(n_genes, pop_size=None, rng=None):
    """
    A random draw from the "Kingman distribution" (continuous time version):
    Time to go from n genes to n-1 genes; i.e. waiting time until two
    lineages coalesce.  This is a random number with an exponential
    distribution with a rate of (n choose 2). Time is in coalescent
    units unless population size is > 1.
    """
    if rng is None:
        rng = GLOBAL_RNG    
    rate = distributions.binomial_coefficient(n_genes, 2)
    tmrca = rng.expovariate(rate) 
    if pop_size is not None and pop_size >= 0:
        return tmrca * pop_size
    else:
        return tmrca

def expected_tmrca(n_genes, pop_size=None, rng=None):
    """
    Expected (mean) value for the TMRCA in coalescent time units unless
    population size > 1
    """
    if rng is None:
        rng = GLOBAL_RNG        
    nc2 = distributions.binomial_coefficient(n_genes, 2)
    tmrca = (float(1)/nc2)
    if pop_size and pop_size >= 0:
        return tmrca * pop_size
    else:
        return tmrca

def coalesce(nodes,
             pop_size=None,
             period=None,
             node_factory=None,
             rng=None):
    """
    `nodes` is a list of DendroPy Nodes representing a sample of
    neutral genes (some, all, or none of these nodes may have
    descendent nodes).

    `pop_size` is the effective population size of a Wright-Fisher
    population in which thes genes are evolving, and must be given for
    correct behaviour if time is generations.
    
    `period` is the time in generations (if pop_size is not None) or
    in populaton/coalescent units (if pop_size is given) that the
    genes have to coalesce. If not given, then the the method will
    continue to coalesce the genes until only one gene is left.

    This function will a draw a coalescence time, `t`, from
    EXP(1/num_genes). If `period` is given and if this time is less
    than `period`, or if `period` is not given, then two nodes are
    selected at random from `nodes`, and coalesced: a new node is
    created, and the two nodes are added as child_nodes to this node with
    an edge length such the the total length from tip to the ancestral
    node is equal to the depth of the deepest child + `t`. The two
    nodes are removed from the list of nodes, and the new node is
    added to it. `t` is then deducted from `period`, and the process
    repeats.

    The function ends and returns the list of nodes once `period` is
    exhausted or if any draw of `t` exceeds `period`, if `period` is
    given or when there is only one node left.

    As each coalescent event occurs, *all* nodes have their edges
    extended to the point of the coalescent event. In the case of
    constrained coalescence, all uncoalesced nodes have their edges
    extended to the end of the period (coalesced nodes have the edges
    fixed by the coalescent event in their ancestor).  Thus multiple
    calls to this method with the same set of nodes will gradually
    'grow' the edges, until all the the nodes coalesce. The edge
    lengths of the nodes passed to this method thus should not be
    modified or reset until the process is complete.
    """

    # idiot-check, because I can be an idiot
    if not nodes:
        return []

    # set the random number generator
    if rng is None:
        rng = GLOBAL_RNG    
    
    # define the function needed to create new coalescence nodes
    if node_factory is None:
        new_node = nodes[0].__class__
    else:
        new_node = node_factory

    # make a shallow copy of the node list
    nodes = list(nodes)

    # start tracking the time remaining
    time_remaining = period

    # If there is no time constraint, we want to continue coalescing
    # until there is only one gene left in the pool. If there is a
    # time constraint, we continue as long as there is time remaining,
    # but we do not control for that here: it is automatically taken
    # care of when the time drawn for the next coalescent event
    # exceeds the time remaining, and triggers a break from the loop
    while len(nodes) > 1:

        # draw a time to coalesce: this will be an exponential random
        # variable with parameter (rate) of BINOMIAL[n_genes 2]
        # multiplied pop_size
        tmrca = time_to_coalescence(len(nodes), pop_size=pop_size, rng=rng)

        # if no time_remaining is given (i.e, we want to coalesce till
        # there is only one gene left) or, if we are working under the
        # constrained coalescent, if the time to the next coalescence
        # event is not longer than the time_remaining
        if time_remaining is None or tmrca <= time_remaining:

            # stretch out the edges of all the nodes to this time
            for node in nodes:
                if node.edge.length is None:
                    node.edge.length = 0.0
                node.edge.length = node.edge.length + tmrca

            # pick two nodes to coalesce at random
            to_coalesce = rng.sample(nodes, 2)

            # create the new ancestor of these nodes
            new_ancestor = new_node()

            # add the nodes as child nodes of the new node, their
            # common ancestor, and set the ancestor's edge length
            new_ancestor.add_child(to_coalesce[0])
            new_ancestor.add_child(to_coalesce[1])
            new_ancestor.edge.length = 0.0

            # remove the nodes that have coalesced from the pool of
            # nodes
            nodes.remove(to_coalesce[0])
            nodes.remove(to_coalesce[1])

            # add the ancestor to the pool of nodes
            nodes.append(new_ancestor)            

            # adjust the time_remaining left to coalesce
            if time_remaining is not None:
                time_remaining = time_remaining - tmrca
        else:
            # the next coalescent event takes place after the period constraint
            break

    # adjust the edge lengths of all the nodes, so they are at the
    # correct height, with the edges 'lining up' at the end of
    # coalescent period
    if time_remaining is not None and time_remaining > 0:
        for node in nodes:
            if node.edge.length is None:
                node.edge.length = 0.0
            node.edge.length = node.edge.length + time_remaining

    # return the list of nodes that have not coalesced
    return nodes

def coalescence_intervals(tree):
    """Returns list of coalescence intervals on `tree`."""
    treecalc.add_depth_to_nodes(tree)
    depths = [n.depth for n in tree.internal_nodes()]
    depths.sort()
    intervals = []
    intervals.append(depths[0])
    for i, d in enumerate(depths[1:]):
        intervals.append(d - depths[i])
    return intervals

def node_waiting_time_pairs(tree):
    """Returns list of tuples of (node, coalescent interval [= time between
    last coalescent event and current node age])"""
    treecalc.add_depth_to_nodes(tree)
    depths = [(n, n.depth) for n in tree.internal_nodes()]    
    depths.sort(lambda x, y: int(x[1] - y[1]))            
    intervals = []
    intervals.append(depths[0])
    for i, d in enumerate(depths[1:]):
        intervals.append( (d[0], d[1] - depths[i][1]) )
    return intervals

def num_genes_waiting_times_pairs(tree):
    """Returns list of tuples of (number of genes in the sample, waiting time
    till coalescent event / size of coalescent interval) for the given tree"""
    nwti = node_waiting_time_pairs(tree)
#     num_genes = len(tree.taxa_block)
    num_genes = len(tree.leaf_nodes())
    num_genes_wt = []
    for n in nwti:
        num_genes_wt.append((num_genes, n[1]))
        num_genes = num_genes - len(n[0].child_nodes()) + 1 
    return num_genes_wt      

def probability_of_coalescent_tree(tree, haploid_pop_size):
    """
    Under the classical neutral coalescent \citep{Kingman1982,
    Kingman1982b}, the waiting times between coalescent events in a
    sample of $k$ alleles segregating in a  population of (haploid) size
    $N_e$ is distributed exponentially with a rate parameter of
    $\frac{{k \choose 2}}{N_e}$:
    
    \begin{align}
         \Pr(T) =  \frac{{k \choose 2}}{N_e} \e{-  \frac{{k \choose 2}}{N_e} T},
    \end{align}
    
    where $T$ is the length of  (chronological) time in which there are
    $k$ alleles in the sample (i.e., for $k$ alleles to coalesce into
    $k-1$ alleles).
    """
    kts = num_genes_waiting_times_pairs(tree)  
#     p = 1.0
#     for kt in kts:
#         print "****"
#         print "p = %s" % p
#         print "k = %s" % kt[0]
#         print "t = %s" % kt[1]
#         print "k2N = %s" % k2N
#         print "e^(-k2N * t) = %s" % math.exp(-k2N * kt[1])
#         k2N = float(distributions.binomial_coefficient(kt[0], 2)) / haploid_pop_size
#         p *=  k2N * math.exp(-k2N * kt[1])
#         print "p' = %s" %  p
    lp = 0.0
    for kt in kts:
        k2N = float(distributions.binomial_coefficient(kt[0], 2)) / haploid_pop_size
        lp =  lp + math.log(k2N) - (k2N * kt[1])
    p = math.exp(lp)
    return p