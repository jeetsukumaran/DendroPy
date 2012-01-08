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
Methods for working with Kingman's n-coalescent framework.
"""

import math

from dendropy.utility import GLOBAL_RNG
from dendropy.utility import probability
from dendropy import dataobject
from dendropy import treecalc

try:
    ### Required for kernel density estimation:
    ###     "Statistics" library
    ### By:
    ###     Michiel de Hoon
    ### From:
    ###     http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/python/Statistics
    ###
    ### Statistics for Python is an extension module, written in ANSI-C, for
    ### the Python scripting language. Currently, this extension module
    ### contains some routines to estimate the probability density function
    ### from a set of random variables.
    ### Statistics for Python was released under the Python License.
    ### Michiel de Hoon
    ### Center for Computational Biology and Bioinformatics,
    ### Columbia University.
    import statistics as de_hoon_lib
    de_hoon_statistics = True
except:
    de_hoon_statistics = False

def discrete_time_to_coalescence(n_genes,
                                 pop_size=None,
                                 rng=None):
    """
    A random draw from the "Kingman distribution" (discrete time version):
    Time to go from n genes to n-1 genes; i.e. waiting time until two
    lineages coalesce.
    `pop_size` is the effective *haploid* population size; i.e., number of
    genes in the population: 2 * N in a diploid population of N individuals, or
    N in a haploid population of N individuals.
    If `pop_size` is 1 or 0 or None, then time is in haploid population units;
    i.e. where 1 unit of time equals 2N generations for a diploid population of
    size N, or N generations for a haploid population of size N. Otherwise time
    is in generations.


    """
    if not pop_size:
        time_units = 1
    else:
        time_units = pop_size * 2
    if rng is None:
        rng = GLOBAL_RNG
    p = float(probability.binomial_coefficient(n_genes, 2)) / time_units
    tmrca = probability.geometric_rv(p)
    return tmrca * time_units

def time_to_coalescence(n_genes,
        pop_size=None,
        haploid=True,
        rng=None):
    """
    A random draw from the "Kingman distribution" (continuous time version):
    Time to go from n genes to n-1 genes; i.e. waiting time until two
    lineages coalesce.  This is a random number with an exponential
    distribution with a rate of (n choose 2).
    `pop_size` is the effective *haploid* population size; i.e., number of gene
    in the population: 2 * N in a diploid population of N individuals,
    or N in a haploid population of N individuals.
    If `pop_size` is 1 or 0 or None, then time is in haploid population units;
    i.e. where 1 unit of time equals 2N generations for a diploid population of
    size N, or N generations for a haploid population of size N. Otherwise time
    is in generations.

    """
    if rng is None:
        rng = GLOBAL_RNG
    if not pop_size:
        time_units = 1
    else:
        time_units = pop_size * 2
    rate = probability.binomial_coefficient(n_genes, 2)
    tmrca = rng.expovariate(rate)
    return tmrca * pop_size

def expected_tmrca(n_genes, pop_size=None):
    """
    Expected (mean) value for the Time to the Most Recent Common Ancestor.
    `n_genes` is the number of genes in the sample.
    `pop_size` is the effective *haploid* population size; i.e., number of gene
    in the population: 2 * N in a diploid population of N individuals,
    or N in a haploid population of N individuals.
    If `pop_size` is 1 or 0 or None, then time is in haploid population units;
    i.e. where 1 unit of time equals 2N generations for a diploid population of
    size N, or N generations for a haploid population of size N. Otherwise time
    is in generations.

    """
    nc2 = probability.binomial_coefficient(n_genes, 2)
    tmrca = (float(1)/nc2)
    return tmrca * pop_size

def coalesce(nodes,
             pop_size=None,
             period=None,
             rng=None,
             use_expected_tmrca=False):
    """
    Returns a list of nodes that have not yet coalesced once `period` is
    exhausted.

    `nodes` is a list of DendroPy Nodes representing a sample of
    neutral genes (some, all, or none of these nodes may have
    descendent nodes).

    `pop_size` is the effective *haploid* population size; i.e., number of gene
    in the population: 2 * N in a diploid population of N individuals,
    or N in a haploid population of N individuals.

    `period` is the time that the genes have to coalesce.  If `pop_size` is 1
    or 0 or None, then time is in haploid population units; i.e. where 1 unit
    of time equals 2N generations for a diploid population of size N, or N
    generations for a haploid population of size N. Otherwise time is in
    generations.

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
    new_node = nodes[0].__class__

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

        if use_expected_tmrca:
            tmrca = expected_tmrca(len(nodes), pop_size=pop_size)
        else:
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

def node_waiting_time_pairs(tree, check_ultrametricity_prec=0.0000001):
    """Returns list of tuples of (node, coalescent interval [= time between
    last coalescent event and current node age])"""
    tree.calc_node_ages(check_prec=check_ultrametricity_prec)
    ages = [(n, n.age) for n in tree.internal_nodes()]
    ages.sort(key=lambda x: x[1])
    intervals = []
    intervals.append(ages[0])
    for i, d in enumerate(ages[1:]):
        nd = d[0]
        prev_nd = ages[i][0]
        intervals.append( (nd, nd.age - prev_nd.age) )
    return intervals

def extract_coalescent_frames(tree, check_ultrametricity_prec=0.0000001):
    """Returns dictionary, with key = number of alleles, and values = waiting time for
    coalescent for the given tree"""
    nwti = node_waiting_time_pairs(tree, check_ultrametricity_prec=check_ultrametricity_prec)
#     num_genes = len(tree.taxon_set)
    num_genes = len(tree.leaf_nodes())
    num_genes_wt = {}
    for n in nwti:
        num_genes_wt[num_genes] = n[1]
        num_genes = num_genes - len(n[0].child_nodes()) + 1

    import sys
    num_alleles_list = sorted(num_genes_wt.keys(), reverse=True)
    return num_genes_wt

def log_probability_of_coalescent_frames(coalescent_frames, haploid_pop_size):
    """
    Under the classical neutral coalescent \citep{Kingman1982,
    Kingman1982b}, the waiting times between coalescent events in a
    sample of $k$ alleles segregating in a  population of (haploid) size
    $N_e$ is distributed exponentially with a rate parameter of
    $\frac{{k \choose 2}}{N_e}$:

         \Pr(T) =  \frac{{k \choose 2}}{N_e} \e{-  \frac{{k \choose 2}}{N_e} T},

    where $T$ is the length of  (chronological) time in which there are
    $k$ alleles in the sample (i.e., for $k$ alleles to coalesce into
    $k-1$ alleles).
    """
    lp = 0.0
    for k, t in coalescent_frames.items():
        k2N = (float(k * (k-1)) / 2) / haploid_pop_size
#         k2N = float(probability.binomial_coefficient(k, 2)) / haploid_pop_size
        lp =  lp + math.log(k2N) - (k2N * t)
    return lp

def log_probability_of_coalescent_tree(tree, haploid_pop_size, check_ultrametricity_prec=0.0000001):
    """
    Wraps up extraction of coalescent frames and reporting of probability.
    """
    return log_probability_of_coalescent_frames(extract_coalescent_frames(tree),
            haploid_pop_size)

if de_hoon_statistics:

    def kl_divergence_coalescent_trees(tree_list, haploid_pop_size):
        """
        Returns KL divergence for coalescent frames found in a collection of
        trees from the theoretical distribution given the specified haploid
        population size.
        """
        allele_waiting_time_dist = {}
        for t in tree_list:
            cf = extract_coalescent_frames(t)
            allele_waiting_time_dist = update_allele_waiting_time_dist(cf, allele_waiting_time_dist)
        return kl_divergence_coalescent_waiting_times(allele_waiting_time_dist, haploid_pop_size)

    def update_allele_waiting_time_dist(coalescent_frames, allele_waiting_time_dist=None):
        """
        `coalescent_frames` is a dictionary with number of alleles as keys and
        a scalar representing the waiting time to a coalescence event given a
        particular number of alleles on a particular tree (as returned by
        `extract_coalescent_frame`. `allele_branch_len_dist` is a dictionary
        with number of alleles as keys and a list of waiting times associated
        with that number of alleles as values. This is simply a convenience
        function that adds the waiting times found in `coalescent_frames`
        to the collection of values tracked in `allele_waiting_time_dist`.
        """
        if allele_waiting_time_dist is None:
            allele_waiting_time_dist = {}
        for k, t in coalescent_frames.items():
            if k not in allele_waiting_time_dist:
                allele_waiting_time_dist[k] = []
            allele_waiting_time_dist[k].append(t)
        return allele_waiting_time_dist

    def kl_divergence_coalescent_waiting_times(allele_waiting_time_dist, haploid_pop_size):
        """
        `allele_branch_len_dist` is a dictionary with number of alleles as keys
        and a list of waiting times associated with that number of alleles as
        values. `haploid_pop_size` is the population size in terms of total numbers
        of genes. This returns a the KL-divergence between the distribution of
        waiting times and the Kingman coalescent distribution.

        D_{\mathrm{KL}}(P\|Q) = \sum_i P(i) \log \frac{P(i)}{Q(i)}.

        """
        d_kl = 0.0
        for k, wts in allele_waiting_time_dist.items():
            p = float(probability.binomial_coefficient(k, 2)) / haploid_pop_size
            for t in wts:
                # Kernel types:
                #
                # 'E' or 'Epanechnikov'
                #     Epanechnikov kernel (default)
                #
                # 'U' or 'Uniform'
                #     Uniform kernel
                #
                # 'T' or 'Triangle'
                #     Triangle kernel
                #
                # 'G' or  'Gaussian'
                #     Gaussian kernel
                #
                # 'B' or 'Biweight'
                #     Quartic/biweight kernel
                #
                # '3' or 'Triweight'
                #     Triweight kernel
                #
                # 'C' or 'Cosine'
                #     Cosine kernel
                q = de_hoon_lib.pdf(wts, [k], kernel = 'Gaussian')
                if q == 0:
                    q = 1e-100
                d_kl += p * math.log(p/q)
        return d_kl

