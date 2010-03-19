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
    ###     http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/python/
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
    p = float(probability.binomial_coefficient(n_genes, 2)) / N
    tmrca = probability.geometric_rv(p)
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
    rate = probability.binomial_coefficient(n_genes, 2)
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
    nc2 = probability.binomial_coefficient(n_genes, 2)
    tmrca = (float(1)/nc2)
    if pop_size and pop_size >= 0:
        return tmrca * pop_size
    else:
        return tmrca

def coalesce(nodes,
             pop_size=None,
             period=None,
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

def node_waiting_time_pairs(tree):
    """Returns list of tuples of (node, coalescent interval [= time between
    last coalescent event and current node age])"""
    tree.add_ages_to_nodes(attr_name='age')
    ages = [(n, n.age) for n in tree.internal_nodes()]
    ages.sort(lambda x, y: int(x[1] - y[1]))
    intervals = []
    intervals.append(ages[0])
    for i, d in enumerate(ages[1:]):
        intervals.append( (d[0], d[1] - ages[i][1]) )
    return intervals

def extract_coalescent_frames(tree):
    """Returns dictionary, with key = number of alleles, and values = waiting time for
    coalescent for the given tree"""
    nwti = node_waiting_time_pairs(tree)
#     num_genes = len(tree.taxon_set)
    num_genes = len(tree.leaf_nodes())
    num_genes_wt = {}
    for n in nwti:
        num_genes_wt[num_genes] = n[1]
        num_genes = num_genes - len(n[0].child_nodes()) + 1
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

def log_probability_of_coalescent_tree(tree, haploid_pop_size):
    """
    Wraps up extraction of coalescent frames and reporting of probability.
    """
    return log_probability_of_coalescent_frames(extract_coalescent_frames(tree), haploid_pop_size)

def num_deep_coalescences_with_fitted_tree(gene_tree, species_tree):
    """
    Given two trees (with splits encoded), this returns the number of gene
    duplications implied by the gene tree reconciled on the species tree, based
    on the algorithm described here:

        Goodman, M. J. Czelnusiniak, G. W. Moore, A. E. Romero-Herrera, and
        G. Matsuda. 1979. Fitting the gene lineage into its species lineage,
        a parsimony strategy illustrated bu cladograms constructed from globin
        sequences. Syst. Zool. 19: 99-113.

        Maddison, W. P. 1997. Gene trees in species dataobject. Syst. Biol. 46:
        523-536.

    Note that for correct results,

        (a) trees must be rooted (i.e., is_rooted = True)
        (b) split masks must have been added as rooted (i.e., when
            encode_splits was called, is_rooted must have been set to True)

    """
    taxa_mask = species_tree.taxon_set.all_taxa_bitmask()

    species_node_gene_nodes = {}
    gene_node_species_nodes = {}

    for gnd in gene_tree.postorder_node_iter():
        gn_children = gnd.child_nodes()
        if len(gn_children) > 0:
            ssplit = 0
            for gn_child in gn_children:
                ssplit = ssplit | gene_node_species_nodes[gn_child].edge.split_bitmask
            sanc = treecalc.mrca(species_tree.seed_node, ssplit, taxa_mask)
            gene_node_species_nodes[gnd] = sanc
            if sanc not in species_node_gene_nodes:
                species_node_gene_nodes[sanc] = []
            species_node_gene_nodes[sanc].append(gnd)
        else:
            gene_node_species_nodes[gnd] = species_tree.find_node(lambda x : x.taxon == gnd.taxon)

    contained_gene_lineages = {}
    for snd in species_tree.postorder_node_iter():
        if snd in species_node_gene_nodes:
            for gnd in species_node_gene_nodes[snd]:
                for gnd_child in gnd.child_nodes():
                    sanc = gene_node_species_nodes[gnd_child]
                    p = sanc
                    while p is not None and p != snd:
                        if p.edge not in contained_gene_lineages:
                            contained_gene_lineages[p.edge] = 0
                        contained_gene_lineages[p.edge] += 1
                        p = p.parent_node

    dc = 0
    for v in contained_gene_lineages.values():
        dc += v - 1

    return dc

def num_deep_coalescences_with_grouping(tree, tax_sets):
    """
    Returns the number of deep coalescences on tree `tree` that would result
    if the taxa in `tax_sets` formed K mutually-exclusive monophyletic groups,
    where K = len(tax_sets)
    `tax_sets` == partition of taxa: list of lists, with the inner lists
    consisting of taxon objects forming a monophyletic group.
    """
    dc_tree = dataobject.Tree()
    dc_tree.taxon_set = dataobject.TaxonSet()

    for t in range(len(tax_sets)):
        dc_tree.taxon_set.append(dataobject.Taxon(label=str(t)))

    def _get_dc_taxon(nd):
        for idx, tax_set in enumerate(tax_sets):
            if nd.taxon in tax_set:
                return dc_tree.taxon_set[idx]
        assert "taxon not found in partition: '%s'" % nd.taxon.label

    src_dc_map = {}
    for snd in tree.postorder_node_iter():
        nnd = dataobject.Node()
        src_dc_map[snd] = nnd
        children = snd.child_nodes()
        if len(children) == 0:
            nnd.taxon = _get_dc_taxon(snd)
        else:
            taxa_set = []
            for cnd in children:
                dc_node = src_dc_map[cnd]
                if len(dc_node.child_nodes()) > 1:
                    nnd.add_child(dc_node)
                else:
                    ctax = dc_node.taxon
                    if ctax is not None and ctax not in taxa_set:
                        taxa_set.append(ctax)
                    del src_dc_map[cnd]
            if len(taxa_set) > 1:
                for t in taxa_set:
                    cnd = dataobject.Node()
                    cnd.taxon = t
                    nnd.add_child(cnd)
            else:
                if len(nnd.child_nodes()) == 0:
                    nnd.taxon = taxa_set[0]
                elif len(taxa_set) == 1:
                    cnd = dataobject.Node()
                    cnd.taxon = taxa_set[0]
                    nnd.add_child(cnd)
    dc_tree.seed_node = nnd
    return len(dc_tree.leaf_nodes()) - len(tax_sets)

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

