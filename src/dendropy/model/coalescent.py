#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
Functions, classes, and methods for working with Kingman's n-coalescent
framework.
"""

import math
import dendropy
from dendropy.utility import GLOBAL_RNG
from dendropy.utility import constants
from dendropy.calculate import probability
from dendropy.calculate import combinatorics
from dendropy.utility.error import ProcessFailedException

###############################################################################
## Calculations and statistics


def discrete_time_to_coalescence(n_genes, pop_size=None, n_to_coalesce=2, rng=None):
    """
    A random draw from the "Kingman distribution" (discrete time version): Time
    to go from ``n_genes`` genes to ``n_genes``-1 genes in a discrete-time
    Wright-Fisher population of ``pop_size`` genes; i.e. waiting time until
    ``n-genes`` lineages coalesce in a population of ``pop_size`` genes.

    Parameters
    ----------

    n_genes : integer
        The number of genes in the sample. Must be greater than or equal to n_to_coalesce.
    pop_size : integer
        The effective *haploid* population size; i.e., number of genes in the
        population: 2 * N in a diploid population of N individuals, or N in a
        haploid population of N individuals.
    n_to_coalesce : integer
        The waiting time that will be returned will be the waiting time for
        this number of genes in the sample to coalesce.
    rng : ``Random`` object
        The random number generator instance.

    Returns
    -------
    k : integer
        A randomly-generated waiting time (in discrete generations) for
        ``n_to_coalesce`` genes to coalesce out of a sample of ``n_genes`` in a
        population of ``pop_size`` genes.

    """
    if not pop_size:
        time_units = 1.0
    else:
        time_units = pop_size
    if rng is None:
        rng = GLOBAL_RNG
    p = pop_size / combinatorics.choose(n_genes, n_to_coalesce)
    tmrca = probability.geometric_rv(p)
    return tmrca * time_units


def time_to_coalescence(n_genes, pop_size=None, n_to_coalesce=2, rng=None):
    r"""
    A random draw from the "Kingman distribution" (discrete time version): Time
    to go from ``n_genes`` genes to ``n_genes``-1 genes in a continuous-time
    Wright-Fisher population of ``pop_size`` genes; i.e. waiting time until
    ``n-genes`` lineages coalesce in a population of ``pop_size`` genes.

    Given the number of gene lineages in a sample, ``n_genes``, and a
    population size, ``pop_size``, this function returns a random number from
    an exponential distribution with rate :math:`\choose(``pop_size``, 2)`.
    ``pop_size`` is the effective *haploid* population size; i.e., number of gene
    in the population: 2 * N in a diploid population of N individuals,
    or N in a haploid population of N individuals. If ``pop_size`` is 1 or 0 or
    None, then time is in haploid population units; i.e. where 1 unit of time
    equals 2N generations for a diploid population of size N, or N generations
    for a haploid population of size N. Otherwise time is in generations.

    The coalescence time, or the waiting time for the coalescence, of two
    gene lineages evolving in a population with haploid size :math:`N` is an
    exponentially-distributed random variable with rate of :math:`N` an
    expectation of :math:`\frac{1}{N}`).
    The waiting time for coalescence of *any* two gene lineages in a sample of
    :math:`n` gene lineages evolving in a population with haploid size :math:`N` is an
    exponentially-distributed random variable with rate of :math:`\choose{N, 2} and
    an expectation of :math:`\frac{1}{\choose{N, 2}}`.

    Parameters
    ----------
    n_genes : integer
        The number of genes in the sample.
    pop_size : integer
        The effective *haploid* population size; i.e., number of genes in the
        population: 2 * N in a diploid population of N individuals, or N in a
        haploid population of N individuals.
    n_to_coalesce : integer
        The waiting time that will be returned will be the waiting time for
        this number of genes in the sample to coalesce.
    rng : ``Random`` object
        The random number generator instance to use.

    Returns
    -------
    k : float
        A randomly-generated waiting time (in continuous time) for
        ``n_to_coalesce`` genes to coalesce out of a sample of ``n_genes`` in a
        population of ``pop_size`` genes.
    """
    if rng is None:
        rng = GLOBAL_RNG
    if not pop_size:
        time_units = 1.0
    else:
        time_units = pop_size
    rate = combinatorics.choose(n_genes, n_to_coalesce)
    tmrca = rng.expovariate(rate)
    return tmrca * time_units


def expected_tmrca(n_genes, pop_size=None, n_to_coalesce=2):
    """
    Expected (mean) value for the Time to the Most Recent Common Ancestor of
    ``n_to_coalesce`` genes in a sample of ``n_genes`` drawn from a population of
    ``pop_size`` genes.

    Parameters
    ----------
    n_genes : integer
        The number of genes in the sample.
    pop_size : integer
        The effective *haploid* population size; i.e., number of genes in the
        population: 2 * N in a diploid population of N individuals, or N in a
        haploid population of N individuals.
    n_to_coalesce : integer
        The waiting time that will be returned will be the waiting time for
        this number of genes in the sample to coalesce.
    rng : ``Random`` object
        The random number generator instance.

    Returns
    -------
    k : float
        The expected waiting time (in continuous time) for ``n_to_coalesce``
        genes to coalesce out of a sample of ``n_genes`` in a population of
        ``pop_size`` genes.

    """
    nc2 = combinatorics.choose(n_genes, n_to_coalesce)
    tmrca = float(1) / nc2
    if pop_size is not None:
        return tmrca * pop_size
    else:
        return tmrca


def coalesce_nodes(
    nodes, pop_size=None, period=None, rng=None, use_expected_tmrca=False
):
    """
    Returns a list of nodes that have not yet coalesced once ``period`` is
    exhausted.

    This function will a draw a coalescence time, ``t``, from an exponential
    distribution with a rate of ``choose(k, 2)``, where ``k`` is the number of
    nodes. If ``period`` is given and if this time is less than ``period``, or if
    ``period`` is not given, then two nodes are selected at random from ``nodes``,
    and coalesced: a new node is created, and the two nodes are added as
    child_nodes to this node with an edge length such the the total length from
    tip to the ancestral node is equal to the depth of the deepest child + ``t``.
    The two nodes are removed from the list of nodes, and the new node is added
    to it. ``t`` is then deducted from ``period``, and the process repeats.

    The function ends and returns the list of nodes once ``period`` is
    exhausted or if any draw of ``t`` exceeds ``period``, if ``period`` is
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

    Parameters
    ----------
    nodes : iterable[|Node|]
        An interable of |Node| objects representing a sample of neutral
        genes (some, all, or none of these nodes may have descendent nodes).
    pop_size : integer
        The effective *haploid* population size; i.e., number of genes in the
        population: 2 * N in a diploid population of N individuals, or N in a
        haploid population of N individuals.
    period : numeric
        The time that the genes have to coalesce. If ``pop_size`` is 1 or 0 or
        None, then time is in haploid population units; i.e. where 1 unit of
        time equals 2N generations for a diploid population of size N, or N
        generations for a haploid population of size N. Otherwise time is in
        generations.
    rng : ``Random`` object
        The random number generator instance to use. If not specified, the
        default RNG will be used.
    use_expected_tmrca : bool
        If |True|, then instead of random times, the *expected* times will be
        used.

    Returns
    -------
    nodes : iterable[|Node|]
        A list of nodes once ``period`` is exhausted or if any draw of ``t``
        exceeds ``period``, if ``period`` is given or when there is only one node
        left.
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


def node_waiting_time_pairs(
    tree, ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION
):
    """
    Returns a list of tuples of (nodes, coalescent interval time) on the tree.
    That is, each element in the list is tuple pair consisting of where: the
    first element of the pair is an internal node representing a coalescent
    event on the tree, and the second element of the pair is the time between
    this coalescence event and the earlier (more recent) one.

    Parameters
    ----------
    tree : |Tree|
        A tree instance.
    ultrametricity_precision : float
        When calculating the node ages, an error will be raised if the tree is
        not ultrametric. This error may be due to floating-point or numerical
        imprecision. You can set the precision of the ultrametricity validation
        by setting the ``ultrametricity_precision`` parameter. E.g., use
        ``ultrametricity_precision=0.01`` for a more relaxed precision, down to
        2 decimal places. Use ``ultrametricity_precision=False`` to disable
        checking of ultrametricity.

    Returns
    -------
    x : list of tuples (node, coalescent interval)
        Returns list of tuples of (node, coalescent interval [= time between
        last coalescent event and current node age])
    """
    tree.calc_node_ages(ultrametricity_precision=ultrametricity_precision)
    ages = [(n, n.age) for n in tree.internal_nodes()]
    ages.sort(key=lambda x: x[1])
    intervals = []
    intervals.append(ages[0])
    for i, d in enumerate(ages[1:]):
        nd = d[0]
        prev_nd = ages[i][0]
        intervals.append((nd, nd.age - prev_nd.age))
    return intervals


def extract_coalescent_frames(
    tree, ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION
):
    """
    Returns a list of tuples describing the coalescent frames on the tree. That
    is, each element in the list is tuple pair consisting of where: the first
    element of the pair is the number of separate lineages remaining on the
    tree at coalescence event, and the second element of the pair is the time
    between this coalescence event and the earlier (more recent) one.

    Parameters
    ----------
    tree : |Tree|
        A tree instance.
    ultrametricity_precision : float
        When calculating the node ages, an error will be raised if the tree is
        not ultrametric. This error may be due to floating-point or numerical
        imprecision. You can set the precision of the ultrametricity validation
        by setting the ``ultrametricity_precision`` parameter. E.g., use
        ``ultrametricity_precision=0.01`` for a more relaxed precision, down to
        2 decimal places. Use ``ultrametricity_precision=False`` to disable
        checking of ultrametricity.

    Returns
    -------
    x : dict
        Returns dictionary, with key = number of alleles, and values = waiting
        time for coalescent for the given tree
    """
    nwti = node_waiting_time_pairs(
        tree, ultrametricity_precision=ultrametricity_precision
    )
    #     num_genes = len(tree.taxon_namespace)
    num_genes = len(tree.leaf_nodes())
    num_genes_wt = {}
    for n in nwti:
        num_genes_wt[num_genes] = n[1]
        num_genes = num_genes - len(n[0].child_nodes()) + 1
    # num_alleles_list = sorted(num_genes_wt.keys(), reverse=True)
    return num_genes_wt


def log_probability_of_coalescent_frames(coalescent_frames, haploid_pop_size):
    r"""
    Under the classical neutral coalescent :math:`\citep{Kingman1982,
    Kingman1982b}`, the waiting times between coalescent events in a
    sample of :math:`k` alleles segregating in a  population of (haploid) size
    :math:`N_e` is distributed exponentially with a rate parameter of
    :math:`\frac{{k \choose 2}}{N_e}`:

        .. math::

            \Pr(T) =  \frac{{k \choose 2}}{N_e} e^{-  \frac{{k \choose 2}}{N_e} T},

    where :math:`T` is the length of  (chronological) time in which there are
    :math:`k` alleles in the sample (i.e., for :math:`k` alleles to coalesce into
    :math:`k-1` alleles).
    """
    lp = 0.0
    for k, t in coalescent_frames.items():
        k2N = (float(k * (k - 1)) / 2) / haploid_pop_size
        #         k2N = float(combinatorics.choose(k, 2)) / haploid_pop_size
        lp = lp + math.log(k2N) - (k2N * t)
    return lp


def log_probability_of_coalescent_tree(
    tree,
    haploid_pop_size,
    ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
):
    """
    Wraps up extraction of coalescent frames and reporting of probability.
    """
    return log_probability_of_coalescent_frames(
        extract_coalescent_frames(tree), haploid_pop_size
    )


###############################################################################
## Tree Simulations


def contained_coalescent_tree(
    containing_tree,
    gene_to_containing_taxon_map,
    edge_pop_size_attr="pop_size",
    default_pop_size=1,
    rng=None,
):
    """
    Returns a gene tree simulated under the coalescent contained within a
    population or species tree.

        ``containing_tree``
            The population or species tree. If ``edge_pop_size_map`` is not None,
            and population sizes given are non-trivial (i.e., >1), then edge
            lengths on this tree are in units of generations. Otherwise edge
            lengths are in population units; i.e. 2N generations for diploid
            populations of size N, or N generations for diploid populations of
            size N.

        ``gene_to_containing_taxon_map``
            A TaxonNamespaceMapping object mapping Taxon objects in the
            ``containing_tree`` TaxonNamespace to corresponding Taxon objects in the
            resulting gene tree.

        ``edge_pop_size_attr``
            Name of attribute of edges that specify population size. By default
            this is "pop_size". If this attribute does not exist,
            ``default_pop_size`` will be used.  The value for this attribute
            should be the haploid population size or the number of genes;
            i.e.  2N for a diploid population of N individuals, or N for a
            haploid population of N individuals. This value determines how
            branch length units are interpreted in the input tree,
            ``containing_tree``.  If a biologically-meaningful value, then branch
            lengths on the ``containing_tree`` are properly read as generations.
            If not (e.g. 1 or 0), then they are in population units, i.e. where
            1 unit of time equals 2N generations for a diploid population of
            size N, or N generations for a haploid population of size N.
            Otherwise time is in generations. If this argument is None, then
            population sizes default to ``default_pop_size``.

        ``default_pop_size``
            Population size to use if ``edge_pop_size_attr`` is None or
            if an edge does not have the attribute. Defaults to 1.

    The returned gene tree will have the following extra attributes:

        ``pop_node_genes``
            A dictionary with nodes of ``containing_tree`` as keys and a list of gene
            tree nodes that are uncoalesced as values.

    Note that this function does very much the same thing as
    ``constrained_kingman_tree()``, but provides a very different API.
    """

    if rng is None:
        rng = GLOBAL_RNG

    gene_tree_taxon_namespace = gene_to_containing_taxon_map.domain_taxon_namespace
    if gene_tree_taxon_namespace is None:
        gene_tree_taxon_namespace = dendropy.TaxonNamespace()
        for gene_taxa in gene_to_containing_taxon_map:
            for taxon in gene_taxa:
                gene_tree_taxon_namespace.add(taxon)
    gene_tree = dendropy.Tree(taxon_namespace=gene_tree_taxon_namespace)
    gene_tree.is_rooted = True

    pop_node_genes = {}
    pop_gene_taxa = gene_to_containing_taxon_map.reverse
    for nd in containing_tree.postorder_node_iter():
        if nd.taxon and nd.taxon in pop_gene_taxa:
            pop_node_genes[nd] = []
            gene_taxa = pop_gene_taxa[nd.taxon]
            for gene_taxon in gene_taxa:
                gene_node = dendropy.Node()
                gene_node.taxon = gene_taxon
                pop_node_genes[nd].append(gene_node)
            # gene_nodes = [dendropy.Node() for i in range(len(gene_taxa))]
            # for gidx, gene_node in enumerate(gene_nodes):
            #    gene_node.taxon = gene_taxa[gidx]
            #    pop_node_genes[nd].append(gene_node)

    for edge in containing_tree.postorder_edge_iter():

        if edge_pop_size_attr and hasattr(edge, edge_pop_size_attr):
            pop_size = getattr(edge, edge_pop_size_attr)
        else:
            pop_size = default_pop_size
        if edge.head_node.parent_node is None:
            if len(pop_node_genes[edge.head_node]) > 1:
                final = coalesce_nodes(
                    nodes=pop_node_genes[edge.head_node],
                    pop_size=pop_size,
                    period=None,
                    rng=rng,
                )
            else:
                final = pop_node_genes[edge.head_node]
            gene_tree.seed_node = final[0]
        else:
            uncoal = coalesce_nodes(
                nodes=pop_node_genes[edge.head_node],
                pop_size=pop_size,
                period=edge.length,
                rng=rng,
            )
            if edge.tail_node not in pop_node_genes:
                pop_node_genes[edge.tail_node] = []
            pop_node_genes[edge.tail_node].extend(uncoal)

    gene_tree.pop_node_genes = pop_node_genes
    return gene_tree


def pure_kingman_tree(taxon_namespace, pop_size=1, rng=None):
    """
    Generates a tree under the unconstrained Kingman's coalescent process.

    Parameters
    ----------
    taxon_namespace: |TaxonNamespace| instance
        A pre-populated |TaxonNamespace| where the contained |Taxon| instances
        represent the genes or individuals sampled from the population.
    pop_size : numeric
        The size of the population from the which the coalescent process is
        sampled.


    Returns
    -------
    t : |Tree|
        A tree sampled from the Kingman's neutral coalescent.

    """
    if rng is None:
        rng = GLOBAL_RNG  # use the global rng by default
    nodes = [dendropy.Node(taxon=t) for t in taxon_namespace]
    seed_node = coalesce_nodes(
        nodes=nodes, pop_size=pop_size, period=None, rng=rng, use_expected_tmrca=False
    )[0]
    tree = dendropy.Tree(taxon_namespace=taxon_namespace, seed_node=seed_node)
    return tree


def pure_kingman_tree_shape(num_leaves, pop_size=1, rng=None):
    """
    Like :func:`dendropy.model.pure_kingman_tree`, but does not assign taxa to tips.

    Parameters
    ----------
    num_leaves : int
        Number of individuals/genes sampled.
    pop_size : numeric
        The size of the population from the which the coalescent process is
        sampled.

    Returns
    -------
    t : |Tree|
        A tree sampled from the Kingman's neutral coalescent.

    """
    if rng is None:
        rng = GLOBAL_RNG  # use the global rng by default
    nodes = [dendropy.Node() for t in range(num_leaves)]
    seed_node = coalesce_nodes(
        nodes=nodes, pop_size=pop_size, period=None, rng=rng, use_expected_tmrca=False
    )[0]
    tree = dendropy.Tree(seed_node=seed_node)
    return tree


def mean_kingman_tree(taxon_namespace, pop_size=1, rng=None):
    """
    Returns a tree with coalescent intervals given by the expected times under
    Kingman's neutral coalescent.
    """
    if rng is None:
        rng = GLOBAL_RNG  # use the global rng by default
    nodes = [dendropy.Node(taxon=t) for t in taxon_namespace]
    seed_node = coalesce_nodes(
        nodes=nodes, pop_size=pop_size, period=None, rng=rng, use_expected_tmrca=True
    )[0]
    tree = dendropy.Tree(taxon_namespace=taxon_namespace, seed_node=seed_node)
    return tree


def constrained_kingman_tree(
    pop_tree,
    gene_tree_list=None,
    rng=None,
    gene_node_label_fn=None,
    gene_sampling_strategy="random_uniform",
    num_genes=None,
    num_genes_attr="num_genes",
    pop_size_attr="pop_size",
    decorate_original_tree=False,
):
    """
    Given a population tree, ``pop_tree`` this will return a *pair of
    trees*: a gene tree simulated on this population tree based on
    Kingman's n-coalescent, and population tree with the additional
    attribute 'gene_nodes' on each node, which is a list of
    uncoalesced nodes from the gene tree associated with the given
    node from the population tree.

    ``pop_tree``: a Tree object.

    ``gene_sampling_strategy``: string
        - "node_attribute": Will expect each leaf of ``pop_tree`` to
          have an attribute, ``num_genes``, that specifies the number
          of genes to be sampled from that population.
        - "fixed_per_population": Will assign ``num_genes`` to each population.
        - "random_uniform": Will assign genes to leaves with
          uniform probability until ``num_genes`` genes have been
          assigned.

    ``pop_size_attr``: string
        The attribute name of the edges of ``pop_tree`` that
        specify the population size. By default it is ``pop_size``. The should
        specify the effective *haploid* population size; i.e., number of gene
        in the population: 2 * N in a diploid population of N individuals,
        or N in a haploid population of N individuals.

    If ``pop_size`` is 1 or 0 or None, then the edge lengths of ``pop_tree`` is
    taken to be in haploid population units; i.e. where 1 unit equals 2N
    generations for a diploid population of size N, or N generations for a
    haploid population of size N. Otherwise the edge lengths of ``pop_tree`` is
    taken to be in generations.

    If ``gene_tree_list`` is given, then the gene tree is added to the
    tree block, and the tree block's taxa block will be used to manage
    the gene tree's ``taxa``.

    ``gene_node_label_fn`` is a function that takes two arguments (a string
    and an integer, respectively, where the string is the containing species
    taxon label and the integer is the gene index) and returns a label for
    the corresponding the gene node.

    if ``decorate_original_tree`` is True, then the list of uncoalesced nodes at
    each node of the population tree is added to the original (input) population
    tree instead of a copy.

    If ``num_genes`` is None, then it will be set to 1 under the
    "node_attribute" strategy (serving as a fallback default for nodes that do
    not spcify ``num_genes_attr``) or the leaf count of ``pop_tree`` under the
    ``random_uniform`` strategy.

    Note that this function does very much the same thing as
    ``contained_coalescent_tree()``,
    but provides a very different API.
    """

    # get our random number generator
    if rng is None:
        rng = GLOBAL_RNG  # use the global rng by default

    if gene_tree_list is not None:
        gtaxa = gene_tree_list.taxon_namespace
    else:
        gtaxa = dendropy.TaxonNamespace()

    if gene_node_label_fn is None:
        gene_node_label_fn = lambda x, y: "%s_%02d" % (x, y)

    # @MAM taking a stab at a reasonable default for num_genes,
    # it may make sense to do something else entirely here
    if num_genes is None:
        if gene_sampling_strategy == "random_uniform":
            num_genes = sum(1 for __ in pop_tree.leaf_node_iter())
        elif gene_sampling_strategy == "node_attribute":
            num_genes = 1
        else:
            num_genes = None

    # we create a set of gene nodes for each leaf node on the population
    # tree, and associate those gene nodes to the leaf by assignment
    # of 'taxon'.
    if gene_sampling_strategy in ("node_attribute", "fixed_per_population"):
        for leaf_count, leaf in enumerate(pop_tree.leaf_node_iter()):
            gene_nodes = []
            if gene_sampling_strategy == "node_attribute":
                node_ngenes = getattr(leaf, num_genes_attr)
            else:
                node_ngenes = num_genes
            for gene_count in range(node_ngenes):
                gene_node = dendropy.Node()
                gene_node.taxon = gtaxa.require_taxon(
                    label=gene_node_label_fn(leaf.taxon.label, gene_count + 1)
                )
                gene_nodes.append(gene_node)
            leaf.gene_nodes = gene_nodes
    elif gene_sampling_strategy == "random_uniform":
        gene_count = 0
        leaves = list(pop_tree.leaf_node_iter())
        while gene_count < num_genes:
            gene_count += 1
            leaf = rng.choice(leaves)
            gene_node = dendropy.Node()
            gene_node.taxon = gtaxa.require_taxon(
                label=gene_node_label_fn(leaf.taxon.label, gene_count)
            )
            try:
                leaf.gene_nodes.append(gene_node)
            except:
                leaf.gene_nodes = [ gene_node ]
    else:
        raise ValueError("Unrecognized strategy '{}'".format(
            gene_sampling_strategy
        ))

    # We iterate through the edges of the population tree in post-order,
    # i.e., visiting child edges before we visit parent edges. For
    # each edge visited, we take the genes found in the child nodes,
    # and run the coalescent simulation on them attacheded by the length
    # of the edge. Any genes that have not yet coalesced at the end of
    # this period are added to the genes of the tail (parent) node of
    # the edge.

    if decorate_original_tree:
        working_poptree = pop_tree
    else:
        # start with a new (deep) copy of the population tree so as to not
        # to change the original tree
        working_poptree = dendropy.Tree(pop_tree)

    # start with a new tree
    gene_tree = dendropy.Tree()
    gene_tree.taxon_namespace = gtaxa
    for edge in working_poptree.postorder_edge_iter():
        if not hasattr(edge.head_node, "gene_nodes"):
            continue
        # if mrca root, run unconstrained coalescent
        if hasattr(edge, pop_size_attr):
            pop_size = getattr(edge, pop_size_attr)
        else:
            # this means all our time will be in population units
            pop_size = 1.0
        if edge.head_node.parent_node is None:
            if len(edge.head_node.gene_nodes) > 1:
                final = coalesce_nodes(
                    nodes=edge.head_node.gene_nodes,
                    pop_size=pop_size,
                    period=None,
                    rng=rng,
                )
            else:
                final = edge.head_node.gene_nodes
            if not final:
                raise ProcessFailedException()
            gene_tree.seed_node = final[0]
        else:
            uncoal = coalesce_nodes(
                nodes=edge.head_node.gene_nodes,
                pop_size=pop_size,
                period=edge.length,
                rng=rng,
            )
            if not hasattr(edge.tail_node, "gene_nodes"):
                edge.tail_node.gene_nodes = []
            edge.tail_node.gene_nodes.extend(uncoal)

    gene_tree.is_rooted = True
    if gene_tree_list is not None:
        gene_tree_list.append(gene_tree)
        return gene_tree, working_poptree
    else:
        return gene_tree, working_poptree
