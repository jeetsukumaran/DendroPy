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
Population genetic simlations.
"""

import random
import copy

from dendropy.utility import GLOBAL_RNG
from dendropy.interop import seqgen
from dendropy.model import discrete
from dendropy.model import coalescent
import dendropy

class FragmentedPopulations(object):

    def __init__(self,
                 div_time_gens,
                 num_desc_pops = 2,
                 mutrate_per_site_per_generation=10e-8,
                 desc_pop_size=10000,
                 use_seq_gen=False,
                 rng=GLOBAL_RNG):
        """
        __init__ arguments:

            - ``div_time_gens`` : generations since divergence,
            - ``num_desc_pops`` : number of descendent populations,
            - ``mutrate_per_site_per_generation`` : sequence mutation rate, per-site per-generation
            - ``desc_diploid_pop_size`` : descendent lineage population size (=N; ancestral pop size = num_desc_pops * N)
            - ``rng`` : random number generator
        """
        self.div_time_gens = div_time_gens
        self.num_desc_pops = num_desc_pops
        self.mutrate_per_site_per_generation = mutrate_per_site_per_generation
        self.desc_pop_size = desc_pop_size
        self.rng = rng
        self.kappa = 1.0
        self.base_freqs=[0.25, 0.25, 0.25, 0.25]
        self.seqgen_path = 'seq-gen'
        self.use_seq_gen = use_seq_gen
        self.gene_tree = None
        self.pop_tree = None
        self.mutation_tree = None

    def _get_theta(self):
        return 4 * self.mutrate_per_gene_per_generation * self.desc_pop_size

    def generate_sequences(self,
                           species_name,
                           samples_per_pop=10,
                           seq_len=2000,
                           use_seq_gen=True):

        self.generate_pop_tree(species_name=species_name, samples_per_pop=samples_per_pop)
        self.generate_gene_tree(species_name=species_name, samples_per_pop=samples_per_pop)
        d = dendropy.DataSet(self.mutation_tree.taxon_namespace)
        if self.use_seq_gen is True:
            sg = seqgen.SeqGen()
            sg.seqgen_path = self.seqgen_path
            sg.num_replicates = 1
            sg.quiet = True
            sg.rng = self.rng
            sg.seq_len = seq_len
            sg.char_model = 'HKY'
            sg.ti_tv = float(self.kappa) / 2
            sg.state_freqs = self.base_freqs
            sg.trees = [self.mutation_tree]
            d = sg.generate_dataset(dataset=d)
        else:
            char_matrix = discrete.hky85_chars(
                    seq_len=seq_len,
                    tree_model=self.mutation_tree,
                    mutation_rate=1.0,
                    kappa=1.0,
                    base_freqs=[0.25, 0.25, 0.25, 0.25],
                    root_states=None,
                    rng=self.rng)
            d.add_char_matrix(char_matrix)
        return d

    def generate_pop_tree(self, species_name, samples_per_pop=10):
        tree_data = { 'sp': species_name, 'divt': self.div_time_gens }
        desc_lineages = []
        for i in xrange(self.num_desc_pops):
            tree_data['id'] = i+1
            desc_lineages.append("%(sp)s%(id)d:%(divt)d" % tree_data)
        tree_string = "(" + (",".join(desc_lineages)) + ("):%d" % 0) #% (self.num_desc_pops * self.desc_pop_size * 10))
        self.pop_tree = dendropy.Tree.get_from_string(tree_string, schema="newick")
        return self.pop_tree

    def generate_gene_tree(self, species_name, samples_per_pop=10):
        """
        Given:
            ``species_name`` : string identifying species/taxon
            ``samples_per_pop`` : number of samples (genes) per population
        Returns:
            DendroPy tree, with branch lengths in generations
        """
        if self.pop_tree is None:
            self.generate_pop_tree(species_name, samples_per_pop=10)
        for idx, leaf in enumerate(self.pop_tree.leaf_iter()):
            if idx == 1:
                # ancestral population = num_desc_pops * desc population
                leaf.parent_node.edge.pop_size = self.num_desc_pops * self.desc_pop_size
            leaf.edge.pop_size = self.desc_pop_size
            leaf.num_genes = samples_per_pop
        self.gene_tree, self.pop_tree = coalescent.constrained_kingman_tree(self.pop_tree,
                                                          gene_node_label_fn=lambda x,y: "%sX%d" % (x,y),
                                                          rng=self.rng)

        self.mutation_tree = copy.deepcopy(self.gene_tree)
        for edge in self.mutation_tree.preorder_edge_iter():
            edge.length = edge.length * self.mutrate_per_site_per_generation
        return self.gene_tree

def pop_gen_tree(tree=None,
                 taxon_namespace=None,
                 ages=None,
                 num_genes=None,
                 pop_sizes=None,
                 num_genes_attr = 'num_genes',
                 pop_size_attr = 'pop_size',
                 rng=None):
    """
    This will simulate and return a tree with edges decorated with
    population sizes and leaf nodes decorated by the number of genes
    (samples or lineages) in each leaf.

    If ``tree`` is given, then this is used as the tree to be decorated.
    Otherwise, a Yule tree is generated based on the given taxon_namespace.
    Either ``tree`` or ``taxon_namespace`` must be given.

    The timing of the divergences can be controlled by specifying a vector of
    ages, ``ages``. This should be sequences of values specifying the ages of the
    first, second, third etc. divergence events, in terms of time from the
    present, specified either in generations (if the ``pop_sizes`` vector is
    given) or population units (if the pop_size vector is not given).
    If an ages vector is given and there are less than num_pops-1 of these,
    then an exception is raised.

    The number of gene lineages per population can be specified through
    the 'num_genes', which can either be an scalar integer or a list.
    If it is an integer, all the population get the same number of
    genes. If it is a list, it must be at least as long as num_pops.

    The population sizes of each edge can be specified using the ``pop_sizes``
    vector, which should be a sequence of values specifying the population
    sizes of the edges in postorder. If the pop_size vector is given, then it
    must be at least as long as there are branches on a tree, i.e. 2 * num_pops
    + 1, otherwise it is an error.  The population size should be the effective
    *haploid* population size; i.e., number of gene copies in the population: 2
    * N in a diploid population of N individuals, or N in a haploid population
    * of N individuals.

    If ``pop_size`` is 1 or 0 or None, then edge lengths of the tree are in
    haploid population units; i.e. where 1 unit of time equals 2N generations
    for a diploid population of size N, or N generations for a haploid
    population of size N. Otherwise edge lengths of the tree are in
    generations.

    This function first generates a tree using a pure-birth model with a
    uniform birth rate of 1.0. If an ages vector is given, it then sweeps
    through the internal nodes, assigning branch lengths such that the
    divergence events correspond to the ages in the vector. If a population
    sizes vector is given, it then visits all the edges in postorder, assigning
    population sizes to the attribute with the name specified in
    'pop_size_attr' (which is persisted as an annotation). During this, if an
    ages vector was *not* given, then the edge lengths are multiplied by the
    population size of the edge so the branch length units will be in
    generations. If an ages vector was given, then it is assumed that the ages
    are already in the proper scale/units.
    """

    # get our random number generator
    if rng is None:
        rng = GLOBAL_RNG # use the global rng by default

    # get a yule tree
    if not tree:
        if taxon_namespace:
            from dendropy.simulate import treesim
            tree = treesim.uniform_pure_birth_tree(taxon_namespace=taxon_namespace,
                                      rng=rng)
        else:
            raise Exception("Either tree or taxon namespace must be given")

    num_pops = len(tree.leaf_nodes())

    # basic idiot-checking
    if ages is not None and len(ages) < (num_pops - 1):
        msg = "Too few ages specified."
        raise Exception(msg)
    if num_genes is not None:
        if isinstance(num_genes, list):
            if len(num_genes) < num_pops:
                msg = "Too few number of gene samples specified"
                raise Exception(msg)
            else:
                samples = num_genes
        else:
            samples = [num_genes for tax in range(num_pops)]
    else:
        samples = None
    if pop_sizes is not None and len(pop_sizes) < (2 * num_pops + 1):
        msg = "Too few population sizes specified."
        raise Exception(msg)

    # set the ages
    if ages is not None:

        # get the internal nodes on the tree in reverse branching
        # order, so that newest nodes are returned first
        nodes = tree.nodes(filter_fn = lambda x : not x.is_leaf())
        nodes.sort(key=lambda x: x.distance_from_root())
        # assign the ages
        for index, node in enumerate(nodes):
            for child in node.child_nodes():
                child.edge.length = ages[index] - child.distance_from_tip()

    # set the gene samples
    if samples is not None:
        for index, leaf in enumerate(tree.leaf_node_iter()):
            setattr(leaf, num_genes_attr, samples[index])
            leaf.annotations.add_bound_attribute(num_genes_attr)

    # set the population sizes
    if pop_sizes is not None:
        index = 0
        for edge in tree.postorder_edge_iter():
            setattr(edge, pop_size_attr, pop_sizes[index])
            edge.annotations.add_bound_attribute(pop_size_attr)
            if ages is None:
                edge.length = edge.length * getattr(edge, pop_size_attr)
            index = index + 1

    return tree

