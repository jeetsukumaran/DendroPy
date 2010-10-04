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
Population genetic simlations.
"""

import StringIO
import random
import copy

from dendropy.utility import GLOBAL_RNG
from dendropy import treesim
from dendropy import seqsim
import dendropy

try:
    from pyseqgen import seqgen
    SEQGEN = True
except:
    from dendropy import seqsim
    SEQGEN = False

class FragmentedPopulations(object):

    def __init__(self,
                 div_time_gens,
                 num_desc_pops = 2,
                 mutrate_per_site_per_generation=10e-8,
                 desc_pop_size=10000,
                 rng=GLOBAL_RNG):
        """
        __init__ arguments:

            - `div_time_gens` : generations since divergence,
            - `num_desc_pops` : number of descendent populations,
            - `mutrate_per_site_per_generation` : sequence mutation rate, per-site per-generation
            - `desc_diploid_pop_size` : descendent lineage population size (=N; ancestral pop size = num_desc_pops * N)
            - `rng` : random number generator
        """
        self.div_time_gens = div_time_gens
        self.num_desc_pops = num_desc_pops
        self.mutrate_per_site_per_generation = mutrate_per_site_per_generation
        self.desc_pop_size = desc_pop_size
        self.rng = rng
        self.kappa = 1.0
        self.base_freqs=[0.25, 0.25, 0.25, 0.25]
        self.seqgen_path = 'seq-gen'
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
        d = dendropy.DataSet(self.mutation_tree.taxon_set)

        if SEQGEN and use_seq_gen:

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

            return d
        else:
            return seqsim.generate_hky_dataset(seq_len=seq_len,
                                                tree_model=self.mutation_tree,
                                                mutation_rate=1.0,
                                                kappa=1.0,
                                                base_freqs=[0.25, 0.25, 0.25, 0.25],
                                                root_states=None,
                                                dataset=d,
                                                rng=self.rng)

    def generate_pop_tree(self, species_name, samples_per_pop=10):
        tree_data = { 'sp': species_name, 'divt': self.div_time_gens }
        desc_lineages = []
        for i in xrange(self.num_desc_pops):
            tree_data['id'] = i+1
            desc_lineages.append("%(sp)s%(id)d:%(divt)d" % tree_data)
        tree_string = "(" + (",".join(desc_lineages)) + ("):%d" % 0) #% (self.num_desc_pops * self.desc_pop_size * 10))
        self.pop_tree = dendropy.Tree(stream=StringIO.StringIO(tree_string), schema="newick")
        return self.pop_tree

    def generate_gene_tree(self, species_name, samples_per_pop=10):
        """
        Given:
            `species_name` : string identifying species/taxon
            `samples_per_pop` : number of samples (genes) per population
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
        self.gene_tree, self.pop_tree = treesim.constrained_kingman(self.pop_tree,
                                                          gene_node_label_func=lambda x,y: "%sX%d" % (x,y),
                                                          rng=self.rng)

        self.mutation_tree = copy.deepcopy(self.gene_tree)
        for edge in self.mutation_tree.preorder_edge_iter():
            edge.length = edge.length * self.mutrate_per_site_per_generation
        return self.gene_tree

