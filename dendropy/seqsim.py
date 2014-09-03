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
DEPRECATED IN DENDROPY 4: USE `dendropy.simulate.charsim` instead.
"""

from dendropy.model.nucleotide import hky_char_matrix as generate_hky_char_matrix

from dendropy.utility import error
error.dendropy_module_migration_warning("dendropy.seqsim", "dendropy.simulate.charsim")

def hky_dataset(seq_len,
        tree_model,
        mutation_rate=1.0,
        kappa=1.0,
        base_freqs=[0.25, 0.25, 0.25, 0.25],
        root_states=None,
        dataset=None,
        rng=None):
    """
    Convenience class to wrap generation of a dataset based on
    the HKY model.
    `seq_len`       : length of sequence (number of characters)
    `tree_model`    : dendropy.Tree object
    `seq_model`    :  dendropy.model.discrete.DiscreteCharacterEvolutionModel object
    `mutation_rate` : mutation *modifier* rate (should be 1.0 if branch lengths
                      on tree reflect true expected number of changes
    `root_states`   : vector of root states (length must equal `seq_len`)
    `dataset`       : a dendropy.DataSet object.
                      if given, the new
                      dendropy.CharacterMatrix object will be added to
                      this (along with a new taxon_namespace if required). Otherwise,
                      a new dendropy.DataSet object will be created.
    `rng`           : random number generator; if not given, `GLOBAL_RNG` will be
                      used
    Returns: a dendropy.DataSet object object.
    """
    from dendropy.model import discrete
    from dendropy.model import nucleotide
    seq_model = nucleotide.Hky85CharacterEvolutionModel(kappa=kappa, base_freqs=base_freqs)
    return discrete.simulate_dataset(seq_len=seq_len,
        tree_model=tree_model,
        seq_model=seq_model,
        mutation_rate=mutation_rate,
        root_states=root_states,
        dataset=dataset,
        rng=rng)
