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
DEPRECATED IN DENDROPY 4: USE `dendropy.model.discrete` or
`dendropy.model.nucleotide` instead.
"""

## legacy support
from dendropy.model.discrete import DiscreteCharacterEvolutionModel as SeqModel
from dendropy.model.discrete import DiscreteCharacterEvolver as SeqEvolver
from dendropy.model.discrete import simulate_discrete_char_dataset as generate_dataset
from dendropy.model.discrete import simulate_discrete_char_matrix as generate_char_matrix
from dendropy.model.nucleotide import Hky85CharacterEvolutionModel as Hky85SeqModel
from dendropy.model.nucleotide import Jc69CharacterEvolutionModel as Jc69CharacterEvolutionModel

from dendropy.utility import error
error.dendropy_module_migration_warning("dendropy.seqmodel", "dendropy.model.discrete")




