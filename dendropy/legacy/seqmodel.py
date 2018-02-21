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
DEPRECATED IN DENDROPY 4: USE `dendropy.model.discrete`.
"""

from dendropy.model import discrete
from dendropy.utility import deprecate

class SeqModel(discrete.DiscreteCharacterEvolutionModel):
    def __init__(self, state_alphabet, rng=None):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: The 'dendropy.seqmodel.SeqModel' class has moved to 'dendropy.model.discrete.DiscreteCharacterEvolutionModel'.",
                old_construct="from dendropy import seqmodel\nm = seqmodel.SeqModel(...)",
                new_construct="from dendropy.model import discrete\nm = discrete.DiscreteCharacterEvolutionModel(...)")
        discrete.DiscreteCharacterEvolutionModel.__init__(
                self,
                state_alphabet=state_alphabet,
                rng=rng)

class Hky85SeqModel(discrete.Hky85):
    def __init__(self, kappa=1.0, base_freqs=None, state_alphabet=None, rng=None):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: The 'dendropy.seqmodel.Hky85SeqModel' class has moved to 'dendropy.model.discrete.Hky85'.",
                old_construct="from dendropy import seqmodel\nm = seqmodel.NucleotideSeqModel(...)",
                new_construct="from dendropy.model import discrete\ndiscrete.Hky85(...)")
        discrete.Hky85.__init__(
                self,
                kappa=kappa,
                base_freqs=base_freqs,
                state_alphabet=state_alphabet,
                rng=rng)

class Jc69SeqModel(discrete.Jc69):
    def __init__(self, state_alphabet=None, rng=None):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: The 'dendropy.seqmodel.Jc69SeqModel' class has moved to 'dendropy.model.discrete.Jc69'.",
                old_construct="from dendropy import seqmodel\nm = seqmodel.NucleotideSeqModel(...)",
                new_construct="from dendropy.model import discrete\ndiscrete.Jc69(...)")
        discrete.Jc69.__init__(
                self,
                state_alphabet=state_alphabet,
                rng=rng)

