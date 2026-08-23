#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
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

