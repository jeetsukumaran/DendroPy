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
This module provides a convenient interface that aggregates, wraps, and/or
implements functions and classes that simulate character data under various
models.
"""

from dendropy.model.continuous import evolve_continuous_char
from dendropy.model.discrete import DiscreteCharacterEvolutionModel
from dendropy.model.discrete import DiscreteCharacterEvolver
from dendropy.model.discrete import simulate_discrete_char_dataset
from dendropy.model.discrete import simulate_discrete_chars
from dendropy.model.discrete import Hky85
from dendropy.model.discrete import Jc69
from dendropy.model.discrete import hky85_chars
