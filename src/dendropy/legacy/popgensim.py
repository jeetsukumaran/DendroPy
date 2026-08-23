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
DEPRECATED IN DENDROPY 4: USE `dendropy.simulate.popgensim` instead.
"""

from dendropy.utility import deprecate
from dendropy.simulate import popgensim

class FragmentedPopulations(popgensim.FragmentedPopulations):

    def __init__(self, *args, **kwargs):
        deprecate.dendropy_deprecation_warning(
                preamble="The 'dendropy.popgensim' module has moved to 'dendropy.simulate.popgensim'.",
                old_construct="from dendropy import popgensim\npopgensim.FragmentedPopulations(...)",
                new_construct="from dendropy.simulate import popgensim\npopgensim.FragmentedPopulations(...)",
                )
        popgensim.FragmentedPopulations.__init__(self, *args, **kwargs)
