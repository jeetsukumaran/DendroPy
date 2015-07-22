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
