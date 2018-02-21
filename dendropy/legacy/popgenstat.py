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
DEPRECATED IN DENDROPY 4: USE `dendropy.calculate.popgenstat` instead.
"""

from dendropy.utility import deprecate
from dendropy.calculate import popgenstat

def num_segregating_sites(char_matrix, ignore_uncertain=True):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.popgenstat' module has moved to 'dendropy.calculate.popgenstat'.",
            old_construct="from dendropy import popgenstat\npopgenstat.num_segregating_sites(...)",
            new_construct="from dendropy.calculate import popgenstat\npopgenstat.num_segregating_sites(...)",
            )
    return popgenstat.num_segregating_sites(char_matrix, ignore_uncertain)

def average_number_of_pairwise_differences(char_matrix, ignore_uncertain=True):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.popgenstat' module has moved to 'dendropy.calculate.popgenstat'.",
            old_construct="from dendropy import popgenstat\npopgenstat.average_number_of_pairwise_differences(...)",
            new_construct="from dendropy.calculate import popgenstat\npopgenstat.average_number_of_pairwise_differences(...)",
            )
    return popgenstat.average_number_of_pairwise_differences(char_matrix=char_matrix, ignore_uncertain=ignore_uncertain)

def nucleotide_diversity(char_matrix, ignore_uncertain=True):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.popgenstat' module has moved to 'dendropy.calculate.popgenstat'.",
            old_construct="from dendropy import popgenstat\npopgenstat.nucleotide_diversity(...)",
            new_construct="from dendropy.calculate import popgenstat\npopgenstat.nucleotide_diversity(...)",
            )
    return popgenstat.nucleotide_diversity(char_matrix=char_matrix, ignore_uncertain=ignore_uncertain)

def tajimas_d(char_matrix, ignore_uncertain=True):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.popgenstat' module has moved to 'dendropy.calculate.popgenstat'.",
            old_construct="from dendropy import popgenstat\npopgenstat.tajimas_d(...)",
            new_construct="from dendropy.calculate import popgenstat\npopgenstat.tajimas_d(...)",
            )
    return popgenstat.tajimas_d(char_matrix=char_matrix, ignore_uncertain=ignore_uncertain)

def wattersons_theta(char_matrix, ignore_uncertain=True):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.popgenstat' module has moved to 'dendropy.calculate.popgenstat'.",
            old_construct="from dendropy import popgenstat\npopgenstat.watterson_theta(...)",
            new_construct="from dendropy.calculate import popgenstat\npopgenstat.watterson_theta(...)",
            )
    return popgenstat.wattersons_theta(char_matrix=char_matrix, ignore_uncertain=ignore_uncertain)

def derived_state_matrix(char_matrix, ancestral_seq=None):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.popgenstat' module has moved to 'dendropy.calculate.popgenstat'.",
            old_construct="from dendropy import popgenstat\npopgenstat.derived_state_matrix(...)",
            new_construct="from dendropy.calculate import popgenstat\npopgenstat.derived_state_matrix(...)",
            )
    return popgenstat.derived_state_matrix(
            char_matrix=char_matrix,
            ancestral_sequence=ancestral_seq)

def unfolded_site_frequency_spectrum(char_matrix, ancestral_seq=None, pad=True):
    deprecate.dendropy_deprecation_warning(
            preamble="The 'dendropy.popgenstat' module has moved to 'dendropy.calculate.popgenstat'.",
            old_construct="from dendropy import popgenstat\npopgenstat.unfolded_site_frequency_spectrum(...)",
            new_construct="from dendropy.calculate import popgenstat\npopgenstat.unfolded_site_frequency_spectrum(...)",
            )
    return popgenstat.unfolded_site_frequency_spectrum(
            char_matrix=char_matrix,
            ancestral_sequence=ancestral_seq,
            pad=pad)

class PopulationPairSummaryStatistics(popgenstat.PopulationPairSummaryStatistics):
    def __init__(self, *args, **kwargs):
        deprecate.dendropy_deprecation_warning(
                preamble="The 'dendropy.popgenstat' module has moved to 'dendropy.calculate.popgenstat'.",
                old_construct="from dendropy import popgenstat\npopgenstat.PopulationPairSummaryStatistics(...)",
                new_construct="from dendropy.calculate import popgenstat\npopgenstat.PopulationPairSummaryStatistics(...)",
                )
        popgenstat.PopulationPairSummaryStatistics.__init__(self, *args, **kwargs)
