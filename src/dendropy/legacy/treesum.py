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
Split calculation and management.
DEPRECATED IN DENDROPY 4: USE `dendropy.calculate.treesum` instead.
"""

from dendropy.calculate import treesum
from dendropy.utility import deprecate

class TreeSummarizer(treesum.TreeSummarizer):
    def __init__(self, **kwargs):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: The 'dendropy.treesum.TreeSummarizer' class has moved to 'dendropy.calculate.treesum.TreeSummarizer'.",
                old_construct="from dendropy import treesum\nm = treesum.TreeSummarizer(...)",
                new_construct="from dendropy.calculate import treesum\nm = treesum.TreeSummarizer(...)")
        treesum.TreeSummarizer.__init__(self, **kwargs)

class TopologyCounter(treesum.TopologyCounter):
    def __init__(self, **kwargs):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: The 'dendropy.treesum.TopologyCounter' class has moved to 'dendropy.calculate.treesum.TopologyCounter'.",
                old_construct="from dendropy import treesum\nm = treesum.TopologyCounter(...)",
                new_construct="from dendropy.calculate import treesum\nm = treesum.TopologyCounter(...)")
        treesum.TopologyCounter.__init__(self, **kwargs)


