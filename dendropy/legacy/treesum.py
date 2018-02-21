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


