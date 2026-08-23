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
Wrappers for interacting with the ETE library.
"""

from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)
import dendropy

DENDROPY_ETE_INTEROPERABILITY = False
try:
    import ete2
    DENDROPY_ETE_INTEROPERABILITY = True
except ImportError:
    _LOG.warn("ete2 not installed: ETE interoperability not available")
else:

    def as_ete_object(o):
        if isinstance(o, ete2.Tree):
            return o
        elif isinstance(o, dendropy.Tree) or isinstance(o, dendropy.Node):
            s = o.as_newick_string() + ";"
#            _LOG.debug(s)
            return ete2.Tree(s)
        elif isinstance(o, list) or isinstance(o, dendropy.TreeList):
            return [as_ete_object(t) for t in o]
        else:
            raise ValueError("Object of type '%s' does not have a native ete2 representation" % type(o))

    def as_dendropy_object(o, taxon_set=None):
        if isinstance(o, dendropy.Tree) or isinstance(o, dendropy.Node) or isinstance(o, dendropy.TreeList):
            return o
        elif isinstance(o, ete2.Tree):
            s = o.write()
#            _LOG.debug(s)
            return dendropy.Tree.get_from_string(s, 'newick', taxon_set=taxon_set)
        elif isinstance(o, list) or isinstance(o, dendropy.TreeList):
            return dendropy.TreeList([as_dendropy_object(t, taxon_set=taxon_set) for t in o], taxon_set=taxon_set)
        else:
            raise ValueError("Object of type '%s' does not have a DendroPy representation" % type(o))

    def show(o):
        if not isinstance(o, dendropy.Tree) \
                and not isinstance(o, dendropy.Node)\
                and not isinstance(o, ete2.Tree):
            raise ValueError("Object of type '%s' cannot be rendered by ETE2" % type(o))
        ete_o = as_ete_object(o)
        ete_o.show()

