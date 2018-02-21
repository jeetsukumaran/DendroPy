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
Scripts and other files that executed, sourced, invoked, or otherwise used by
various DendroPy entities.
"""

import os

def filepath(filename):
    try:
        import pkg_resources
        # note that this creates a temporary file with the contents of the
        # filename if the package is in an egg
        filepath = pkg_resources.resource_filename("dendropy", "utility/libexec/{}".format(filename))
        # print("-->{}".format(filepath))
    except:
        filepath = os.path.normpath(os.path.join(os.path.dirname(__file__), filename))
    return filepath

