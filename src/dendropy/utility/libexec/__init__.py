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

