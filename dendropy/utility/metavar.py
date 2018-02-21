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
Environmental variables controlling DendroPy behavior (mostly for
development/testing usage).
"""

FAIL_INCOMPLETE_TESTS_ENVAR        = "DENDROPY_FAIL_INCOMPLETE_TESTS"
LOGGING_LEVEL_ENVAR                = "DENDROPY_LOGGING_LEVEL"
LOGGING_FORMAT_ENVAR               = "DENDROPY_LOGGING_FORMAT"
DENDROPY_PAUP_PATH_ENVAR           = "DENDROPY_PAUP_EXECUTABLE_PATH"
DENDROPY_RSCRIPT_PATH_ENVAR        = "DENDROPY_RSCRIPT_EXECUTABLE_PATH"

# error: Turn the warning into an exception.
# ignore: Discard the warning.
# always: Always emit a warning.
# default: Print the warning the first time it is generated from each location.
# module: Print the warning the first time it is generated from each module.
# once: Print the warning the first time it is generated.
DEPRECATION_WARNING_FILTER         = "DENDROPY_DEPRECATION_WARNINGS"
