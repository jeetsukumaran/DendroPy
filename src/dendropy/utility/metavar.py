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
