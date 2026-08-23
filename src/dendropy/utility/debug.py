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
Various data structures.
"""

import inspect
import sys

def get_calling_code_info(stack_level):
    frame = inspect.stack()[stack_level]
    filename = inspect.getfile(frame[0])
    lineno = inspect.getlineno(frame[0])
    return filename, lineno

def dump_stack(out=None):
    if out is None:
        out = sys.stderr
    for frame, filename, line_num, func, source_code, source_index in inspect.stack()[2:]:
        if source_code is None:
            out.write("{}: {}\n".format(filename, line_num))
        else:
            out.write("{}: {}: {}\n".format(filename, line_num, source_code[source_index].strip()))

