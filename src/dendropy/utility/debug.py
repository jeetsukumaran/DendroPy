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

